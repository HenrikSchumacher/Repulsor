#pragma once

namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename ClusterTree_T_,
        typename T1, typename T2, int q_flag,
        bool symmetricQ_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_, bool density_flag_
    >
    class TP0_Kernel_VF : public FMM_Kernel_VF<
        S_DOM_DIM_,T_DOM_DIM_,ClusterTree_T_,symmetricQ_,
        energy_flag_,diff_flag_,metric_flag_,density_flag_
    >
    {
    private:
        
        using Base_T = FMM_Kernel_VF<
            S_DOM_DIM_,T_DOM_DIM_,ClusterTree_T_,symmetricQ_,
            energy_flag_,diff_flag_,metric_flag_,density_flag_
        >;
        
    public:
        
        using ClusterTree_T      = ClusterTree_T_;
        
        using Real               = typename Base_T::Real;
        using SReal              = typename Base_T::SReal;
        using ExtReal            = typename Base_T::ExtReal;
        using Int                = typename Base_T::Int;
        using LInt               = typename Base_T::LInt;
        
        using Configurator_T     = typename Base_T::Configurator_T;
        using ValueContainer_T   = typename Base_T::ValueContainer_T;
        using Values_T           = typename ValueContainer_T::Values_T;

        using Base_T::AMB_DIM;
        using Base_T::PROJ_DIM;
        using Base_T::S_DOM_DIM;
        using Base_T::T_DOM_DIM;
        using Base_T::S_COORD_DIM;
        using Base_T::T_COORD_DIM;
        using Base_T::T_DATA_DIM;
        using Base_T::S_DATA_DIM;
        using Base_T::energy_flag;
        using Base_T::diff_flag;
        using Base_T::metric_flag;
        using Base_T::density_flag;
        
        static constexpr  Int ROWS      = 1 + AMB_DIM;
        static constexpr  Int COLS      = 1 + AMB_DIM;
        static constexpr  Int BLOCK_NNZ = 1 + 2 * AMB_DIM;
        static constexpr  Int DIAG_NNZ  = ROWS * COLS;
        
        static constexpr Real half      = Scalar::Half<Real>;
        
        using Base_T::S;
        using Base_T::T;
        
        using Base_T::zero;
        using Base_T::one;
        using Base_T::two;
        using Base_T::symmetricQ;
        
    public:
        
        TP0_Kernel_VF(
            mref<Configurator_T> conf, const Int thread_,
            const Real theta_, const Int max_refinement_, const T1 q_half_, const T2 p_half_
        )
        :   Base_T               (conf, thread_, theta_, max_refinement_ )
        ,   q                    (two*q_half_)
        ,   q_half               (q_half_    )
        ,   q_half_minus_1       (q_half-1   )
        ,   p                    (two*p_half_)
        ,   p_half               (p_half_    )
        ,   minus_p_half         (-p_half    )
        ,   minus_p_half_minus_1 (-p_half-1  )
        ,   s                    ((p - S_DOM_DIM) / q)
        ,   s_exp                ( -Scalar::Half<Real> * (Scalar::Two<Real> * (s - Scalar::One<Real>) + S_DOM_DIM))
        {}
        
        TP0_Kernel_VF( mref<TP0_Kernel_VF> other, const Int thread_ )
        :   Base_T               (other, thread_            )
        ,   q                    (other.q                   )
        ,   q_half               (other.q_half              )
        ,   q_half_minus_1       (other.q_half_minus_1      )
        ,   p                    (other.p                   )
        ,   p_half               (other.p_half              )
        ,   minus_p_half         (other.minus_p_half        )
        ,   minus_p_half_minus_1 (other.minus_p_half_minus_1)
        ,   s                    ((p - S_DOM_DIM) / q)
        ,   s_exp                ( -Scalar::Half<Real> * (Scalar::Two<Real> * (s - Scalar::One<Real>) + S_DOM_DIM))
        {}
        

        // Default constructor
        TP0_Kernel_VF() = delete;
        // Destructor
        virtual ~TP0_Kernel_VF() override = default;
        // Copy constructor
        TP0_Kernel_VF( const TP0_Kernel_VF & other ) = default;
        // Copy assignment operator
        TP0_Kernel_VF & operator=( const TP0_Kernel_VF & other ) = default;
        // Move constructor
        TP0_Kernel_VF( TP0_Kernel_VF && other ) = default;
        // Move assignment operator
        TP0_Kernel_VF & operator=( TP0_Kernel_VF && other ) = default;
        
    protected:
        
        using Base_T::metric_data;
        
        using Base_T::S_diag;
        using Base_T::T_diag;
        
        using Base_T::DX;
        using Base_T::DY;
        
        using Base_T::a;
        using Base_T::x;
        using Base_T::x_buffer;
        using Base_T::P;
        
        using Base_T::b;
        using Base_T::y;
        using Base_T::y_buffer;
        using Base_T::Q;
        
        using Base_T::S_Tree;
        using Base_T::T_Tree;
        
        using Base_T::lambda;
        using Base_T::mu;
        
        const Real q;
        const T1   q_half;
        const T1   q_half_minus_1;
        
        const Real p;
        const T2   p_half;
        const T2   minus_p_half;
        const T2   minus_p_half_minus_1;
        
        const Real s;
        const Real s_exp;
        
        Tiny::Vector<BLOCK_NNZ,Real,Int> ij_block;
        Tiny::Vector< DIAG_NNZ,Real,Int> ii_block;
        Tiny::Vector< DIAG_NNZ,Real,Int> jj_block;
        
#include "../../FMM/FMM_Kernel_Common.hpp"
#include "../../FMM/FMM_Kernel_VF_Common.hpp"
        
        
    protected:
        
        TOOLS_FORCE_INLINE Real compute()
        {
            Real x    [AMB_DIM]  = {};
            Real y    [AMB_DIM]  = {};
            Real v    [AMB_DIM]  = {};
            Real Pv   [AMB_DIM]  = {};
            Real Qv   [AMB_DIM]  = {};
            Real dEdv [AMB_DIM]  = {};
            Real V    [PROJ_DIM] = {};
            
            Real r2        = zero;
            Real rCosPhi_2 = zero;
            Real rCosPsi_2 = zero;
            
            const Real w = S_Tree.Weight() * T_Tree.Weight();
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                x[i] = lambda[0] * x_buffer[AMB_DIM*0 +i];
                y[i] = mu    [0] * y_buffer[AMB_DIM*0 +i];

                for( Int ii = 1; ii < S_DOM_DIM+1; ++ii )
                {
                    x[i] += lambda[ii] * x_buffer[AMB_DIM*ii +i];
                }
                
                for( Int ii = 1; ii < T_DOM_DIM+1; ++ii )
                {
                    y[i] += mu    [ii] * y_buffer[AMB_DIM*ii +i];
                }
                
                v[i] = y[i] - x[i];
                r2  += v[i] * v[i];
            }
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    const Int k = lin_k<AMB_DIM>(i,j);
                    Pv[i] += P[k] * v[j];
                    Qv[i] += Q[k] * v[j];
                    if( j >= i )
                    {
                        V[k] = ( Scalar::One<Real> + static_cast<Real>(i!=j) )*v[i]*v[j];
                    }
                }
                rCosPhi_2 += v[i] * Pv[i];
                rCosPsi_2 += v[i] * Qv[i];
            }
            
            // TODO: It probably suffices to multiply w only once against r_minus_p_minus_2, because the latter is multiplied to almost everything.
            
            
            Real result = 0;
            
            if constexpr ( energy_flag || density_flag || diff_flag )
            {
                // r^{-p-2}
                const Real r_minus_p_minus_2 = Power<Real,T2>( r2, minus_p_half_minus_1 );
                // |y-x|^-p
                const Real r_minus_p = r_minus_p_minus_2 * r2;
                
                Real rCosPhi_q_minus_2;
                Real rCosPsi_q_minus_2;
                Real rCosPhi_q;
                Real rCosPsi_q;
                
                if constexpr ( q_flag >= 2 )
                {
                    // |P*(y-x)|^{q-2}
                    rCosPhi_q_minus_2 = Power<Real,T1>( Abs(rCosPhi_2), q_half_minus_1 );
                    
                    // |Q*(y-x)|^{q-2}
                    rCosPsi_q_minus_2 = Power<Real,T1>( Abs(rCosPsi_2), q_half_minus_1 );
                    
                    // |P*(y-x)|^q
                    rCosPhi_q = rCosPhi_q_minus_2 * rCosPhi_2;
                    // |Q*(y-x)|^q
                    rCosPsi_q = rCosPsi_q_minus_2 * rCosPsi_2;
                }
                else
                {
                    // Does not make any sense unless q == 2.
                    rCosPhi_q_minus_2 = one;
                    // Does not make any sense unless q == 2.
                    rCosPsi_q_minus_2 = one;
                    
                    // |P*(y-x)|^q
                    rCosPhi_q = Power<Real,T1>( Abs(rCosPhi_2), q_half );
                    // |Q*(y-x)|^q
                    rCosPsi_q = Power<Real,T1>( Abs(rCosPsi_2), q_half );
                }
                
                const Real Num = ( rCosPhi_q + rCosPsi_q );

                // E = ( |P*(y-x)|^q + |Q*(y-x)|^q) / |y-x|^p
                const Real E = Num * r_minus_p;
                
                if constexpr ( energy_flag )
                {
                    // TODO: Is there a reason why this->symmetry_factor does not show up?
                    result = a * E * b;
                }
                
                if constexpr ( density_flag )
                {
                    DX[0] += b * rCosPhi_q * r_minus_p;
                    DY[0] += a * rCosPsi_q * r_minus_p;
                }
                
                if constexpr ( diff_flag )
                {
                    // factor = q / |y-x|^p
                    const Real factor = ( q_flag == 0 ? zero : q * r_minus_p );
                    // K_xy = q * |P*(y-x)|^(q-2) / |y-x|^p
                    const Real K_xy   = ( q_flag == 0 ? zero : factor * rCosPhi_q_minus_2 );
                    // K_yx = q * |Q*(y-x)|^(q-2) / |y-x|^p
                    const Real K_yx   = ( q_flag == 0 ? zero : factor * rCosPsi_q_minus_2 );
                    // H    = p * ( |P*(y-x)|^q + |Q*(y-x)|^q) / |y-x|^(p+2)
                    const Real H = p * r_minus_p_minus_2 * Num;
                    
                    Real dEdvx = zero;
                    Real dEdvy = zero;
                    
                    const Real wa = w * a;
                    const Real wb = w * b;
                    
                    for( Int i = 0; i < AMB_DIM; ++i )
                    {
                        dEdv[i] = ( q_flag == 0 ? zero : K_xy * Pv[i] + K_yx * Qv[i]) - H * v[i];
                        dEdvx  += dEdv[i] * x[i];
                        dEdvy  += dEdv[i] * y[i];
                        
                        for( Int ii = 0; ii < S_DOM_DIM+1; ++ii )
                        {
                            DX[1+AMB_DIM*ii+i] -= lambda[ii] * wb * dEdv[i];
                        }
                        
                        for( Int ii = 0; ii < T_DOM_DIM+1; ++ii )
                        {
                            DY[1+AMB_DIM*ii+i] +=     mu[ii] * wa * dEdv[i];
                        }
                    }
                    
                    DX[0] += wb * ( E - half * factor * rCosPhi_q + dEdvx );
                    DY[0] += wa * ( E - half * factor * rCosPsi_q - dEdvy );
                    
                    const Real half_w_b_K_xy = half * wb * K_xy;
                    const Real half_w_a_K_yx = half * wa * K_yx;
                    
                    for( Int k = 0; k < PROJ_DIM; ++k )
                    {
                        DX[1+S_COORD_DIM+k] += half_w_b_K_xy * V[k];
                        DY[1+T_COORD_DIM+k] += half_w_a_K_yx * V[k];
                    }
                }
            }
            
            if constexpr ( metric_flag )
            {
/*
// ij_block
//              /                               \
//              |  a[0]     0       0       0   |
//              |                               |
//              |   0      a[1]     0       0   |
//              |                               |
//              |   0       0      a[1]     0   |
//              |                               |
//              |   0       0       0      a[1] |
//              \                               /
//
// ii_block
//              /                               \
//              | -a[0]     0       0       0   |
//              |                               |
//              |   0     -a[1]     0       0   |
//              |                               |
//              |   0       0     -a[1]     0   |
//              |                               |
//              |   0       0       0     -a[1] |
//              \                               /
//
// jj_block
//              /                               \
//              | -a[0]     0       0       0   |
//              |                               |
//              |   0     -a[1]     0       0   |
//              |                               |
//              |   0       0     -a[1]     0   |
//              |                               |
//              |   0       0       0     -a[1] |
//              \                               /
 */
                
                const Real r4 = r2 * r2;

                // The following line makes up approx 2/3 of this function's runtime! This is why we avoid pow as much as possible and replace it with Power.;
                // I got it down to this single call to pow. We might want to generate a lookup table for it...;
                // The factor of (-2.) is here, because we assemble the _metric_, not the kernel.;
                const Real a_1 = w * - Scalar::Two<Real> * Power(r2, s_exp);
                
                const Real a_0 = w * Scalar::Half<Real> * (rCosPhi_2 + rCosPsi_2) / r4 * a_1;
                
                const Real b_over_a   = b/a;
                const Real a_over_b   = a/b;
                
                ij_block[0] +=   a_0;
                ii_block[0] -=   b_over_a * a_0;
                jj_block[0] -=   a_over_b * a_0;
                
                ij_block[1] +=   a_1;
                ii_block[1] -=   b_over_a * a_1;
                jj_block[1] -=   a_over_b * a_1;
                
            }
            
            return result;
        }
        
        
    public:
        
        TOOLS_FORCE_INLINE void LoadS( const Int i_global )
        {
            this->loadS( i_global );
            
            if constexpr ( metric_flag )
            {
                ii_block.SetZero();
            }
        }
        
        TOOLS_FORCE_INLINE void WriteS( const Int i_global )
        {
            this->writeS( i_global );
            
            if constexpr ( metric_flag )
            {
                ii_block.AddTo( &S_diag[DIAG_NNZ * i_global] );
            }
        }
        
        TOOLS_FORCE_INLINE void LoadT( const Int j_global )
        {
            this->loadT( j_global );
            
            if constexpr ( metric_flag )
            {
                ij_block.SetZero();
                
                jj_block.SetZero();
            }
        }

        TOOLS_FORCE_INLINE void WriteT( const Int j_global )
        {
            this->writeT( j_global );
            
            if constexpr ( metric_flag )
            {
                jj_block.AddTo( &T_diag[DIAG_NNZ * j_global] );
            }
        }
        
    public:
        
        std::string ClassName() const
        {
            return "TP0_Kernel_VF<"
            + ToString(S_DOM_DIM) + ","
            + ToString(T_DOM_DIM) + ","
//            + S.ClassName() + ","
            + "...,"
            + TypeName<T1> + ","
            + TypeName<T2> + ","
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(metric_flag) + ","
            + ToString(density_flag) +
            ">";
        }
    };

} // namespace Repulsor
