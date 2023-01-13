#pragma once

#define BASE FMM_Kernel_VF<                                             \
        S_DOM_DIM_,T_DOM_DIM_,ClusterTree_T_,is_symmetric_,             \
        energy_flag_,diff_flag_,metric_flag_                            \
    >

namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename ClusterTree_T_,
        typename T1, typename T2,
        bool is_symmetric_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_
    >
    class TP_Kernel_VF : public BASE
    {
    private:
        
        using Base_T = BASE;
        
    public:
        
        using ClusterTree_T      = ClusterTree_T_;
        
        using Real               = typename Base_T::Real;
        using SReal              = typename Base_T::SReal;
        using ExtReal            = typename Base_T::ExtReal;
        using Int                = typename Base_T::Int;
        using LInt               = typename Base_T::LInt;
        
        using Configurator_T     = typename Base_T::Configurator_T;
        using Values_T           = typename Base_T::Values_T;
        using ValueContainer_T   = typename Base_T::ValueContainer_T;
        
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
        
        static constexpr  Int ROWS      = 1 + AMB_DIM;
        static constexpr  Int COLS      = 1 + AMB_DIM;
        static constexpr LInt BLOCK_NNZ = 1 + 2 * AMB_DIM;
        static constexpr LInt DIAG_NNZ  = ROWS * COLS;
        
        using Base_T::zero;
        using Base_T::one;
        using Base_T::two;
        using Base_T::is_symmetric;
        
    public:
        
        TP_Kernel_VF() = delete;
        
        TP_Kernel_VF(
            Configurator_T & conf,
            const Real theta_, const Int max_refinement_, const T1 q_half_, const T2 p_half_
        )
        :   Base_T               (conf, theta_, max_refinement_ )
        ,   q                    (two*q_half_)
        ,   q_half               (q_half_    )
        ,   q_half_minus_1       (q_half-1   )
        ,   p                    (two*p_half_)
        ,   p_half               (p_half_    )
        ,   minus_p_half         (-p_half    )
        ,   minus_p_half_minus_1 (-p_half-1  )
        {}
        
        TP_Kernel_VF( TP_Kernel_VF & other )
        :   Base_T               (other                     )
        ,   q                    (other.q                   )
        ,   q_half               (other.q_half              )
        ,   q_half_minus_1       (other.q_half_minus_1      )
        ,   p                    (other.p                   )
        ,   p_half               (other.p_half              )
        ,   minus_p_half         (other.minus_p_half        )
        ,   minus_p_half_minus_1 (other.minus_p_half_minus_1)
        {}
        
        ~TP_Kernel_VF() = default;
        
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
        
        using Base_T::lin_k;
        
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
        
        
        Real ij_block [BLOCK_NNZ]  = {};
        Real ii_block [ROWS][COLS] = {{}};
        Real jj_block [ROWS][COLS] = {{}};
        
    protected:
        
        force_inline Real compute()
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
                    const Int k = lin_k[i][j];
                    Pv[i] += P[k] * v[j];
                    Qv[i] += Q[k] * v[j];
                    if( j >= i )
                    {
                        V[k] = ( static_cast<Real>(1) + static_cast<Real>(i!=j) )*v[i]*v[j];
                    }
                }
                rCosPhi_2 += v[i] * Pv[i];
                rCosPsi_2 += v[i] * Qv[i];
            }
            
            // TODO: It probably suffices to multiply w only once against r_minus_p_minus_2, because the latter is multiplied to almost everything.
            
            // |P*(y-x)|^{q-2}
            const Real rCosPhi_q_minus_2 = MyMath::pow<Real,T1>( fabs(rCosPhi_2), q_half_minus_1);
            // |Q*(y-x)|^{q-2}
            const Real rCosPsi_q_minus_2 = MyMath::pow<Real,T1>( fabs(rCosPsi_2), q_half_minus_1);
            // r^{2-p}
            const Real r_minus_p_minus_2 = MyMath::pow<Real,T2>( r2, minus_p_half_minus_1 );
            // |y-x|^-p
            const Real r_minus_p = r_minus_p_minus_2 * r2;
            // |P*(y-x)|^q
            const Real rCosPhi_q = rCosPhi_q_minus_2 * rCosPhi_2;
            // |Q*(y-x)|^q
            const Real rCosPsi_q = rCosPsi_q_minus_2 * rCosPsi_2;
            
            const Real Num = ( rCosPhi_q + rCosPsi_q );

            const Real E = Num * r_minus_p;

            if constexpr ( diff_flag || metric_flag )
            {
                // Needed for both the differential and the metric.
                
                const Real factor = q_half * r_minus_p;         // = q / |y-x|^p
                
                const Real K_xy = factor * rCosPhi_q_minus_2;   // = |P*(y-x)|^(q-2) / |y-x|^p
                const Real K_yx = factor * rCosPsi_q_minus_2;   // = |Q*(y-x)|^(q-2) / |y-x|^p

                if constexpr ( diff_flag )
                {
                    const Real H = - p * r_minus_p_minus_2 * Num;
                    
                    Real dEdvx = zero;
                    Real dEdvy = zero;
                    
                    const Real wa = w * a;
                    const Real wb = w * b;
                    
                    for( Int i = 0; i < AMB_DIM; ++i )
                    {
                        dEdv[i] = two * ( K_xy * Pv[i] + K_yx * Qv[i] ) + H * v[i];
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
                    
                    DX[0] += wb * ( E - factor * rCosPhi_q + dEdvx );
                    DY[0] += wa * ( E - factor * rCosPsi_q - dEdvy );
                    
                    const Real w_b_K_xy = wb * K_xy;
                    const Real w_a_K_yx = wa * K_yx;
                    
                    for( Int k = 0; k < PROJ_DIM; ++k )
                    {
                        DX[1+S_COORD_DIM+k] += w_b_K_xy * V[k];
                        DY[1+T_COORD_DIM+k] += w_a_K_yx * V[k];
                    }
                }
                
                if constexpr ( metric_flag )
                {
                    const Real w_K_xy = w * K_xy;
                    const Real w_K_yx = w * K_yx;
                    
/*
// ij_block
//    /                                                                 \
//    |   - K_xy - K_yx    K_yx * v[0]    K_yx * v[1]    K_yx * v[2]    |
//    |                                                                 |
//    |   - K_xy * v[0]         0              0              0         |
//    |                                                                 |
//    |   - K_xy * v[1]         0              0              0         |
//    |                                                                 |
//    |   - K_xy * v[2]         0              0              0         |
//    \                                                                 /
//
//    This are 1 + 2 * AMB_DIM nonzero values.
//    It is tempting to compress this into 2 + AMB_DIM values.
//    But it did not appear to be faster...
//
//
// ii block
//    /                                                                                 \
//    |     K_xy + K_yx      K_xy * v[0]         K_xy * v[1]         K_xy * v[2]        |
//    |                                                                                 |
//    |     K_xy * v[0]   K_xy * v[0] * v[0]  K_xy * v[0] * v[1]  K_xy * v[0] * v[2]    |
//    |                                                                                 |
//    |     K_xy * v[1]   K_xy * v[1] * v[0]  K_xy * v[1] * v[1]  K_xy * v[1] * v[2]    |
//    |                                                                                 |
//    |     K_xy * v[2]   K_xy * v[2] * v[0]  K_xy * v[2] * v[1]  K_xy * v[2] * v[2]    |
//    \                                                                                 /
//
// jj block
//
//    /                                                                                 \
//    |     K_xy + K_yx     - K_yx * v[0]       - K_yx * v[1]       - K_yx * v[2]       |
//    |                                                                                 |
//    |   - K_yx * v[0]   K_yx * v[0] * v[0]  K_yx * v[0] * v[1]  K_yx * v[0] * v[2]    |
//    |                                                                                 |
//    |   - K_yx * v[1]   K_yx * v[1] * v[0]  K_yx * v[1] * v[1]  K_yx * v[1] * v[2]    |
//    |                                                                                 |
//    |   - K_yx * v[2]   K_yx * v[2] * v[0]  K_yx * v[2] * v[1]  K_yx * v[2] * v[2]    |
//    \                                                                                 /
//
*/
                    
                    const Real w_K_sym = w_K_xy + w_K_yx;
                    
                    const Real b_over_a = b/a;
                    const Real a_over_b = a/b;
                    
                    const Real b_over_a_w_K_xy = b_over_a * w_K_xy;
                    const Real a_over_b_w_K_yx = a_over_b * w_K_yx;
                    
                    ij_block[0]    -=   w_K_sym;
                    ii_block[0][0] +=   b_over_a * w_K_sym;
                    jj_block[0][0] +=   a_over_b * w_K_sym;
                    
                    for( Int i = 1; i < COLS; ++i )
                    {
                        ij_block[i]         +=   w_K_yx * v[i-1];
                        ij_block[AMB_DIM+i] -=   w_K_xy * v[i-1];
                        
                        // store only upper triangle
                        ii_block[0][i] +=   b_over_a_w_K_xy * v[i-1];
                        jj_block[0][i] -=   a_over_b_w_K_yx * v[i-1];
                        
                        // store also lower triangle
                        ii_block[i][0] +=   b_over_a_w_K_xy * v[i-1];
                        jj_block[i][0] -=   a_over_b_w_K_yx * v[i-1];
                    }
                    
                    for( Int i = 1; i < ROWS; ++i )
                    {
                        {
                            const Real vv = v[i-1] * v[i-1];
                            ii_block[i][i] += b_over_a_w_K_xy * vv;
                            jj_block[i][i] += a_over_b_w_K_yx * vv;
                        }
                        for( Int j = i+1; j < COLS; ++j )
                        {
                            const Real vv = v[i-1] * v[j-1];
                            
                            // store only upper triangle
                            ii_block[i][j] += b_over_a_w_K_xy * vv;
                            jj_block[i][j] += a_over_b_w_K_yx * vv;
                            
                            // store also lower triangle
                            ii_block[j][i] += b_over_a_w_K_xy * vv;
                            jj_block[j][i] += a_over_b_w_K_yx * vv;
                        }
                    }
                }
                
            }
            
            if constexpr ( energy_flag )
            {
                return  a * E * b;
            }
            else
            {
                return zero;
            }
        }
        
        
#include "../../FMM/FMM_Kernel_Common.hpp"
        
// Now load the actual Compute method.
#include "../../FMM/FMM_Kernel_VF_Common.hpp"
        
    public:
        
        force_inline void LoadS( const Int i_global )
        {
            this->loadS( i_global );
            
            if constexpr ( metric_flag )
            {
                zerofy_buffer<DIAG_NNZ>( &ii_block[0][0] );
            }
        }
        
        force_inline void WriteS( const Int i_global )
        {
            this->writeS( i_global );
            
            if constexpr ( metric_flag )
            {
                add_to_buffer<DIAG_NNZ>( &ii_block[0][0], &S_diag[DIAG_NNZ * i_global] );
            }
        }
        
        force_inline void LoadT( const Int j_global )
        {
            this->loadT( j_global );
            
            if constexpr ( metric_flag )
            {
                zerofy_buffer<BLOCK_NNZ>( &ij_block[0] );
                
                zerofy_buffer<DIAG_NNZ>( &jj_block[0][0] );
            }
        }

        
        force_inline void WriteT( const Int j_global )
        {
            this->writeT( j_global );
            
            if constexpr ( metric_flag )
            {
                add_to_buffer<DIAG_NNZ>( &jj_block[0][0], &T_diag[DIAG_NNZ * j_global] );
            }
        }
        
        std::string ClassName() const
        {
            return "TP_Kernel_VF<"
            + ToString(S_DOM_DIM) + ","
            + ToString(T_DOM_DIM) + ","
            + this->S.ClassName() + ","
            + TypeName<T1>::Get() + ","
            + TypeName<T2>::Get() + ","
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(metric_flag) +
            ">";
        }
        
    }; // class TP_Kernel_VF

} // namespace Repulsor

#undef BASE
