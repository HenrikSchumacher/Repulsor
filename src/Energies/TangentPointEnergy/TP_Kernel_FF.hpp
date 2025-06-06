#pragma once

namespace Repulsor
{
    template<
        typename ClusterTree_T_,
        typename T1, typename T2, int q_flag,
        bool symmetricQ_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_, bool density_flag_
    >
    class TP_Kernel_FF : public FMM_Kernel_FF<
        ClusterTree_T_,symmetricQ_,
        energy_flag_,diff_flag_,metric_flag_, density_flag_
    >
    {
    private:
        
        using Base_T = FMM_Kernel_FF<
            ClusterTree_T_,symmetricQ_,
            energy_flag_,diff_flag_,metric_flag_, density_flag_
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
        using Base_T::T_DATA_DIM;
        using Base_T::S_DATA_DIM;
        
        static constexpr  Int ROWS      = 1 + AMB_DIM;
        static constexpr  Int COLS      = 1 + AMB_DIM;
        static constexpr LInt BLOCK_NNZ = 1 + 2 * AMB_DIM;
        static constexpr LInt DIAG_NNZ  = ROWS * COLS;
        
        static constexpr Real half      = 0.5;
        
        using Base_T::zero;
        using Base_T::one;
        using Base_T::two;
        using Base_T::symmetricQ;
        using Base_T::energy_flag;
        using Base_T::diff_flag;
        using Base_T::metric_flag;
        using Base_T::density_flag;
        
    public:
        
        TP_Kernel_FF() = delete;
        
        TP_Kernel_FF( mref<Configurator_T> conf, const Int thread_, const T1 q_half_, const T2 p_half_ )
        :   Base_T               (conf, thread_ )
        ,   q                    (two*q_half_   )
        ,   q_half               (q_half_       )
        ,   q_half_minus_1       (q_half-1      )
        ,   p                    (two*p_half_   )
        ,   p_half               (p_half_       )
        ,   minus_p_half         (-p_half       )
        ,   minus_p_half_minus_1 (-p_half-1     )
        {}
        
        TP_Kernel_FF( mref<TP_Kernel_FF> other, const Int thread_ )
        :   Base_T               (other, thread_                )
        ,   q                    (other.q                       )
        ,   q_half               (other.q_half                  )
        ,   q_half_minus_1       (other.q_half_minus_1          )
        ,   p                    (other.p                       )
        ,   p_half               (other.p_half                  )
        ,   minus_p_half         (other.minus_p_half            )
        ,   minus_p_half_minus_1 (other.minus_p_half_minus_1    )
        {}
        
        ~TP_Kernel_FF() = default;
        
#include "../../FMM/FMM_Kernel_Common.hpp"
#include "../../FMM/FMM_Kernel_FF_Common.hpp"
        
    protected:

        using Base_T::metric_data;
        
        using Base_T::S_diag;
        using Base_T::T_diag;
        
        using Base_T::DX;
        using Base_T::DY;
        
        using Base_T::a;
        using Base_T::x;
        using Base_T::P;
        
        using Base_T::b;
        using Base_T::y;
        using Base_T::Q;
        
        const Real q;
        const T1   q_half;
        const T1   q_half_minus_1;
        
        const Real p;
        const T2   p_half;
        const T2   minus_p_half;
        const T2   minus_p_half_minus_1;
        
        Tiny::Vector<BLOCK_NNZ,Real,Int> ij_block;
        Tiny::Matrix<ROWS,COLS,Real,Int> ii_block;
        Tiny::Matrix<ROWS,COLS,Real,Int> jj_block;
        
    public:
        
        TOOLS_FORCE_INLINE Real Compute( const LInt k_global )
        {
            TOOLS_MAKE_FP_FAST()
            
            Real v    [AMB_DIM ] = {};
            Real Pv   [AMB_DIM ] = {};
            Real Qv   [AMB_DIM ] = {};
            Real dEdv [AMB_DIM ] = {};
            Real V    [PROJ_DIM] = {};
            
            Real r2        = zero;
            Real rCosPhi_2 = zero;
            Real rCosPsi_2 = zero;
            
            for( Int l = 0; l < AMB_DIM; ++l )
            {
                v[l] = y[l] - x[l];
                r2  += v[l] * v[l];
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
                        V[k] = ( one + static_cast<Real>(i!=j) )*v[i]*v[j];
                    }
                }
                rCosPhi_2 += v[i] * Pv[i];
                rCosPsi_2 += v[i] * Qv[i];
            }
            
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
                // Does not make any sense.
                rCosPhi_q_minus_2 = one;
                // Does not make any sense.
                rCosPsi_q_minus_2 = one;
                
                // |P*(y-x)|^q
                rCosPhi_q = Power<Real,T1>( Abs(rCosPhi_2), q_half );
                // |Q*(y-x)|^q
                rCosPsi_q = Power<Real,T1>( Abs(rCosPsi_2), q_half );
            }
            
            const Real Num = ( rCosPhi_q + rCosPsi_q );

            // E = ( |P*(y-x)|^q + |Q*(y-x)|^q) / |y-x|^p
            const Real E = Num * r_minus_p;

            if constexpr ( density_flag )
            {
                DX[0] += b * rCosPhi_q * r_minus_p;
                DY[0] += a * rCosPsi_q * r_minus_p;
            }
            
            if constexpr ( diff_flag || metric_flag )
            {
                // factor = q / |y-x|^p
                const Real factor = ( q_flag == 0 ? zero : q * r_minus_p );
                // K_xy = q * |P*(y-x)|^(q-2) / |y-x|^p
                const Real K_xy   = ( q_flag == 0 ? zero : factor * rCosPhi_q_minus_2 );
                // K_yx = q * |Q*(y-x)|^(q-2) / |y-x|^p
                const Real K_yx   = ( q_flag == 0 ? zero : factor * rCosPsi_q_minus_2 );
                // H    = p * ( |P*(y-x)|^q + |Q*(y-x)|^q) / |y-x|^(p+2)
                const Real H = p * r_minus_p_minus_2 * Num;
                
                if constexpr ( diff_flag )
                {
                    Real dEdvx = zero;
                    Real dEdvy = zero;
                    
                    for( Int i = 0; i < AMB_DIM; ++i )
                    {
                        dEdv[i] = ( q_flag == 0 ? zero : K_xy * Pv[i] + K_yx * Qv[i]) - H * v[i];
                        
                        dEdvx   += dEdv[i] * x[i];
                        dEdvy   += dEdv[i] * y[i];
                        
                        DX[1+i] -= b * dEdv[i];
                        DY[1+i] += a * dEdv[i];
                    }
                    
                    DX[0] += b * ( E - half * factor * rCosPhi_q + dEdvx );
                    DY[0] += a * ( E - half * factor * rCosPsi_q - dEdvy );
                    
                    const Real half_b_K_xy = half * b * K_xy;
                    const Real half_a_K_yx = half * a * K_yx;
                    
                    for( Int k = 0; k < PROJ_DIM; ++k )
                    {
                        DX[1+AMB_DIM+k] += half_b_K_xy * V[k];
                        DY[1+AMB_DIM+k] += half_a_K_yx * V[k];
                    }
                }
                
                if constexpr ( metric_flag )
                {
/*
// ij_block
//    /                                                                 \
//    |       - C_xy       K_yx * v[0]    K_yx * v[1]    K_yx * v[2]    |
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
//    |         C_xy          K_xy * v[0]         K_xy * v[1]         K_xy * v[2]       |
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
//    |        C_xy         - K_yx * v[0]       - K_yx * v[1]       - K_yx * v[2]       |
//    |                                                                                 |
//    |   - K_yx * v[0]   K_yx * v[0] * v[0]  K_yx * v[0] * v[1]  K_yx * v[0] * v[2]    |
//    |                                                                                 |
//    |   - K_yx * v[1]   K_yx * v[1] * v[0]  K_yx * v[1] * v[1]  K_yx * v[1] * v[2]    |
//    |                                                                                 |
//    |   - K_yx * v[2]   K_yx * v[2] * v[0]  K_yx * v[2] * v[1]  K_yx * v[2] * v[2]    |
//    \                                                                                 /
*/
                    const Real C_xy = K_xy + K_yx + H;
                    
                    const Real b_over_a   = b/a;
                    const Real a_over_b   = a/b;
                    
                    const Real b_over_a_K_xy = b_over_a * K_xy;
                    const Real a_over_b_K_yx = a_over_b * K_yx;
                    
                    ij_block[0]     = - C_xy;
                    ii_block[0][0] +=   b_over_a * C_xy;
                    jj_block[0][0]  =   a_over_b * C_xy;
                    
                    for( Int i = 1; i < COLS; ++i )
                    {
                        ij_block[i]         =   K_yx * v[i-1];
                        ij_block[AMB_DIM+i] = - K_xy * v[i-1];
                        
                        // store only upper triangle
                        ii_block[0][i] +=   b_over_a_K_xy * v[i-1];
                        jj_block[0][i]  = - a_over_b_K_yx * v[i-1];
                        
                        // store also lower triangle
                        ii_block[i][0] +=   b_over_a_K_xy * v[i-1];
                        jj_block[i][0]  = - a_over_b_K_yx * v[i-1];
                    }
                    
                    for( Int i = 1; i < ROWS; ++i )
                    {
                        {
                            const Real vv = v[i-1] * v[i-1];
                            ii_block[i][i] += b_over_a_K_xy * vv;
                            jj_block[i][i]  = a_over_b_K_yx * vv;
                        }
                        for( Int j = i+1; j < COLS; ++j )
                        {
                            const Real vv = v[i-1] * v[j-1];
                            
                            // store only upper triangle
                            ii_block[i][j] += b_over_a_K_xy * vv;
                            jj_block[i][j]  = a_over_b_K_yx * vv;
                            
                            // store also lower triangle
                            ii_block[j][i] += b_over_a_K_xy * vv;
                            jj_block[j][i]  = a_over_b_K_yx * vv;
                        }
                    }
                    
                    ij_block.Write( &metric_data[BLOCK_NNZ * k_global] );
                }
            }

            if constexpr ( energy_flag )
            {
                return  this->symmetry_factor * a * E * b;
            }
            else
            {
                return zero;
            }
        }
        
        
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
                // We can do an overwrite here.
//                ij_block.SetZero();
                
                // We can do an overwrite here.
//                jj_block.SetZero();
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
            return std::string("TP_Kernel_FF<")
            + "...,"
            + TypeName<T1> + ","
            + TypeName<T2> + ","
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(metric_flag) +  ","
            + ToString(density_flag) +
            ">" + this->ThreadString();
        }
    };

} // namespace Repulsor
