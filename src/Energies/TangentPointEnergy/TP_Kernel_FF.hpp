#pragma once

#define BASE FMM_Kernel_FF<ClusterTree_T_,is_symmetric_,energy_flag_,diff_flag_,metric_flag_>

namespace Repulsor
{
    template<
        typename ClusterTree_T_,
        typename T1, typename T2,
        bool is_symmetric_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_
    >
    class TP_Kernel_FF : public BASE
    {
    private:
        
        using Base_T = BASE;
        
    public:
        
        using ClusterTree_T      = ClusterTree_T_;
        using Values_T           = typename Base_T::Values_T;
        using ValueContainer_T   = typename Base_T::ValueContainer_T;
        
        using Real               = typename Base_T::Real;
        using SReal              = typename Base_T::SReal;
        using ExtReal            = typename Base_T::ExtReal;
        using Int                = typename Base_T::Int;
        using LInt               = typename Base_T::LInt;
        
        using Configurator_T     = typename Base_T::Configurator_T;
        
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
        using Base_T::is_symmetric;
        using Base_T::energy_flag;
        using Base_T::diff_flag;
        using Base_T::metric_flag;
        
    public:
        
        TP_Kernel_FF() = delete;
        
        TP_Kernel_FF( Configurator_T & conf, const Int thread_, const T1 q_half_, const T2 p_half_ )
        :   Base_T               (conf, thread_ )
        ,   q                    (two*q_half_   )
        ,   q_half               (q_half_       )
        ,   q_half_minus_1       (q_half-1      )
        ,   p                    (two*p_half_   )
        ,   p_half               (p_half_       )
        ,   minus_p_half         (-p_half       )
        ,   minus_p_half_minus_1 (-p_half-1     )
        {}
        
        TP_Kernel_FF( TP_Kernel_FF & other, const Int thread_ )
        :   Base_T               (other, thread_            )
        ,   q                    (other.q                   )
        ,   q_half               (other.q_half              )
        ,   q_half_minus_1       (other.q_half_minus_1      )
        ,   p                    (other.p                   )
        ,   p_half               (other.p_half              )
        ,   minus_p_half         (other.minus_p_half        )
        ,   minus_p_half_minus_1 (other.minus_p_half_minus_1)
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
        
        using Base_T::lin_k;
        
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
        
    public:
        
        force_inline Real Compute( const LInt k_global )
        {
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
                    const Int k = lin_k[i][j];
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

            // E = ( |P*(y-x)|^q + |Q*(y-x)|^q) / |y-x|^p
            const Real E = Num * r_minus_p;

            if constexpr ( diff_flag || metric_flag )
            {
                // factor = q / |y-x|^p
                const Real factor = q * r_minus_p;
                // K_xy = q * |P*(y-x)|^(q-2) / |y-x|^p
                const Real K_xy = factor * rCosPhi_q_minus_2;
                // K_yx = q * |Q*(y-x)|^(q-2) / |y-x|^p
                const Real K_yx = factor * rCosPsi_q_minus_2;

                if constexpr ( diff_flag )
                {
                    // H = - p * ( |P*(y-x)|^q + |Q*(y-x)|^q) / |y-x|^(p+2)
                    const Real H = - p * r_minus_p_minus_2 * Num;
                    
                    Real dEdvx = zero;
                    Real dEdvy = zero;
                    
                    for( Int l = 0; l < AMB_DIM; ++l )
                    {
                        dEdv[l] = K_xy * Pv[l] + K_yx * Qv[l] + H * v[l];
                        dEdvx += dEdv[l] * x[l];
                        dEdvy += dEdv[l] * y[l];
                        
                        DX[1+l] -= b * dEdv[l];
                        DY[1+l] += a * dEdv[l];
                    }
                    
                    DX[0] +=  b * ( E - half * factor * rCosPhi_q + dEdvx );
                    DY[0] +=  a * ( E - half * factor * rCosPsi_q - dEdvy );
                    
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
*/
                    
                    const Real K_sym = K_xy + K_yx;
                    
                    const Real b_over_a   = b/a;
                    const Real a_over_b   = a/b;
                    
                    const Real b_over_a_K_xy = b_over_a * K_xy;
                    const Real a_over_b_K_yx = a_over_b * K_yx;
                    
                    ij_block[0]     = - K_sym;
                    ii_block[0][0] +=   b_over_a * K_sym;
                    jj_block[0][0]  =   a_over_b * K_sym;
                    
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
                    
                    copy_buffer<BLOCK_NNZ>( &ij_block[0], &metric_data[BLOCK_NNZ * k_global] );
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
                // We can do an overwrite here.
//                zerofy_buffer<BLOCK_NNZ>( &ij_block[0] );
                
                // We can do an overwrite here.
//                zerofy_buffer<DIAG_NNZ>( &jj_block[0][0] );
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
        
    public:
        
        std::string ClassName() const
        {
            return "TP_Kernel_FF<"
            + this->S.ClassName() + ","
            + TypeName<T1> + ","
            + TypeName<T2> + ","
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(metric_flag) + ","
            ">";
        }
    };

} // namespace Repulsor

#undef BASE
