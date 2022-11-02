#pragma once

#define CLASS TP0_Kernel_FF
#define BASE  FMM_Kernel_FF<BlockClusterTree_T_,energy_flag_,diff_flag_,metric_flag_>

namespace Repulsor
{
    template<
        int S_DOM_DIM, int T_DOM_DIM,
        typename BlockClusterTree_T_, typename T1, typename T2,
        bool energy_flag_, bool diff_flag_, bool metric_flag_
    >
    class CLASS : public BASE
    {
    public:
        
        using BlockClusterTree_T = typename BASE::BlockClusterTree_T;
        
        using ClusterTree_T      = typename BASE::ClusterTree_T;
        using Values_T           = typename BASE::Values_T;
        using ValueContainer_T   = typename BASE::ValueContainer_T;
        
        using Real               = typename BASE::Real;
        using SReal              = typename BASE::SReal;
        using ExtReal            = typename BASE::ExtReal;
        using Int                = typename BASE::Int;
        using LInt               = typename BASE::LInt;
        
        using Configurator_T     = typename BASE::Configurator_T;
        
        using BASE::AMB_DIM;
        using BASE::PROJ_DIM;
        using BASE::T_DATA_DIM;
        using BASE::S_DATA_DIM;
        
        static constexpr Int ROWS      = 1 + AMB_DIM;
        static constexpr Int COLS      = 1 + AMB_DIM;
        static constexpr Int BLOCK_NNZ = 2;
        static constexpr Int DIAG_NNZ  = 2;


        using BASE::S;
        using BASE::T;
        
        using BASE::zero;
        using BASE::one;
        using BASE::two;
        using BASE::is_symmetric;
        using BASE::energy_flag;
        using BASE::diff_flag;
        using BASE::metric_flag;
        
    public:
        
        CLASS() = delete;
        
        CLASS( Configurator_T & conf, const T1 q_half_, const T2 p_half_ )
        :   BASE                 (conf       )
        ,   q                    (two*q_half_)
        ,   q_half               (q_half_    )
        ,   q_half_minus_1       (q_half-1   )
        ,   p                    (two*p_half_)
        ,   p_half               (p_half_    )
        ,   minus_p_half         (-p_half    )
        ,   minus_p_half_minus_1 (-p_half-1  )
        ,   s                    ((p - S_DOM_DIM) / q)
        ,   s_exp                (static_cast<Real>(-0.5) * (static_cast<Real>(2) * (s - static_cast<Real>(1)) + S_DOM_DIM))
        {}
        
        CLASS( const CLASS & other )
        :   BASE                 (other                     )
        ,   q                    (other.q                   )
        ,   q_half               (other.q_half              )
        ,   q_half_minus_1       (other.q_half_minus_1      )
        ,   p                    (other.p                   )
        ,   p_half               (other.p_half              )
        ,   minus_p_half         (other.minus_p_half        )
        ,   minus_p_half_minus_1 (other.minus_p_half_minus_1)
        ,   s                    ((p - S_DOM_DIM) / q)
        ,   s_exp                (static_cast<Real>(-0.5) * (static_cast<Real>(2) * (s - static_cast<Real>(1)) + S_DOM_DIM))
        {}
        
        ~CLASS() = default;
        
    protected:

        using BASE::metric_data;
        
        using BASE::S_data;
        using BASE::S_D_data;
        
        using BASE::T_data;
        using BASE::T_D_data;
        
        using BASE::S_diag;
        using BASE::T_diag;
        
        using BASE::DX;
        using BASE::DY;
        
        using BASE::a;
        using BASE::x;
        using BASE::P;
        
        using BASE::b;
        using BASE::y;
        using BASE::Q;
        
        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
        const Real q;
        const T1   q_half;
        const T1   q_half_minus_1;
        
        const Real p;
        const T2   p_half;
        const T2   minus_p_half;
        const T2   minus_p_half_minus_1;
        
        const Real s;
        const Real s_exp;
        
        Real ij_block [BLOCK_NNZ] = {};
        Real ii_block [BLOCK_NNZ] = {};
        Real jj_block [BLOCK_NNZ] = {};
     
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
            
            Real result = 0;
            
            if constexpr ( energy_flag || diff_flag )
            {
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
                
                if constexpr ( energy_flag )
                {
                    result = this->symmetry_factor * a * E * b;
                }
                
                if constexpr ( diff_flag )
                {
                    // Needed for both the differential and the metric.
                    
                    const Real factor = q_half * r_minus_p;         // = q / |y-x|^p
                    
                    const Real K_xy = factor * rCosPhi_q_minus_2;   // = |P*(y-x)|^(q-2) / |y-x|^p
                    const Real K_yx = factor * rCosPsi_q_minus_2;   // = |Q*(y-x)|^(q-2) / |y-x|^p
                    
                    const Real H = - p * r_minus_p_minus_2 * Num;
                    
                    Real dEdvx = zero;
                    Real dEdvy = zero;
                    
                    for( Int l = 0; l < AMB_DIM; ++l )
                    {
                        dEdv[l] = static_cast<Real>(2) * ( K_xy * Pv[l] + K_yx * Qv[l] ) + H * v[l];
                        dEdvx += dEdv[l] * x[l];
                        dEdvy += dEdv[l] * y[l];
                        
                        DX[1+l] -= b * dEdv[l];
                        DY[1+l] += a * dEdv[l];
                    }
                    
                    DX[0] +=  b * ( E - factor * rCosPhi_q + dEdvx );
                    DY[0] +=  a * ( E - factor * rCosPsi_q - dEdvy );
                    
                    const Real  b_K_xy =  b * K_xy;
                    const Real  a_K_yx =  a * K_yx;
                    
                    for( Int k = 0; k < PROJ_DIM; ++k )
                    {
                        DX[1+AMB_DIM+k] +=  b_K_xy * V[k];
                        DY[1+AMB_DIM+k] +=  a_K_yx * V[k];
                    }
                    
                }
            }
            
            if constexpr ( metric_flag )
            {
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
    
                const Real r4 = r2 * r2;

                // The following line makes up approx 2/3 of this function's runtime! This is why we avoid pow as much as possible and replace it with MyMath::pow.;
                // I got it down to this single call to pow. We might want to generate a lookup table for it...;
                // The factor of (-2.) is here, because we assemble the _metric_, not the kernel.;
                const Real a_1 = static_cast<Real>(-2) * MyMath::pow(r2, s_exp);
                
                const Real a_0 = static_cast<Real>(0.5) * (rCosPhi_2 + rCosPsi_2) / r4 * a_1;
                
                const Real b_over_a = b/a;
                const Real a_over_b = a/b;
                
                ij_block[0]  =   a_0;
                ii_block[0] -=   b_over_a * a_0;
                jj_block[0]  = - a_over_b * a_0;
                
                ij_block[1]  =   a_1;
                ii_block[1] -=   b_over_a * a_1;
                jj_block[1]  = - a_over_b * a_1;
                
                copy_buffer( &ij_block[0], &metric_data[BLOCK_NNZ * k_global], BLOCK_NNZ );
            }
            
            return result;
        }
        
        
        force_inline void LoadS( const Int i_global )
        {
            this->loadS( i_global );
            
            if constexpr ( metric_flag )
            {
                zerofy_buffer( &ii_block[0][0], DIAG_NNZ );
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
//                zerofy_buffer( &ij_block[0], BLOCK_NNZ );
                
                // We can do an overwrite here.
//                zerofy_buffer( &jj_block[0][0], DIAG_NNZ );
            }
        }

        force_inline void WriteT( const Int j_global )
        {
            if constexpr ( metric_flag )
            {
                add_to_buffer<DIAG_NNZ>( &jj_block[0][0], &T_diag[DIAG_NNZ * j_global] );
            }
        }

        std::string className() const
        {
            return TO_STD_STRING(CLASS)+"<"
            + ToString(S_DOM_DIM) + ","
            + ToString(T_DOM_DIM) + ","
            + S.ClassName() + ","
            + TypeName<T1>::Get() + ","
            + TypeName<T2>::Get() + ","
            + ToString(is_symmetric) + ","
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(metric_flag) + ","
            ">";
        }
    };

} // namespace Repulsor

#undef BASE
#undef CLASS

