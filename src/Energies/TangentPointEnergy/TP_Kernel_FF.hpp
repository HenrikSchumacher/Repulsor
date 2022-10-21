#pragma once

#define CLASS TP_Kernel_FF
#define BASE  FMM_Kernel_FF<ClusterTree_T_,is_symmetric_,energy_flag_,diff_flag_,hess_flag_,metric_flag_,prec_flag_>

namespace Repulsor
{
    template<
        typename ClusterTree_T_, typename T1, typename T2,
        bool is_symmetric_,
        bool energy_flag_, bool diff_flag_, bool hess_flag_, bool metric_flag_, bool prec_flag_
    >
    class CLASS : public BASE
    {
    public:
        
        using ClusterTree_T = ClusterTree_T_;
        
        using Real    = typename ClusterTree_T::Real;
        using Int     = typename ClusterTree_T::Int;
        using SReal   = typename ClusterTree_T::SReal;
        using ExtReal = typename ClusterTree_T::ExtReal;
        using typename BASE::Configurator_T;
        
        using BASE::AMB_DIM;
        using BASE::PROJ_DIM;
        using BASE::T_DATA_DIM;
        using BASE::S_DATA_DIM;
        
        static constexpr Int METRIC_NNZ = 2 + AMB_DIM;
        static constexpr Int PREC_NNZ   = 1;

        using BASE::S;
        using BASE::T;
        
        using BASE::zero;
        using BASE::one;
        using BASE::two;
        using BASE::is_symmetric;
        using BASE::energy_flag;
        using BASE::diff_flag;
        using BASE::hess_flag;
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
        {}
        
        virtual ~CLASS() = default;
        
    protected:

        using BASE::metric_data;
        using BASE::prec_data;
        
        using BASE::S_data;
        using BASE::S_D_data;
        
        using BASE::T_data;
        using BASE::T_D_data;
        
        using BASE::DX;
        using BASE::DY;
        
        using BASE::a;
        using BASE::x;
        using BASE::P;
        
        using BASE::b;
        using BASE::y;
        using BASE::Q;
        
        using BASE::S_ID;
        using BASE::T_ID;
        
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
        
    public:
        
        virtual force_inline Real compute( const Int block_ID ) override
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
                
                if constexpr ( metric_flag )
                {
                    Real * restrict const m_vals = &metric_data[ METRIC_NNZ * block_ID ];
    //                Real * restrict const p_vals = &prec_values  [ PREC_NNZ   * block_ID ];
                    
//                     The metric block looks like this for AMB_DIM == 3:
//
//                          /                                                                 \
//                          |   - K_xy - K_yx     Kyx * v[0]     Kyx * v[1]     Kyx * v[2]    |
//                          |                                                                 |
//                          |   - K_xy * v[0]         0              0              0         |
//                          |                                                                 |
//                          |   - K_xy * v[1]         0              0              0         |
//                          |                                                                 |
//                          |   - K_xy * v[2]         0              0              0         |
//                          |                                                                 |
//                          \                                                                 /
//
//                     This are 1 + 2 * AMB_DIM nonzero values.
//                     With sparse matrix multiplication the bottleneck will be the bandwidth.
//                     Hence instead, we store this data in only 2 + AMB_DIM values by just storing
//
//                              K_xy, K_yx, v[0], ..., v[AMB_DIM-1]!
//
//                     We just must not forget to recompute these values in the matrix multiplication kernel.
                    
                    m_vals[0] = K_xy;
                    m_vals[1] = K_yx;
                    
                    for( Int l = 0; l < AMB_DIM; ++l )
                    {
                        m_vals[2+l] = v[l];
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
        
    public:
        
        virtual Int MetricNonzeroCount() const override
        {
            return METRIC_NNZ;
        }

        virtual Int PreconditionerNonzeroCount() const override
        {
            return PREC_NNZ;
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        
    private:
        
        std::string className() const
        {
            return TO_STD_STRING(CLASS)+"<"
            + S.ClassName() + ","
            + TypeName<T1>::Get() + ","
            + TypeName<T2>::Get() + ","
            + ToString(is_symmetric) + ","
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(hess_flag) + ","
            + ToString(metric_flag) + ","
            ">";
        }
    };

} // namespace Repulsor

#undef BASE
#undef CLASS

