#pragma once

#define CLASS TP_Kernel_FarField
#define BASE  FMM_Kernel_FarField<ClusterTree_T_,is_symmetric_,energy_flag,diff_flag,hess_flag,metric_flag>

namespace Repulsor
{
    template<
        typename ClusterTree_T_, typename T1, typename T2,
        bool is_symmetric_,
        bool energy_flag, bool diff_flag, bool hess_flag, bool metric_flag
    >
    class CLASS : public BASE
    {
    public:
        
        using ClusterTree_T = ClusterTree_T_;
        
        using Real    = typename ClusterTree_T::Real;
        using Int     = typename ClusterTree_T::Int;
        using SReal   = typename ClusterTree_T::SReal;
        using ExtReal = typename ClusterTree_T::ExtReal;
        
        using BASE::AMB_DIM;
        using BASE::PROJ_DIM;
        using BASE::T_DATA_DIM;
        using BASE::S_DATA_DIM;

        using BASE::S;
        using BASE::T;
        
        using BASE::zero;
        using BASE::one;
        using BASE::two;
        using BASE::is_symmetric;
        
    public:
        
        CLASS() = delete;
        
        CLASS( const ClusterTree_T & S_, const ClusterTree_T & T_, const T1 q_half_, const T2 p_half_ )
        :   BASE( S_, T_ )
        ,   q                    (two*q_half_)
        ,   q_half               (q_half_    )
        ,   q_half_minus_1       (q_half-1   )
        ,   p                    (two*p_half_)
        ,   p_half               (p_half_    )
        ,   minus_p_half         (-p_half    )
        ,   minus_p_half_minus_1 (-p_half-1  )
        {}
        
        CLASS( const CLASS & other )
        :   BASE( other )
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
        
        virtual force_inline Real compute() override
        {
            if( S_ID == T_ID )
            {
                eprint("!!!!");
            }
            
            Real v    [AMB_DIM ] = {};
            Real Pv   [AMB_DIM ] = {};
            Real Qv   [AMB_DIM ] = {};
            Real dEdv [AMB_DIM ] = {};
            Real V    [PROJ_DIM] = {};
            
            Real r2        = zero;
            Real rCosPhi_2 = zero;
            Real rCosPsi_2 = zero;
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
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
                        V[k] = ( one + static_cast<Real>(i!=j) )*v[i]*v[j];
                    }
                }
                rCosPhi_2 += v[i] * Pv[i];
                rCosPsi_2 += v[i] * Qv[i];
            }
            
            const Real rCosPhi_q_minus_2 = MyMath::pow<Real,T1>( fabs(rCosPhi_2), q_half_minus_1);
            const Real rCosPsi_q_minus_2 = MyMath::pow<Real,T1>( fabs(rCosPsi_2), q_half_minus_1);
            const Real r_minus_p_minus_2 = MyMath::pow<Real,T2>( r2, minus_p_half_minus_1 );

            const Real r_minus_p = r_minus_p_minus_2 * r2;
            const Real rCosPhi_q = rCosPhi_q_minus_2 * rCosPhi_2;
            const Real rCosPsi_q = rCosPsi_q_minus_2 * rCosPsi_2;
            const Real Num = ( rCosPhi_q + rCosPsi_q );

            const Real E = Num * r_minus_p;

            if constexpr ( diff_flag )
            {
                const Real factor = q_half * r_minus_p;
                const Real F = factor * rCosPhi_q_minus_2;
                const Real G = factor * rCosPsi_q_minus_2;
                const Real H = - p * r_minus_p_minus_2 * Num;
                
                Real dEdvx = zero;
                Real dEdvy = zero;
                
                for( Int i = 0; i < AMB_DIM; ++i )
                {
                    dEdv[i] = static_cast<Real>(2) * ( F * Pv[i] + G * Qv[i] ) + H * v[i];
                    dEdvx += dEdv[i] * x[i];
                    dEdvy += dEdv[i] * y[i];
                    
                    DX[1+i] -= b * dEdv[i];
                    DY[1+i] += a * dEdv[i];
                }
                
                DX[0] +=  b * ( E - factor * rCosPhi_q + dEdvx );
                DY[0] +=  a * ( E - factor * rCosPsi_q - dEdvy );
                
                const Real  bF =  b * F;
                const Real  aG =  a * G;
                
                for( Int k = 0; k < PROJ_DIM; ++k )
                {
                    DX[1+AMB_DIM+k] +=  bF * V[k];
                    DY[1+AMB_DIM+k] +=  aG * V[k];
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

