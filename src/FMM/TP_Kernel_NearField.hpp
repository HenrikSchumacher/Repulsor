#pragma once

#define CLASS TP_Kernel_NearField
#define  BASE  FMM_Kernel_NearField<S_DOM_DIM_,T_DOM_DIM_,ClusterTree_T_,energy_flag,diff_flag,hess_flag,metric_flag>

namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename ClusterTree_T_, typename T1, typename T2,
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
        using BASE::S_DOM_DIM;
        using BASE::T_DOM_DIM;
        using BASE::S_COORD_DIM;
        using BASE::T_COORD_DIM;
        using BASE::T_DATA_DIM;
        using BASE::S_DATA_DIM;
        
        using BASE::S_Tree;
        using BASE::T_Tree;
        
        using BASE::zero;
        using BASE::one;
        using BASE::two;
    public:
        
        CLASS() = delete;
        
        CLASS( const ClusterTree_T & S_, const ClusterTree_T & T_, const T1 alphahalf_, const T2 betahalf_ )
        :   BASE( S_, T_ )
        ,   alpha(two*alphahalf_)
        ,   alphahalf(alphahalf_)
        ,   alphahalf_minus_1(alphahalf-1)
        ,   beta(two*betahalf_)
        ,   betahalf(betahalf_)
        ,   minus_betahalf(-betahalf)
        ,   minus_betahalf_minus_1(-betahalf-1)
        {}
        
        CLASS( const CLASS & other )
        :   BASE(other)
        ,   alpha(other.alpha)
        ,   alphahalf(other.alphahalf)
        ,   alphahalf_minus_1(other.alphahalf_minus_1)
        ,   beta(other.beta)
        ,   betahalf(other.betahalf)
        ,   minus_betahalf(other.minus_betahalf)
        ,   minus_betahalf_minus_1(other.minus_betahalf_minus_1)
        {}
        
        virtual ~CLASS() = default;
        
    protected:
        
        using BASE::S_near;
        using BASE::S_D_near;
        
        using BASE::T_near;
        using BASE::T_D_near;
        
        using BASE::DX;
        using BASE::DY;
        
        using BASE::a;
        using BASE::x;
        using BASE::p;
        
        using BASE::b;
        using BASE::y;
        using BASE::q;
        
        using BASE::S_ID;
        using BASE::T_ID;
        
        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
        using BASE::S_scale;
        using BASE::T_scale;
        
        const Real alpha;
        const T1 alphahalf;
        const T1 alphahalf_minus_1;
        
        const Real beta;
        const T2 betahalf;
        const T2 minus_betahalf;
        const T2 minus_betahalf_minus_1;
        
    public:
        
        virtual Real Compute() override
        {
            Real v    [AMB_DIM ] = {};
            Real Pv   [AMB_DIM ] = {};
            Real Qv   [AMB_DIM ] = {};
            Real dEdv [AMB_DIM ] = {};
            Real V    [PROJ_DIM] = {};
            
            Real r2       = zero;
            Real rCosPhi2 = zero;
            Real rCosPsi2 = zero;
            
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
                    Pv[i] += p[k] * v[j];
                    Qv[i] += q[k] * v[j];
                    if( j >= i )
                    {
                        V[k] = ( one + static_cast<Real>(i!=j) )*v[i]*v[j];
                    }
                }
                rCosPhi2 += v[i] * Pv[i];
                rCosPsi2 += v[i] * Qv[i];
            }
            
            const Real rCosPhiAlphaMinus2 = MyMath::pow<Real,T1>( fabs(rCosPhi2), alphahalf_minus_1);
            const Real rCosPsiAlphaMinus2 = MyMath::pow<Real,T1>( fabs(rCosPsi2), alphahalf_minus_1);
            const Real rMinusBetaMinus2   = MyMath::pow<Real,T2>( r2, minus_betahalf_minus_1 );

            const Real rMinusBeta   = rMinusBetaMinus2 * r2;
            const Real rCosPhiAlpha = rCosPhiAlphaMinus2 * rCosPhi2;
            const Real rCosPsiAlpha = rCosPsiAlphaMinus2 * rCosPsi2;
            const Real Num = ( rCosPhiAlpha + rCosPsiAlpha );

            const Real E = Num * rMinusBeta;

            if constexpr ( diff_flag )
            {
                const Real factor = alphahalf * rMinusBeta;
                const Real F = factor * rCosPhiAlphaMinus2;
                const Real G = factor * rCosPsiAlphaMinus2;
                const Real H = - beta * rMinusBetaMinus2 * Num;
                
                Real dEdvx = zero;
                Real dEdvy = zero;
                
                for( Int i = 0; i < AMB_DIM; ++i )
                {
                    dEdv[i] = two * ( F * Pv[i] + G * Qv[i] ) + H * v[i];
                    dEdvx += dEdv[i] * x[i];
                    dEdvy += dEdv[i] * y[i];
                    
                    for( Int ii = 0; ii < S_DOM_DIM+1; ++ii )
                    {
                        DX[1+AMB_DIM*ii+i] -=    S_scale *  b * dEdv[i];
                    }
                    
                    for( Int ii = 0; ii < T_DOM_DIM+1; ++ii )
                    {
                        DY[1+AMB_DIM*ii+i] +=    T_scale *  a * dEdv[i];
                    }
                }
                
                DX[0] +=  b * ( E - factor * rCosPhiAlpha + dEdvx );
                DY[0] +=  a * ( E - factor * rCosPsiAlpha - dEdvy );
                
                const Real  bF =  b * F;
                const Real  aG =  a * G;
                
                for( Int k = 0; k < PROJ_DIM; ++k )
                {
                    DX[1+S_COORD_DIM+k] +=  bF * V[k];
                    DY[1+T_COORD_DIM+k] +=  aG * V[k];
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
            return TO_STD_STRING(CLASS)+"<"+this->tree_string+">";
        }

    };

} // namespace Repulsor

#undef BASE
#undef CLASS

