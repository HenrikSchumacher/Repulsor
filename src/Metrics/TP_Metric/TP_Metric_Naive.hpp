#pragma once

#define CLASS TP_Metric_Naive
#define BASE  Metric_Naive<DOM_DIM,AMB_DIM,Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:
        
        using ClusterTree_T = ClusterTree<AMB_DIM,Real,Int,SReal,ExtReal>;

        CLASS() {}
        
        CLASS(  const ClusterTree_T & S_,
                const ClusterTree_T & T_,
                const Real alpha,
                const Real beta,
                const ExtReal weight = static_cast<Real>(1)
        )
        :
        BASE(
            S_,T_,
            *std::make_unique<TP_Metric__NFK<DOM_DIM,AMB_DIM,Real,Int,SReal>>( alpha, beta ),
            weight
        )
        {
            ptic(ClassName());
            
            this->RequireMetrics();
            
            ptoc(ClassName());
        }
        
        virtual ~CLASS() override = default;
                
    public:
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    };
    
}// namespace Repulsor

#undef BASE
#undef CLASS
