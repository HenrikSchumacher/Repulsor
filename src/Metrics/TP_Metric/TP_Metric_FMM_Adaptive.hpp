#pragma once

#define CLASS TP_Metric_FMM_Adaptive
#define BASE  Metric_FMM_Adaptive<DOM_DIM,AMB_DIM,DEGREE,Real,Int,SReal,ExtReal>

namespace Repulsion
{
    template<int DOM_DIM, int AMB_DIM, int DEGREE, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:
                                     
        using BlockClusterTree_T = BlockClusterTree<AMB_DIM,Real,Int,SReal,ExtReal>;

        CLASS() : BASE() {};
        
        CLASS(
            const BlockClusterTree_T & bct_,
            const Real alpha_,
            const Real beta_,
            const ExtReal weight = static_cast<ExtReal>(1),
            const AdaptivitySettings & settings_ = AdaptivitySettings()
        )
        :   BASE(
                bct_,
                *std::make_unique<
                    TP_Metric__NFK          <DOM_DIM,AMB_DIM,Real,Int,SReal>
                > (alpha_,beta_),
                *std::make_unique<
                    TP_Metric__NFK_Adaptive <DOM_DIM,AMB_DIM,Real,Int,SReal>
                >(alpha_,beta_,settings_),
                *std::make_unique<
                    TP_Metric__FFK_FMM      <DOM_DIM,AMB_DIM,DEGREE,Real,Int>
                >(alpha_,beta_),
                weight
            )
        {
            this->RequireMetrics();
        }
        
        virtual ~CLASS() override = default;
                
    public:
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+ToString(DEGREE)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    };
    
}// namespace Repulsion

#undef BASE
#undef CLASS
