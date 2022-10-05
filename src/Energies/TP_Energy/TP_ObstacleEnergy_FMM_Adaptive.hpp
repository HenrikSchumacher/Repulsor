#pragma once

#define CLASS TP_ObstacleEnergy_FMM_Adaptive
#define BASE  ObstacleEnergy_FMM_Adaptive<DOM_DIM1,DOM_DIM2,AMB_DIM,DEGREE,Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM1, int DOM_DIM2, int AMB_DIM, int DEGREE, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:

        CLASS(
            const Real alpha_,
            const Real beta_,
            const ExtReal weight = static_cast<ExtReal>(1),
            const AdaptivitySettings & settings_ = AdaptivitySettings()
        )
        :   BASE(
            *Make_TP_Energy__NFK         <DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Int,SReal>( alpha_, beta_ ),
            *Make_TP_Energy__NFK_Adaptive<DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Int,SReal>( alpha_, beta_, settings_),
            *Make_TP_Energy__FFK_FMM     <AMB_DIM,DEGREE,Real,Int>( alpha_, beta_ ),
            weight
            )
        ,   alpha(alpha_)
        ,   beta(beta_)
        {}

        virtual ~CLASS() override = default;
                
    protected:
        
        const Real alpha;
        const Real beta;
        
    public:
        
        Real GetAlpha() const
        {
            return alpha;
        }
        
        Real GetBeta() const
        {
            return beta;
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM1)+","+ToString(DOM_DIM2)+","+ToString(AMB_DIM)+","+ToString(DEGREE)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    };
    
}// namespace Repulsor

#undef BASE
#undef CLASS
