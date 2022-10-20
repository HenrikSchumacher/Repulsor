#pragma once

#define CLASS AreaObstacle_FMM_Adaptive
#define BASE  ObstacleEnergy_FMM_Adaptive<DOM_DIM1,DOM_DIM2,AMB_DIM,DEGREE,Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM1, int DOM_DIM2, int AMB_DIM, int DEGREE, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:

        explicit CLASS(
            const ExtReal weight = static_cast<ExtReal>(1),
            const AdaptivitySettings & set_ = AdaptivitySettings()
        )
        :
        BASE(
            *std::make_unique<Area__NFK          <DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Int,SReal>>(),
            *std::make_unique<Area__NFK_Adaptive <DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Int,SReal>>(set_),
            *std::make_unique<Area__FFK_FMM      <AMB_DIM,DEGREE,Real,Int>>(),
            weight
        )
        {}

        virtual ~CLASS() override = default;
                
    public:
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM1)+","+ToString(DOM_DIM2)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    };
    
}// namespace Repulsor

#undef BASE
#undef CLASS
