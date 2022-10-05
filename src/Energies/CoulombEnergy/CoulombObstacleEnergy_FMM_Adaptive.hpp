#pragma once

#define CLASS CoulombObstacleEnergy_FMM_Adaptive
#define BASE  ObstacleEnergy_FMM_Adaptive<DOM_DIM1,DOM_DIM2,AMB_DIM,DEGREE,Real,Int,SReal,ExtReal>

namespace Repulsion
{
    template<int DOM_DIM1, int DOM_DIM2, int AMB_DIM, int DEGREE, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:

        explicit CLASS(
            const Real beta_,
            ExtReal weight = static_cast<ExtReal>(1),
            const AdaptivitySettings & set_ = AdaptivitySettings()
        )
        :   BASE(
                 *Make_CoulombEnergy__NFK          <DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Int,SReal> (beta_),
                 *Make_CoulombEnergy__NFK_Adaptive <DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Int,SReal> (beta_,set_),
                 *Make_CoulombEnergy__FFK_FMM      <AMB_DIM,DEGREE,Real,Int>(beta_),
                 weight
            )
        {}

        virtual ~CLASS() override = default;
              
    public:

        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM1)+","+ToString(DOM_DIM2)+","+ToString(AMB_DIM)+","+ToString(DEGREE)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    };
    
}// namespace Repulsion

#undef BASE
#undef CLASS
