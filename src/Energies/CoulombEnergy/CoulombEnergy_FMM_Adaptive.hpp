#pragma once

#define CLASS CoulombEnergy_FMM_Adaptive
#define BASE  Energy_FMM_Adaptive<DOM_DIM,AMB_DIM,DEGREE,Real,Int,SReal,ExtReal>

namespace Repulsion
{
    template<int DOM_DIM, int AMB_DIM, int DEGREE, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:

        explicit CLASS(
            const Real beta_,
            const ExtReal weight = static_cast<ExtReal>(1),
            const AdaptivitySettings & set_ = AdaptivitySettings()
        )
        :
        BASE(
            *Make_CoulombEnergy__NFK          <DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>(beta_),
            *Make_CoulombEnergy__NFK_Adaptive <DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>(beta_,set_),
            *Make_CoulombEnergy__FFK_FMM      <AMB_DIM,DEGREE,Real,Int>(beta_),
            weight
        )
        {
//            wprint("Far field of CoulombEnergy_FMM_Adaptive is deactivated.");
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
