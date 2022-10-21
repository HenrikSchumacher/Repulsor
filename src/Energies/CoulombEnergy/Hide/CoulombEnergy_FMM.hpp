#pragma once

#define CLASS CoulombEnergy_FMM
#define BASE  Energy_FMM<DOM_DIM,AMB_DIM,DEGREE,Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, int DEGREE, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:

        explicit CLASS(
            const Real beta_,
            const ExtReal weight = static_cast<ExtReal>(1)
        )
        :   BASE(
            *Make_CoulombEnergy__NFK    <DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>(beta_),
            *Make_CoulombEnergy__FFK_FMM <AMB_DIM,DEGREE,Real,Int>(beta_),
            weight
                 )
        {}

        virtual ~CLASS() override = default;
                
    public:
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+ToString(DEGREE)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    };
    
}// namespace Repulsor

#undef BASE
#undef CLASS
