#pragma once

#define CLASS TrivialEnergy_FMM_Adaptive
#define BASE  Energy_FMM_Adaptive<DOM_DIM,AMB_DIM,DEGREE,Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, int DEGREE, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:

        explicit CLASS(
            ExtReal weight = static_cast<ExtReal>(1),
            const AdaptivitySettings & settings_ = AdaptivitySettings()
        )
        :
        BASE(
            *std::make_unique<TrivialEnergy__NFK<DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>>(),
            *std::make_unique<TrivialEnergy__NFK_Adaptive<DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>>(settings_),
            *std::make_unique<TrivialEnergy__FFK_FMM<AMB_DIM,DEGREE,Real,Int>>(),
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
