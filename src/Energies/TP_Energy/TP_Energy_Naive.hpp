#pragma once

#define CLASS TP_Energy_Naive
#define BASE  Energy_Naive<DOM_DIM,AMB_DIM,Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:
        
        CLASS() {};
        
        CLASS(
            const Real alpha,
            const Real beta,
            const ExtReal weight = static_cast<ExtReal>(1)
        )
        :   BASE(
                 *Make_TP_Energy__NFK<DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>( alpha, beta ),
            weight
            )
        {}

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
