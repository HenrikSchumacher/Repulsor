#pragma once

#define CLASS SquaredDistanceMatrix_FMM_Adaptive
#define BASE  Metric_FMM<DOM_DIM,AMB_DIM,DEGREE,Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, int DEGREE, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:
                                     
        using BlockClusterTree_T = BlockClusterTree<AMB_DIM,Real,Int,SReal,ExtReal>;

        CLASS() {}
        
        explicit CLASS(
            const BlockClusterTree_T & bct_,
            const ExtReal weight = static_cast<ExtReal>(1)
        )
        :
        BASE(
            bct_,
            *std::make_unique<
                SquaredDistanceMatrix__NFK_Adaptive<DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>
            >(),
            *std::make_unique<
                SquaredDistanceMatrix__FFK_FMM<AMB_DIM,DEGREE,Real,Int>
            >(),
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
    
}// namespace Repulsor

#undef BASE
#undef CLASS
