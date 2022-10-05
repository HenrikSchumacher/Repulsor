#pragma once

#define CLASS TP_SingularMetric_FMM
#define BASE  Metric_FMM<DOM_DIM,AMB_DIM,DEGREE,Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, int DEGREE, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:

        using BlockClusterTree_T = typename BASE::BlockClusterTree_T;
        
        CLASS() : BASE() {}
        
        CLASS(
            const BlockClusterTree_T & bct_,
            const Real alpha_,
            const Real beta_,
            const ExtReal weight = static_cast<ExtReal>(1)
        )
        :   BASE(
                bct_,
                *std::make_unique<
                 TP_SingularMetric__NFK   <DOM_DIM,AMB_DIM,       Real,Int,SReal>
                    >(alpha_,beta_),
                 *std::make_unique<
                 TP_SingularMetric__FFK_FMM<DOM_DIM,AMB_DIM,DEGREE,Real,Int>
                 >(alpha_,beta_),
                weight
            )
        {
            this->RequireMetrics();
        }
        
        virtual ~CLASS() override = default;
                
        virtual void ApplyNearFieldKernel( const Real factor, const KernelType type ) const override
        {
            ptic(ClassName()+"::ApplyNearFieldKernel ("+KernelTypeName[type]+")" );
            
            eprint(ClassName()+"::ApplyNearFieldKernel is yet to be implemented.");
//            bct->Near().Multiply_DenseMatrix(
//                factor,
//                near_values.find(type)->second.data(),
//                bct->GetT().PrimitiveInputBuffer().data(),
//                static_cast<Real>(0),
//                bct->GetS().PrimitiveOutputBuffer().data(),
//                bct->GetT().BufferDimension()
//            );
            
            ptoc(ClassName()+"::ApplyNearFieldKernel ("+KernelTypeName[type]+")" );
        }
        
        virtual void ApplyFarFieldKernel( const Real factor, const KernelType type ) const override
        {
            ptic(ClassName()+"::ApplyFarFieldKernel ("+KernelTypeName[type]+")");
            
            eprint(ClassName()+"::ApplyFarFieldKernel is yet to be implemented.");
            
//            bct->Far().Multiply_DenseMatrix(
//                factor,
//                far_values.find(type)->second.data(),
//                bct->GetT().ClusterInputBuffer().data(),
//                static_cast<Real>(0),
//                bct->GetS().ClusterOutputBuffer().data(),
//                bct->GetT().BufferDimension()
//            );
            
            ptoc(ClassName()+"::ApplyFarFieldKernel ("+KernelTypeName[type]+")");
        }
        
    public:
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+ToString(DEGREE)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    };
    
}// namespace Repulsor

#undef BASE
#undef CLASS
