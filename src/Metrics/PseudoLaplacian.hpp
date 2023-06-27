#pragma once

#include "PseudoLaplacian/PseudoLaplacian_Traversor.hpp"

#define CLASS PseudoLaplacian

#define BASE  MetricDimRestricted<                                          \
        SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>,        \
        high_order ? OperatorType::HighOrder : OperatorType::LowOrder       \
    >

#define MESH  SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>
#define BESH  SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>
#define ROOT  MetricBase<SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>>

namespace Repulsor
{
    template<typename Mesh_T, bool high_order> class CLASS {};
    
    template<int DOM_DIM, int AMB_DIM,  typename Real, typename Int, typename LInt, typename SReal, typename ExtReal, bool high_order>
    class CLASS<MESH,high_order> : public BASE
    {
    private:
        
        using Base_T = BASE;
        
    public:
        
        using Mesh_T                  = typename Base_T::Mesh_T;
        using BlockClusterTree_T      = typename Mesh_T::BlockClusterTree_T;
        
        using TangentVector_T         = typename Base_T::TangentVector_T;
        using CotangentVector_T       = typename Base_T::CotangentVector_T;
        using ValueContainer_T        = typename Base_T::ValueContainer_T;
        using Values_T                = typename ValueContainer_T::Values_T;
        
        using Base_T::MetricValues;
        
        explicit CLASS( const Real s_ )
        :   Base_T ()
        ,   s ( s_ )
        {}
        
        virtual ~CLASS() = default;
        
    protected:
        
        const Real s;
        
    public:
        
        virtual ValueContainer_T compute_metric( const Mesh_T & M ) const override
        {
            ValueContainer_T metric_values;
            
            PseudoLaplacian_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,high_order>
                traversor( M.GetBlockClusterTree(), metric_values, s );

            (void)traversor.Compute();
            
            return metric_values;
        }
        
        virtual void multiply_metric(
            const Mesh_T & M,
            const bool VF_flag, const bool NF_flag, const bool FF_flag
        ) const override
        {
            PseudoLaplacian_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,high_order>
                traversor ( M.GetBlockClusterTree(), MetricValues(M), s );
            
            (void)traversor.MultiplyMetric(VF_flag,NF_flag,FF_flag);
        }
        
    public:
        
        
        std::string className() const
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">("+ToString(s)+")";
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
    };
    
} // namespace Repulsor


#include "PseudoLaplacian/PseudoLaplacian_Factory.hpp"

#undef BESH
#undef MESH
#undef ROOT
#undef BASE
#undef CLASS
