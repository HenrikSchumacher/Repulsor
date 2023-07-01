#pragma once

#include "../Energies/TangentPointEnergy/TP_Traversor.hpp"

#define CLASS TangentPointMetric
#define BASE  MetricDimRestricted<SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>,OperatorType::MixedOrder>
#define MESH  SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>
#define BESH  SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>
#define ROOT  MetricBase<SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>>

namespace Repulsor
{
    template<typename Mesh_T> class CLASS {};
    
    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename LInt, typename SReal, typename ExtReal>
    class CLASS<MESH> : public BASE
    {
    public:
        
        using Mesh_T                  = typename BASE::Mesh_T;
        using BlockClusterTree_T      = typename Mesh_T::BlockClusterTree_T;
        
        using ValueContainer_T        = typename BASE::ValueContainer_T;
        using TangentVector_T         = typename BASE::TangentVector_T;
        using CotangentVector_T       = typename BASE::CotangentVector_T;
        using Values_T                = typename ValueContainer_T::Values_T;
        using BASE::MetricValues;
        
        CLASS( const Real q_, const Real p_ )
        :   BASE ()
        ,   q    ( static_cast<Real>(q_) )
        ,   p    ( static_cast<Real>(p_) )
        {}
        
        virtual ~CLASS() = default;
        
    protected:
        
        const Real q;
        const Real p;
        
    public:
        
        virtual ValueContainer_T compute_metric( const Mesh_T & restrict M ) const override
        {
            ValueContainer_T metric_values;
            
            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,false,false,true>
                traversor( M.GetBlockClusterTree(), metric_values, q, p );

            (void)traversor.Compute();

            return metric_values;
        }
        
        virtual void multiply_metric(
            const Mesh_T & restrict M,
            const bool VF_flag, const bool NF_flag, const bool FF_flag
        ) const override
        {
            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,false,false,true>
                traversor ( M.GetBlockClusterTree(), MetricValues(M), q, p );
            
            (void)traversor.MultiplyMetric(VF_flag,NF_flag,FF_flag);
        }
        
    public:
        
        
        std::string className() const
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">("+ToString(q)+","+ToString(p)+")";
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
    };
    
} // namespace Repulsor


#include "../Energies/TangentPointEnergy/TP_Factory.hpp"

#undef BESH
#undef MESH
#undef ROOT
#undef BASE
#undef CLASS
