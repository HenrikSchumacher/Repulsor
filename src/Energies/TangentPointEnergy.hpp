#pragma once

#include "TangentPointEnergy/TP_Traversor.hpp"

#define CLASS TangentPointEnergy
#define BASE  EnergyDimRestricted<SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>>
#define MESH  SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>
#define BESH  SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>
#define ROOT  EnergyBase<SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>>

namespace Repulsor
{
    template<typename Mesh_T> class CLASS {};
    
    template<int DOM_DIM,int AMB_DIM, typename Real, typename Int, typename LInt, typename SReal, typename ExtReal>
    class CLASS<MESH> : public BASE
    {
    public:
        
        using Base_T                 = BASE;
        using Mesh_T                 = typename Base_T::Mesh_T;
        using MeshBase_T             = typename Mesh_T::Base_T;
        using Root_T                 = typename Base_T::Root_T;
        
        using BlockClusterTree_T     = typename Mesh_T::BlockClusterTree_T;
        using ClusterTree_T          = typename BlockClusterTree_T::ClusterTree_T;
        
        using ValueContainer_T       = typename Base_T::ValueContainer_T;
        using TangentVector_T        = typename Base_T::TangentVector_T;
        using CotangentVector_T      = typename Base_T::CotangentVector_T;
        
        CLASS( const Real q_, const Real p_ )
        :   Base_T ()
        ,   q ( q_ )
        ,   p ( p_ )
        {}
        
        virtual ~CLASS() = default;
        
    protected:
        
        const Real q;
        const Real p;
        
    public:
        
        virtual ExtReal value( const Mesh_T & M ) const override
        {
            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,true,false,false>
                traversor ( M.GetBlockClusterTree(), this->metric_values, q, p );
            
            return traversor.Compute();
        }
        
        virtual void differential( const Mesh_T & M ) const override
        {
            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,false,true,false>
                traversor( M.GetBlockClusterTree(), this->metric_values, q, p );
            
            (void)traversor.Compute();
        }
        
    public:
        
        std::string className() const
        {
            return std::string("TangentPointEnergy")+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">("+ToString(q)+","+ToString(p)+")";
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
    };
    
} // namespace Repulsor


#include "TangentPointEnergy/TP_Factory.hpp"

#undef BESH
#undef MESH
#undef ROOT
#undef BASE
#undef CLASS
