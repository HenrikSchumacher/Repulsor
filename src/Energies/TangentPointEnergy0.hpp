#pragma once

#include "TangentPointEnergy0/TP0_Traversor.hpp"

#define CLASS TangentPointEnergy0
#define BASE  EnergyDimRestricted<SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>>
#define MESH  SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>
#define BESH  SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>
#define ROOT  EnergyBase<SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>>

namespace Repulsor
{
    template<typename Mesh_T> class CLASS {};
    
    template<int DOM_DIM,int AMB_DIM,typename Real, typename Int, typename LInt, typename SReal, typename ExtReal>
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
        
        using Base_T::Value;
        using Base_T::Differential;
        
        CLASS( const Real p_ )
        :   Base_T ()
        ,   q ( p_ )
        ,   p ( Scalar::Two<Real> * p_ )
        {}
        
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
            TP0_Traversor<
                DOM_DIM,DOM_DIM,BlockClusterTree_T,
                true,false,false,false
            >
            traversor( M.GetBlockClusterTree(), this->metric_values, q, p );
            
            return traversor.Compute();
        }
        
        virtual ExtReal differential( const Mesh_T & M ) const override
        {
            TP0_Traversor<
                DOM_DIM,DOM_DIM,BlockClusterTree_T,true,true,false,false
            >
            traversor( M.GetBlockClusterTree(), this->metric_values, q, p );
            
            return traversor.Compute();
        }
        
        virtual ExtReal density( const Mesh_T & M ) const override
        {
            TP0_Traversor<
                DOM_DIM,DOM_DIM,BlockClusterTree_T,
                true,false,false,true
            >
            traversor( M.GetBlockClusterTree(), this->metric_values, q, p );
            
            return traversor.Compute();
        }
        
    public:
        
        std::string className() const
        {
            return TOOLS_TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">("+ToString(q)+","+ToString(p)+")";
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
