#pragma once

#include "TangentPointEnergy/TP_Kernel_FF.hpp"
#include "TangentPointEnergy/TP_Kernel_NF.hpp"
#include "TangentPointEnergy/TP_Kernel_VF.hpp"

#include "TangentPointEnergy/TP_Traversor.hpp"

#define CLASS TangentPointObstacleEnergy
#define BASE  EnergyDimRestricted<SimplicialMesh<DOM_DIM_S,AMB_DIM,Real,Int,LInt,SReal,ExtReal>>
#define MESH  SimplicialMesh<DOM_DIM_S,AMB_DIM,Real,Int,LInt,SReal,ExtReal>
#define BESH  SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>
#define ROOT  EnergyBase<SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>>

namespace Repulsor
{
    
    template<typename Mesh_T,int DOM_DIM_T> class CLASS {};
    
    template<
        int DOM_DIM_S, int DOM_DIM_T, int AMB_DIM,
        typename Real, typename Int, typename LInt,
        typename SReal, typename ExtReal
    >
    class CLASS<MESH, DOM_DIM_T> : public BASE
    {
    public:
        
        using Base_T              = BASE;
        using Mesh_T              = typename Base_T::Mesh_T;
        using MeshBase_T          = typename Mesh_T::Base_T;
        using Root_T              = typename Base_T::Root_T;
        
        using ClusterTree_T       = typename Mesh_T::ClusterTree_T;
        
        // We have to explicitly allow an unsymmetric BCT.
        using BlockClusterTree_T  = BlockClusterTree<AMB_DIM,Real,Int,LInt,SReal,ExtReal,false>;
        
        
        using ValueContainer_T    = typename Base_T::ValueContainer_T;
        using TangentVector_T     = typename Base_T::TangentVector_T;
        using CotangentVector_T   = typename Base_T::CotangentVector_T;
        
        using Base_T::Value;
        using Base_T::Differential;
        
        CLASS( const Real q_, const Real p_ )
        :   Base_T ()
        ,   q      ( static_cast<Real>(q_) )
        ,   p      ( static_cast<Real>(p_) )
        {}
        
        virtual ~CLASS() = default;
        
    protected:
        
        const Real q;
        const Real p;
        
    public:
        
        virtual ExtReal value( const Mesh_T & M ) const override
        {
            if( M.InPersistentCacheQ("Obstacle") )
            {
                TP_Traversor<
                    DOM_DIM_S,DOM_DIM_T,BlockClusterTree_T,
                    true,false,false,false
                >
                traversor( M.GetObstacleBlockClusterTree(), this->metric_values, q, p );
                
                return traversor.Compute();
            }
            else
            {
                return Scalar::Zero<ExtReal>;
            }
        }
        
        virtual ExtReal differential( const Mesh_T & M ) const override
        {
            if( M.InPersistentCacheQ("Obstacle") )
            {
                TP_Traversor<
                    DOM_DIM_S,DOM_DIM_T,BlockClusterTree_T,
                    true,true,false,false
                >
                traversor( M.GetObstacleBlockClusterTree(), this->metric_values, q, p );
                
                return traversor.Compute();
            }
            else
            {
                return Scalar::Zero<ExtReal>;
            }
        }
        
        virtual ExtReal density( const Mesh_T & M ) const override
        {
            if( M.InPersistentCacheQ("Obstacle") )
            {
                TP_Traversor<
                    DOM_DIM_S,DOM_DIM_T,BlockClusterTree_T,
                    true,false,false,true
                >
                traversor( M.GetObstacleBlockClusterTree(), this->metric_values, q, p );
                
                return traversor.Compute();
            }
            else
            {
                return Scalar::Zero<ExtReal>;
            }
        }
        
    public:
        
        
        std::string className() const
        {
            return TOOLS_TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM_S)+","+ToString(DOM_DIM_T)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+"("+ToString(q)+","+ToString(p)+")";
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
    };
    
} // namespace Repulsor


//#include "TangentPointEnergy/Make_TangentPointObstacleEnergy.hpp"
#include "TangentPointEnergy/TPO_Factory.hpp"

#undef BESH
#undef MESH
#undef ROOT
#undef BASE
#undef CLASS
