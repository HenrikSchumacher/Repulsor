#pragma once

#include "TangentPointEnergy0/TP0_Traversor.hpp"

#define CLASS TangentPointEnergy0
#define BASE  EnergyDimRestricted<DOM_DIM,AMB_DIM,Real_,Int_,SReal_,ExtReal_>
#define ROOT  EnergyBase<Real_,Int_,SReal_,ExtReal_>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, typename Real_, typename Int_, typename SReal_, typename ExtReal_>
    class CLASS : public BASE
    {
    public:
        
        using Real    = Real_;
        using Int     = Int_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        
        using Mesh_T                  = typename BASE::Mesh_T;
        using BlockClusterTree_T      = typename Mesh_T::BlockClusterTree_T;
        using ClusterTree_T           = typename BlockClusterTree_T::ClusterTree_T;
        
        using Values_T                = typename BASE::Values_T;
        using ValueContainer_T        = typename BASE::ValueContainer_T;
        using TangentVector_T         = typename BASE::TangentVector_T;
        using CotangentVector_T       = typename BASE::CotangentVector_T;
        
        CLASS( const Real q_, const Real p_ )
        :   BASE ()
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
            // Create some dummies.
            ValueContainer_T metric_values;
            
            TP0_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,true,false,false>
                traversor( M.GetBlockClusterTree(), metric_values, q, p );
            
            return traversor.Compute();
        }
        
        virtual void differential( const Mesh_T & M ) const override
        {
            // Create some dummies.
            ValueContainer_T metric_values;
            
            TP0_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,false,true,false>
                traversor( M.GetBlockClusterTree(), metric_values, q, p );
            
            (void)traversor.Compute();
        }
        
    public:
        
        std::string className() const
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+","+TypeName<ExtReal>+">("+ToString(q)+","+ToString(p)+")";
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
    };
    
} // namespace Repulsor


#include "TangentPointEnergy/TP_Factory.hpp"

#undef ROOT  
#undef BASE
#undef CLASS
