#pragma once

#include "TangentPointEnergy/TP_Kernel_NF.hpp"

#include "TangentPointEnergy/TP_AllPairs_Traversor.hpp"

#define CLASS TangentPointEnergy_AllPairs
#define BASE  EnergyDimRestricted<DOM_DIM,AMB_DIM,Real,Int,SReal,ExtReal>
#define ROOT  EnergyBase<Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:
        
        using Mesh_T                  = typename BASE::Mesh_T;
        using BlockClusterTree_T      = typename Mesh_T::BlockClusterTree_T;
        using ClusterTree_T           = typename BlockClusterTree_T::ClusterTree_T;
        
        using Values_T                = typename BASE::Values_T;
        using ValueContainer_T        = typename BASE::ValueContainer_T;
        using TangentVector_T         = typename BASE::TangentVector_T;
        using CotangentVector_T       = typename BASE::CotangentVector_T;
        
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
        
        virtual ExtReal value( const Mesh_T & M ) const override
        {
            // Create some dummies.
            ValueContainer_T metric_values;
            
            TP_AllPairs_Traversor<DOM_DIM,DOM_DIM,ClusterTree_T,true,true,false,false>
            traversor( M.GetClusterTree(), M.GetClusterTree(), metric_values, q, p );
            
            return traversor.Compute();
        }
        
        virtual void differential( const Mesh_T & M ) const override
        {
            // Create some dummies.
            ValueContainer_T metric_values;
            
            TP_AllPairs_Traversor<DOM_DIM,DOM_DIM,ClusterTree_T,true,false,true,false>
                traversor( M.GetClusterTree(), M.GetClusterTree(), metric_values, q, p );
            
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



