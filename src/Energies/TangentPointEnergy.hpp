#pragma once

#include "TangentPointEnergy/TP_Kernel_FF.hpp"
#include "TangentPointEnergy/TP_Kernel_NF.hpp"
#include "TangentPointEnergy/TP_Kernel_VF.hpp"

#include "TangentPointEnergy/TP_Kernel_FF_MultiplyMetric.hpp"
#include "TangentPointEnergy/TP_Kernel_NF_MultiplyMetric.hpp"
#include "TangentPointEnergy/TP_Kernel_VF_MultiplyMetric.hpp"

#include "TangentPointEnergy/TP_Traversor.hpp"

#define CLASS TangentPointEnergy
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
        
        using Values_T                = typename BASE::Values_T;
        using ValueContainer_T        = typename BASE::ValueContainer_T;
        using TangentVector_T         = typename BASE::TangentVector_T;
        using CotangentVector_T       = typename BASE::CotangentVector_T;
        
        using BASE::MetricValues;
        
        CLASS( const Real weight_, const Real q_, const Real p_ )
        :   BASE ( weight_               )
        ,   q    ( static_cast<Real>(q_) )
        ,   p    ( static_cast<Real>(p_) )
        {}
        
        virtual ~CLASS() = default;
        
    protected:
        
        const Real q;
        const Real p;
        
    public:
        
        virtual ExtReal compute( const Mesh_T & M ) const override
        {
            // Create some dummies.
            ValueContainer_T metric_values;
            ValueContainer_T prec_values;
            
            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,true,true,false,false,false>
            traversor( M.GetBlockClusterTree(), metric_values, prec_values, q, p );
            
            return traversor.Compute();
        }
        
        virtual ExtReal value( const Mesh_T & M ) const override
        {
            // Create some dummies.
            ValueContainer_T metric_values;
            ValueContainer_T prec_values;
            
            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,true,false,false,false,false>
                traversor( M.GetBlockClusterTree(), metric_values, prec_values, q, p );
            
            return traversor.Compute();
        }
        
        virtual void differential( const Mesh_T & M ) const override
        {
            // Create some dummies.
            ValueContainer_T metric_values;
            ValueContainer_T prec_values;
            
            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,false,true,false,false,false>
                traversor( M.GetBlockClusterTree(), metric_values, prec_values, q, p );
            
            (void)traversor.Compute();
        }
        
        virtual ValueContainer_T compute_metric( const Mesh_T & M ) const override
        {
            ValueContainer_T metric_values;
            
            // Create some dummies.
            ValueContainer_T prec_values;
        
            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,false,false,false,true,false>
                traversor( M.GetBlockClusterTree(), metric_values, prec_values, q, p );
            
            (void)traversor.Compute();
            
            return metric_values;
        }
        
        virtual void multiply_metric( const Mesh_T & M ) const override
        {
            // Create some dummies.
            ValueContainer_T prec_values;
        
            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,false,false,false,true,false>
                traversor( M.GetBlockClusterTree(), MetricValues(M), prec_values, q, p );
            
            (void)traversor.MultiplyMetric();
        }
        
    public:
        
        
        std::string className() const
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+"("+ToString(q)+","+ToString(p)+")";
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
    };
    
} // namespace Repulsor


#include "TangentPointEnergy/Make_TangentPointEnergy.hpp"

#undef ROOT  
#undef BASE
#undef CLASS


