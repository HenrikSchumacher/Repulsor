#pragma once

#define CLASS TangentPointMetric0
#define BASE  MetricDimRestricted<DOM_DIM,AMB_DIM,Real,Int,SReal,ExtReal,OperatorType::MixedOrder>
#define ROOT  MetricBase<Real,Int,SReal,ExtReal>

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
        
        CLASS( const Real q_, const Real p_ )
        :   BASE (                       )
        ,   q    ( static_cast<Real>(q_) )
        ,   p    ( static_cast<Real>(p_) )
        {}
        
        virtual ~CLASS() = default;
        
    protected:
        
        const Real q;
        const Real p;
        
    public:
        
        virtual ValueContainer_T compute_metric( const Mesh_T & M ) const override
        {
            ValueContainer_T metric_values;
         
            TP0_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,false,false,true>
                traversor( M.GetBlockClusterTree(), metric_values, q, p );
         
            (void)traversor.Compute();
         
            return metric_values;
        }
        
        virtual void multiply_metric(
            const Mesh_T & M,
            const bool VF_flag, const bool NF_flag, const bool FF_flag
        ) const override
        {
            TP0_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,false,false,true>
                traversor( M.GetBlockClusterTree(), MetricValues(M), q, p );
            
            (void)traversor.MultiplyMetric(VF_flag,NF_flag,FF_flag);
        }
        
    public:
        
        
        std::string className() const
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+"> ("+ToString(q)+","+ToString(p)+")";
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
    };
    
} // namespace Repulsor


#include "../Energies/TangentPointEnergy/TP_Factory.hpp"

#undef ROOT
#undef BASE
#undef CLASS



