#pragma once

#include "TangentPointEnergy/TP_Kernel_FF.hpp"
#include "TangentPointEnergy/TP_Kernel_NF.hpp"
#include "TangentPointEnergy/TP_Kernel_VF.hpp"
#include "TangentPointEnergy/TP_Traversor.hpp"


#define CLASS TangentPoint
#define BASE  EnergyDimRestricted<DOM_DIM,AMB_DIM,Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:
        
        using Mesh_T                  = typename BASE::Mesh_T;
        using BlockClusterTree_T      = typename Mesh_T::BlockClusterTree_T;
        
        using ValueContainer_T        = typename BASE::ValueContainer_T;
        using Differential_T          = typename BASE::Differential_T;
        using TangentVector_T         = typename BASE::TangentVector_T;
        

        
        template<typename T_q, typename T_p>
        CLASS( const T_q q_, const T_p p_ )
        :   q( static_cast<Real>(q_) )
        ,   p( static_cast<Real>(p_) )
        {}
        
        virtual ~CLASS() = default;
        
    protected:
        
        const Real q;
        const Real p;
        
    public:
        
        virtual ExtReal compute( const Mesh_T & M ) const override
        {
            // Create some dummies.
            std::array<ValueContainer_T,3> metric_values;
            std::array<ValueContainer_T,3> prec_values;
            
            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,false,true,false,false>
                traversor( M.GetBlockClusterTree(), metric_values, prec_values, p, q );
            
            return traversor.Compute();
        }
        
        virtual ExtReal value( const Mesh_T & M ) const override
        {
            // Create some dummies.
            std::array<ValueContainer_T,3> metric_values;
            std::array<ValueContainer_T,3> prec_values;
            
            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,true,false,false,false>
                traversor( M.GetBlockClusterTree(), metric_values, prec_values, p, q );
            
            return traversor.Compute();
        }

        virtual void differential( const Mesh_T & M ) const override
        {
            // Create some dummies.
            std::array<ValueContainer_T,3> metric_values;
            std::array<ValueContainer_T,3> prec_values;
            
            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,false,true,false,false>
                traversor( M.GetBlockClusterTree(), metric_values, prec_values, p, q );
            
            (void)traversor.Compute();
        }
        
        
    public:
        
        
        static std::string className()
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
    };
    
}

#undef CLASS
