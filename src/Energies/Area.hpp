#pragma once

#include "Area/Area_Kernel_FF.hpp"
#include "Area/Area_Kernel_NF.hpp"
#include "Area/Area_Kernel_VF.hpp"
#include "Area/Area_Traversor.hpp"


#define CLASS Area
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
        
        using ValueContainer_T        = typename BASE::ValueContainer_T;
        using Differential_T          = typename BASE::Differential_T;
        using TangentVector_T         = typename BASE::TangentVector_T;
        
        
        
        CLASS( const Real weight_ )
        :   BASE ( weight_               )
        {}
        
        virtual ~CLASS() = default;
        
    protected:
        
    public:
        
        virtual ExtReal compute( const Mesh_T & M ) const override
        {
            // Create some dummies.
            std::array<ValueContainer_T,3> metric_values;
            std::array<ValueContainer_T,3> prec_values;
            
            Area_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,true,true,false,false>
            traversor( M.GetBlockClusterTree(), metric_values, prec_values );
            
            return traversor.Compute();
        }
        
        virtual ExtReal value( const Mesh_T & M ) const override
        {
            // Create some dummies.
            std::array<ValueContainer_T,3> metric_values;
            std::array<ValueContainer_T,3> prec_values;
            
            Area_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,true,false,false,false>
                traversor( M.GetBlockClusterTree(), metric_values, prec_values );
            
            return traversor.Compute();
        }
        
        virtual void differential( const Mesh_T & M ) const override
        {
            // Create some dummies.
            std::array<ValueContainer_T,3> metric_values;
            std::array<ValueContainer_T,3> prec_values;
            
            Area_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,false,true,false,false>
                traversor( M.GetBlockClusterTree(), metric_values, prec_values );
            
            (void)traversor.Compute();
        }
        
        
    public:
        
        
        std::string className() const
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+"("+")";
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
    };
    
} // namespace Repulsor


#include "Area/Make_Area.hpp"

#undef ROOT
#undef BASE
#undef CLASS



