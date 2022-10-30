#pragma once

#define CLASS MetricDimRestricted
#define BASE  MetricBase<Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:
        
        using MeshBase_T         = typename BASE::MeshBase_T;

        using TangentVector_T    = typename BASE::TangentVector_T;
        using CotangentVector_T  = typename BASE::CotangentVector_T;
        
        using Values_T           = typename BASE::Values_T;
        using ValueContainer_T   = typename BASE::ValueContainer_T;
        
        using Mesh_T             = SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,SReal,ExtReal>;
        
        CLASS() = default;

        virtual ~CLASS() override = default;
        
        // Do a down cast and delegate implementation further to descendant class.
        ValueContainer_T compute_metric( const MeshBase_T & M ) const override
        {
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                        
            if( Q != nullptr )
            {
                return compute_metric(*Q);
            }
            else
            {
                eprint(ClassName()+"::differential: Input could not be downcast to compatible type. Doing nothing.");
                
                return ValueContainer_T();
            }
        }
        
        // Actual implementation to be specified by descendants.
        virtual ValueContainer_T compute_metric( const Mesh_T & M ) const = 0;

        // Actual implementation to be specified by descendants.
        virtual void multiply_metric(
            const MeshBase_T & M,
            const bool VF_flag, const bool NF_flag, const bool FF_flag
        ) const override
        {
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                        
            if( Q != nullptr )
            {
                multiply_metric(*Q, VF_flag, NF_flag, FF_flag );
            }
            else
            {
                eprint(ClassName()+"::multiply_metric: Input could not be downcast to compatible type. Doing nothing.");
            }
        }
        
        // Actual implementation to be specified by descendants.
        virtual void multiply_metric(
            const Mesh_T & M,
            const bool VF_flag, const bool NF_flag, const bool FF_flag
        ) const = 0;

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

}// namespace Repulsor

#undef BASE
#undef CLASS
