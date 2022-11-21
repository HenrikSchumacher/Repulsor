#pragma once

#define CLASS MetricDimRestricted
#define BASE  MetricBase<Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<
        int DOM_DIM, int AMB_DIM,
        typename Real, typename Int, typename SReal, typename ExtReal,
        OperatorType op_type
    >
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
        
        ValueContainer_T & MetricValues( const MeshBase_T & M ) const override
        {
            ptic(ClassName()+"::MetricValues");
            if( !M.IsCached(ClassName()+"::MetricValues"))
            {
                std::any thing ( std::move(compute_metric(M)) );
                
                M.SetCache( ClassName()+"::MetricValues", thing );
            }

            ptoc(ClassName()+"::MetricValues");

            auto & result = std::any_cast<ValueContainer_T &>(
                  M.GetCache(ClassName()+"::MetricValues")
            );
            
            return result;
        }
    
    protected:
        
        
        ValueContainer_T compute_metric( const MeshBase_T & M ) const
        {
            // Do a down cast and delegate implementation further to descendant class.
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
        
    public:
        
        // Actual implementation to be specified by descendants.
        virtual ValueContainer_T compute_metric( const Mesh_T & M ) const = 0;
        
    public:
        
        void MultiplyMetric(
            const MeshBase_T &             M,
            const ExtReal                  alpha,
            const ExtReal * restrict const X,
            const ExtReal                  beta,
                  ExtReal * restrict const Y,
            const Int                      rhs_count,
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const override
        {
            ptic(ClassName()+"::MultiplyMetric");
            auto & S = M.GetBlockClusterTree().GetS();
            auto & T = M.GetBlockClusterTree().GetT();

            T.Pre( X, rhs_count, op_type );
            
            S.RequireBuffers( T.BufferDimension() ); // Tell the S-side what it has to expect.
            
            // Do a down cast and delegate implementation further to descendant class.
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                        
            if( Q != nullptr )
            {
                multiply_metric(*Q, VF_flag, NF_flag, FF_flag );
            }
            else
            {
                eprint(ClassName()+"::multiply_metric: Input could not be downcast to compatible type. Doing nothing.");
            }

            S.Post( Y, alpha, beta, op_type );
            
            ptoc(ClassName()+"::MultiplyMetric");
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
