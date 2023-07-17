#pragma once

namespace Repulsor
{
        
    template<typename Mesh_T, OperatorType op_type> class MetricDimRestricted {};

    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename LInt, typename SReal, typename ExtReal, OperatorType op_type>
    class MetricDimRestricted<SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>, op_type>
        : public MetricBase<SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>>
    {
    public:
        
        using Mesh_T     = SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>;
        using MeshBase_T = typename Mesh_T::Base_T;
        
        using Base_T     = MetricBase<MeshBase_T>;
        using Root_T     = MetricBase<MeshBase_T>;
    

        using TangentVector_T    = typename Base_T::TangentVector_T;
        using CotangentVector_T  = typename Base_T::CotangentVector_T;
        using ValueContainer_T   = typename Base_T::ValueContainer_T;
        
        MetricDimRestricted() = default;

        virtual ~MetricDimRestricted() override = default;
        
        ValueContainer_T & MetricValues( const MeshBase_T & restrict M ) const override
        {
            std::string tag ( ClassName()+"::MetricValues" );

            ptic(tag);
            if( !M.InCacheQ(tag))
            {
                M.SetCache( tag, compute_metric(M) );
            }
            
            ValueContainer_T & restrict result = std::any_cast<ValueContainer_T &>( M.GetCache(tag) );
            
            ptoc(tag);
            
            return result;
        }
        
        ValueContainer_T & MetricValues( const Mesh_T & restrict M ) const
        {
            std::string tag ( ClassName()+"::MetricValues" );

            ptic(tag);
            if( !M.InCacheQ(tag))
            {
                M.SetCache( tag, compute_metric(M) );
            }
            
            ValueContainer_T & restrict result = std::any_cast<ValueContainer_T &>( M.GetCache(tag) );
            
            ptoc(tag);
            
            return result;
        }
        
    
    protected:
        
        ValueContainer_T compute_metric( const MeshBase_T & restrict M ) const
        {
            // Do a down cast and delegate implementation further to descendant class.
            const Mesh_T * restrict Q = dynamic_cast<const Mesh_T *>(&M);
                        
            if( Q != nullptr )
            {
                return compute_metric(*Q);
            }
            else
            {
                eprint(ClassName()+"::compute_metric: Input could not be downcast to compatible type. Doing nothing.");
                
                return ValueContainer_T();
            }
        }
        
        
        // Actual implementation to be specified by descendants.
        virtual ValueContainer_T compute_metric( const Mesh_T & restrict M ) const = 0;
        
        
    public:
        
        void MultiplyMetric(
            cref<MeshBase_T> M,
            cref<ExtReal> alpha, cptr<ExtReal> X,
            cref<ExtReal> beta,  mptr<ExtReal> Y,
            const Int  rhs_count,
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const override
        {
            // Do a down cast and delegate implementation further to descendant class.
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                        
            if( Q != nullptr )
            {
                MultiplyMetric(*Q, alpha, X, beta, Y, rhs_count, VF_flag, NF_flag, FF_flag );
            }
            else
            {
                eprint(ClassName()+"::MultiplyMetric: Input could not be downcast to compatible type. Doing nothing.");
            }
        }
        
        void MultiplyMetric(
            cref<Mesh_T> M,
            cref<ExtReal> alpha, cptr<ExtReal> X,
            cref<ExtReal> beta,  mptr<ExtReal> Y,
            const Int  rhs_count,
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const
        {
            ptic(ClassName()+"::MultiplyMetric");
            auto & restrict S = M.GetBlockClusterTree().GetS();
            auto & restrict T = M.GetBlockClusterTree().GetT();

            T.Pre( X, rhs_count, op_type );
            
            S.RequireBuffers( T.BufferDimension() ); // Tell the S-side what it has to expect.
            
            multiply_metric( M, VF_flag, NF_flag, FF_flag );

            S.Post( Y, alpha, beta, op_type );
            
            ptoc(ClassName()+"::MultiplyMetric");
        }
        
        
        // Actual implementation to be specified by descendants.
        virtual void multiply_metric(
            cref<Mesh_T> M,
            const bool VF_flag, const bool NF_flag, const bool FF_flag
        ) const = 0;

    public:
        
        static std::string className()
        {
            return "MetricDimRestricted<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        
    }; // class MetricDimRestricted

}// namespace Repulsor
