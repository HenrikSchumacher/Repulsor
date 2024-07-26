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
        using CG_T               = typename Base_T::CG_T;
        
        MetricDimRestricted() = default;

        virtual ~MetricDimRestricted() override = default;
        
        mref<ValueContainer_T> MetricValues( cref<MeshBase_T> M ) const override
        {
            std::string tag ( ClassName()+"::MetricValues" );

            ptic(tag);
            if( !M.InCacheQ(tag))
            {
                M.SetCache( tag, compute_metric(M) );
            }
            ptoc(tag);
            
            return M.template GetCache<ValueContainer_T>(tag);
        }
        
        mref<ValueContainer_T> MetricValues( cref<Mesh_T> M ) const
        {
            std::string tag ( ClassName()+"::MetricValues" );

            ptic(tag);
            
            if( !M.InCacheQ(tag))
            {
                M.SetCache( tag, compute_metric(M) );
            }
            
            ptoc(tag);
            
            return M.template GetCache<ValueContainer_T>(tag);
        }
        
    
    protected:
        
        ValueContainer_T compute_metric( cref<MeshBase_T> M ) const
        {
            // Do a down cast and delegate implementation further to descendant class.
            cptr<Mesh_T> Q = dynamic_cast<const Mesh_T *>(&M);
                        
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
        virtual ValueContainer_T compute_metric( cref<Mesh_T> M ) const = 0;
        
        using Base_T::iter;
        using Base_T::rel_residuals;
        
        
    public:
        
        void MultiplyMetric(
            cref<MeshBase_T> M,
            cref<ExtReal> alpha, cptr<ExtReal> X,
            cref<ExtReal> beta,  mptr<ExtReal> Y,
            const Int  nrhs,
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const override
        {
            // Do a down cast and delegate implementation further to descendant class.
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                        
            if( Q != nullptr )
            {
                MultiplyMetric(*Q, alpha, X, beta, Y, nrhs, VF_flag, NF_flag, FF_flag );
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
            const Int  nrhs,
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const
        {
            ptic(ClassName()+"::MultiplyMetric");
            
            const auto & restrict S = M.GetBlockClusterTree().GetS();
            const auto & restrict T = M.GetBlockClusterTree().GetT();

            T.Pre( X, nrhs, op_type );
            
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
        
        virtual void MultiplyPreconditioner(
            cref<MeshBase_T> M, cptr<ExtReal> X, mptr<ExtReal> Y, const Int nrhs
        ) const override
        {
            // Do a down cast and delegate implementation further to descendant class.
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                        
            if( Q != nullptr )
            {
                MultiplyPreconditioner(*Q, X, Y, nrhs );
            }
            else
            {
                eprint(ClassName()+"::MultiplyPreconditioner: Input could not be downcast to compatible type. Doing nothing.");
            }
        }

        
        void MultiplyPreconditioner(
            cref<Mesh_T> M, cptr<ExtReal> X, mptr<ExtReal> Y, const Int nrhs
        ) const
        {
            ptic(ClassName()+"::MultiplyPreconditioner");

            multiply_preconditioner(M, X, Y, nrhs );
            
            ptoc(ClassName()+"::MultiplyPreconditioner");
        }
        
        
        virtual void multiply_preconditioner(
            cref<Mesh_T> M, cptr<ExtReal> X, mptr<ExtReal> Y, const Int nrhs
        ) const = 0;
        
        
        virtual void Solve(
            cref<MeshBase_T> M, cptr<ExtReal> B, mptr<ExtReal> X, const Int  nrhs,
            const Int  max_iter,
            const Real tolerance
        ) const override
        {
            // Do a down cast and delegate implementation further to descendant class.
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                        
            if( Q != nullptr )
            {
                Solve(*Q, B, X, nrhs, max_iter, tolerance );
            }
            else
            {
                eprint(ClassName()+"::Solve: Input could not be downcast to compatible type. Doing nothing.");
            }
        }
        
        void Solve(
            cref<Mesh_T> M, cptr<ExtReal> B, mptr<ExtReal> X, const Int nrhs,
            const Int  max_iter,
            const Real tolerance
        ) const
        {
            ptic(ClassName()+"::Solve");
            // The operator for the metric.
            auto A = [&M,nrhs,this]( cptr<ExtReal> X_, mptr<ExtReal> Y_ )
            {
                // Forward operator.
                MultiplyMetric( M, Scalar::One<ExtReal>, X_, Scalar::Zero<ExtReal>, Y_, nrhs );
            };

            // The operator for the preconditioner.
            auto P = [&M,nrhs,this]( cptr<ExtReal> X_, mptr<ExtReal> Y_ )
            {
                MultiplyPreconditioner( M, X_, Y_, nrhs );
            };
            
            
            if( nrhs == AMB_DIM )
            {
                CG_T CG ( M.VertexCount(), max_iter, nrhs, M.ThreadCount() );
                
                CG( A, P, B, nrhs, X, nrhs, tolerance );
                
                iter          = CG.IterationCount();
                rel_residuals = CG.RelativeResiduals();
            }
            else
            {
                CG_T CG ( M.VertexCount(), max_iter, nrhs, M.ThreadCount() );
                
                CG( A, P, B, nrhs, X, nrhs, tolerance );
                
                iter          = CG.IterationCount();
                rel_residuals = CG.RelativeResiduals();
            }
            
            ptoc(ClassName()+"::Solve");
        }
        
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
