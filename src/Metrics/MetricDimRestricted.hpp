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
        
        mref<ValueContainer_T> MetricValues( cref<MeshBase_T> M ) const override
        {
            std::string tag ( ClassName()+"::MetricValues" );

            TOOLS_PTIC(tag);
            if( !M.InCacheQ(tag))
            {
                M.SetCache( tag, compute_metric(M) );
            }
            TOOLS_PTOC(tag);
            
            return M.template GetCache<ValueContainer_T>(tag);
        }
        
        mref<ValueContainer_T> MetricValues( cref<Mesh_T> M ) const
        {
            std::string tag ( ClassName()+"::MetricValues" );

            TOOLS_PTIC(tag);
            
            if( !M.InCacheQ(tag))
            {
                M.SetCache( tag, compute_metric(M) );
            }
            
            TOOLS_PTOC(tag);
            
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
            const ExtReal alpha, cptr<ExtReal> X, const Int ldX,
            const ExtReal beta,  mptr<ExtReal> Y, const Int ldY,
            const Int nrhs,
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const override
        {
            // Do a down cast and delegate implementation further to descendant class.
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                        
            if( Q != nullptr )
            {
                MultiplyMetric(*Q,
                    alpha, X, ldX,
                    beta,  Y, ldY,
                    nrhs, VF_flag, NF_flag, FF_flag
                );
            }
            else
            {
                eprint(ClassName()+"::MultiplyMetric: Input could not be downcast to compatible type. Doing nothing.");
            }
        }
        
        void MultiplyMetric(
            cref<Mesh_T> M,
            const ExtReal alpha, cptr<ExtReal> X, const Int ldX,
            const ExtReal beta,  mptr<ExtReal> Y, const Int ldY,
            const Int  nrhs,
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const
        {
            std::string tag = ClassName()+"::MultiplyMetric";
            
            TOOLS_PTIC(tag);
            
            if( nrhs != AMB_DIM )
            {
                NrhsError( tag, nrhs );
                
                TOOLS_PTOC(tag);
                
                return;
            }
            
            const auto & restrict S = M.GetBlockClusterTree().GetS();
            const auto & restrict T = M.GetBlockClusterTree().GetT();

            T.Pre( X, ldX, nrhs, op_type );
            
            S.RequireBuffers( T.BufferDim() ); // Tell the S-side what it has to expect.
            
            multiply_metric( M, VF_flag, NF_flag, FF_flag );

            S.Post( alpha, beta, Y, ldY, op_type );
            
            TOOLS_PTOC(tag);
        }
    
        
        
        // Actual implementation to be specified by descendants.
        virtual void multiply_metric(
            cref<Mesh_T> M,
            const bool VF_flag, const bool NF_flag, const bool FF_flag
        ) const = 0;

        
    public:
        
        virtual void MultiplyPreconditioner(
            cref<MeshBase_T> M, 
            const ExtReal alpha, cptr<ExtReal> X, const Int ldX,
            const ExtReal beta,  mptr<ExtReal> Y, const Int ldY,
            const Int nrhs
        ) const override
        {
            // Do a down cast and delegate implementation further to descendant class.
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                        
            if( Q != nullptr )
            {
                MultiplyPreconditioner(*Q, alpha, X, ldX, beta, Y, ldY, nrhs );
            }
            else
            {
                eprint(ClassName()+"::MultiplyPreconditioner: Input could not be downcast to compatible type. Doing nothing.");
            }
        }

        
        void MultiplyPreconditioner(
            cref<Mesh_T> M,
            const ExtReal alpha, cptr<ExtReal> X, const Int ldX,
            const ExtReal beta,  mptr<ExtReal> Y, const Int ldY,
            const Int nrhs
        ) const
        {
            std::string tag = ClassName()+"::MultiplyPreconditioner";
            
            TOOLS_PTIC(tag);

            if( nrhs != AMB_DIM )
            {
                NrhsError( tag, nrhs );
                
                TOOLS_PTOC(tag);
                
                return;
            }
            
            multiply_preconditioner(M,
                alpha, X, ldX,
                beta,  Y, ldY,
                nrhs
            );
            
            TOOLS_PTOC(tag);
        }
        
        
        virtual void multiply_preconditioner(
            cref<Mesh_T> M, 
            const ExtReal alpha, cptr<ExtReal> X, const Int ldX,
            const ExtReal beta,  mptr<ExtReal> Y, const Int ldY,
            const Int nrhs
        ) const = 0;
        
        
        virtual void Solve(
            cref<MeshBase_T> M, 
            const ExtReal alpha, cptr<ExtReal> B, const Int ldB,
            const ExtReal beta,  mptr<ExtReal> X, const Int ldX,
            const Int  nrhs,
            const Int  max_iter,
            const Real tolerance
        ) const override
        {
            // Do a down cast and delegate implementation further to descendant class.
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                        
            if( Q != nullptr )
            {
                Solve(*Q, 
                    alpha, B, ldB,
                    beta,  X, ldX,
                    nrhs, max_iter, tolerance
                );
            }
            else
            {
                eprint(ClassName()+"::Solve: Input could not be downcast to compatible type. Doing nothing.");
            }
        }
        
        void Solve(
            cref<Mesh_T> M, 
            const ExtReal alpha, cptr<ExtReal> B, const Int ldB,
            const ExtReal beta,  mptr<ExtReal> X, const Int ldX,
            const Int  nrhs,
            const Int  max_iter,
            const Real tolerance
        ) const
        {
            std::string tag = ClassName()+"::Solve";
            
            TOOLS_PTIC(tag);

            if( nrhs != AMB_DIM )
            {
                NrhsError( tag, nrhs );
                
                TOOLS_PTOC(tag);
                
                return;
            }
            
            auto A = [&M,nrhs,this]( cptr<Real> X_, mptr<Real> Y_ )
            {
                // Forward operator.
                MultiplyMetric( M, 
                    Scalar::One <Real>, X_, nrhs,
                    Scalar::Zero<Real>, Y_, nrhs,
                    nrhs
                );
            };

            // The operator for the preconditioner.
            auto P = [&M,nrhs,this]( cptr<Real> X_, mptr<Real> Y_ )
            {
                MultiplyPreconditioner( M, 
                   Scalar::One <Real>, X_, nrhs,
                   Scalar::Zero<Real>, Y_, nrhs,
                   nrhs
                );
            };
            
            ConjugateGradient<AMB_DIM,Real,Int,false,false> solver (
                M.VertexCount(), max_iter, nrhs, M.ThreadCount()
            );
            
            solver( A, P, alpha, B, ldB, beta, X, ldX, tolerance );
            

//            GMRES<AMB_DIM,Real,Int,Side::Left,false,false> solver (
//                M.VertexCount(), max_iter, nrhs, M.ThreadCount()
//            );
//            
//            solver( A, P, alpha, B, ldB, beta, X, ldX, tolerance, 10 );
            
            iter          = solver.IterationCount();
            rel_residuals = solver.RelativeResiduals();
            
            TOOLS_PTOC(tag);
        }
        
        
    private:
        
        void NrhsError( const std::string & tag, Int nrhs ) const
        {
            eprint( tag + ": nrhs = " + ToString(nrhs)+ " != " + ToString(AMB_DIM) + " = AMB_DIM. The current implementation only accepts nrhs = AMB_DIM. (After all, this is a metric on the space of infinitesimal displacement in " + ToString(AMB_DIM) + "-dimensional Euclidean space.) Doing nothing. The argument nrhs is there only for compatibility reasons. Please set nrhs = " + ToString(AMB_DIM) + ".");
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
