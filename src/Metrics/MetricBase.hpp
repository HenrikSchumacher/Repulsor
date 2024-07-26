#pragma once

namespace Repulsor
{
    template<typename Mesh_T> class MetricBase {};
    
    template< typename Real, typename Int, typename LInt, typename SReal, typename ExtReal >
    class MetricBase<SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>>
    {
    public:
        
        using MeshBase_T             = SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>;
        using BlockClusterTreeBase_T = typename MeshBase_T::BlockClusterTreeBase_T;
        using TangentVector_T        = typename MeshBase_T::TangentVector_T;
        using CotangentVector_T      = typename MeshBase_T::CotangentVector_T;
        using ValueContainer_T       = MetricValueContainer<Real,LInt>;
        using CG_T                   = ConjugateGradient<VarSize,Real,Int>;
        
        MetricBase() = default;

        virtual ~MetricBase() = default;
        
    protected:
        
        mutable Int iter = 0;
        mutable CG_T::RealVector_T rel_residuals;
        
    public:
        
        virtual mref<ValueContainer_T> MetricValues( cref<MeshBase_T> M ) const = 0;

        // Computes Y = alpha * A.X + beta * Y
        virtual void MultiplyMetric(
            cref<MeshBase_T> M,
            const ExtReal alpha, cptr<ExtReal> X, const Int ldX,
            const ExtReal beta,  mptr<ExtReal> Y, const Int ldY,
            const Int  nrhs,
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const = 0;
        
        // Computes Y = alpha * A.X + beta * Y
        void MultiplyMetric(
            cref<MeshBase_T> M,
            const ExtReal alpha, cref<TangentVector_T>   X,
            const ExtReal beta,  mref<CotangentVector_T> Y,
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const
        {
            this->MultiplyMetric( M, 
                alpha, X.data(), X.Dimension(1),
                beta,  Y.data(), Y.Dimension(1),
                X.Dimension(1), VF_flag, NF_flag, FF_flag
            );
        }
        
        // Computes Y = alpha * P.X + beta * Y
        virtual void MultiplyPreconditioner(
            cref<MeshBase_T> M, 
            const ExtReal alpha, cptr<ExtReal> X, const Int ldX,
            const ExtReal beta,  mptr<ExtReal> Y, const Int ldY,
            const Int nrhs
        ) const = 0;

        
        // Computes X = alpha * A^{-1}.B + beta * X
        virtual void Solve(
            cref<MeshBase_T> M, 
            const ExtReal alpha, cptr<ExtReal> B, const Int ldB,
            const ExtReal beta,  mptr<ExtReal> X, const Int ldX,
            const Int  nrhs,
            const Int  max_iter,
            const Real tolerance
        ) const = 0;
        
    
        ExtReal MetricNorm( cref<MeshBase_T> M,
            cptr<ExtReal> X, const Int ldX, const Int nrhs
        ) const
        {
            // Not optimized; only meant for test purposes.
            mref<Tensor2<Real,Int>> X_buf = M.XBuffer( nrhs );
            mref<Tensor2<Real,Int>> Y_buf = M.YBuffer( nrhs );
            
            X_buf.Read( X, ldX );
            
            cptr<Real> X_ = X_buf.data();
            mptr<Real> Y_ = Y_buf.data();
            
            MultiplyMetric(M,
                Scalar::One <ExtReal>, X_, nrhs,
                Scalar::Zero<ExtReal>, Y_, nrhs,
                nrhs
            );
            
            return Sqrt( dot_buffers( X_, Y_, M.VertexCount() * nrhs, M.ThreadCount()) );
        }
        
        ExtReal PreconditionerNorm( cref<MeshBase_T> M,
            cptr<ExtReal> X, const Int ldX, const Int nrhs
        ) const
        {
            // Not optimized; only meant for test purposes.
            
            mref<Tensor2<Real,Int>> X_buf = M.XBuffer( nrhs );
            mref<Tensor2<Real,Int>> Y_buf = M.YBuffer( nrhs );
            
            X_buf.Read( X, ldX );
            
            cptr<Real> X_ = X_buf.data();
            mptr<Real> Y_ = Y_buf.data();
            
            MultiplyPreconditioner(M,
                Scalar::One <ExtReal>, X_, nrhs,
                Scalar::Zero<ExtReal>, Y_, nrhs,
                nrhs
            );
            
            return Sqrt( dot_buffers( X_, Y_, M.VertexCount() * nrhs, M.ThreadCount()) );
        }
        
        Int CG_IterationCount() const
        {
            return iter;
        }
        
        CG_T::RealVector_T CG_RelativeResiduals() const
        {
            return rel_residuals;
        }
        
    public:

        static std::string className()
        {
            return std::string("MetricBase<")+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }
        
        virtual std::string ClassName() const
        {
            return className();
        }
        
    }; // class MetricBase

}// namespace Repulsor
