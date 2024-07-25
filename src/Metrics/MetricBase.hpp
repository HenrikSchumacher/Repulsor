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
        
        MetricBase() = default;

        virtual ~MetricBase() = default;
        
    protected:
        
        mutable Int iter = 0;
        mutable ConjugateGradient<VarSize,Real,Int>::RealVector_T rel_residuals;
        
    public:
        
        virtual mref<ValueContainer_T> MetricValues( cref<MeshBase_T> M ) const = 0;

        virtual void MultiplyMetric(
            cref<MeshBase_T> M,
            cref<ExtReal> alpha, cptr<ExtReal> X, const Int ldX,
            cref<ExtReal> beta,  mptr<ExtReal> Y, const Int ldY,
            const Int  rhs_count,
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const = 0;
        
        void MultiplyMetric(
            cref<MeshBase_T> M,
            cref<ExtReal> alpha,  cref<TangentVector_T>   X,
            cref<ExtReal> beta,   mref<CotangentVector_T> Y,
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
        
        virtual void MultiplyPreconditioner(
            cref<MeshBase_T> M, 
            cref<ExtReal> alpha, cptr<ExtReal> X, const Int ldX,
            cref<ExtReal> beta,  mptr<ExtReal> Y, const Int ldY,
            const Int rhs_count
        ) const = 0;

        virtual void Solve(
            cref<MeshBase_T> M, cptr<ExtReal> B, mptr<ExtReal> X, const Int  rhs_count,
            const Int  max_iter,
            const Real tolerance
        ) const = 0;
        
    
        ExtReal FrobeniusNorm( cref<MeshBase_T> M,
            cptr<ExtReal> X, const Int ldX,
            const Int rhs_count
        ) const
        {
            mref<Tensor2<Real,Int>> X_buf = M.XBuffer();
            mref<Tensor2<Real,Int>> Y_buf = M.YBuffer();
            
            X_buf.RequireSize( M.VertexCount(), rhs_count );
            Y_buf.RequireSize( M.VertexCount(), rhs_count );
            
            X_buf.Read( X, ldX );
            
            cptr<Real> X_ = X_buf.data();
            mptr<Real> Y_ = Y_buf.data();
            
            MultiplyMetric(M,
                Scalar::One <ExtReal>, X_, rhs_count,
                Scalar::Zero<ExtReal>, Y_, rhs_count,
                rhs_count
            );
            
            return Sqrt(
                dot_buffers( X_, Y_, M.VertexCount() * rhs_count, M.ThreadCount())
            );
        }
        
        Int CG_IterationCount() const
        {
            return iter;
        }
        
        ConjugateGradient<VarSize,Real,Int>::RealVector_T CG_RelativeResiduals() const
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
