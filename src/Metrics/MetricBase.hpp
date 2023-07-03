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
        
    public:
        
        virtual ValueContainer_T & MetricValues( const MeshBase_T & restrict M ) const = 0;

        virtual void MultiplyMetric(
            const MeshBase_T & restrict M,
            const ExtReal alpha, ptr<ExtReal> X,
            const ExtReal beta,  mut<ExtReal> Y,
            const Int  rhs_count,
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const = 0;
        
        void MultiplyMetric(
            const MeshBase_T & restrict M,
            const ExtReal alpha,
            const TangentVector_T & restrict X,
            const ExtReal beta,
                  CotangentVector_T & restrict Y,
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const
        {
            this->MultiplyMetric( M, alpha, X.data(), beta, Y.data(), X.Dimension(1), VF_flag, NF_flag, FF_flag );
        }
        
//        virtual void MultiplyPreconditioner(
//            const MeshBase_T & restrict M,
//            const ExtReal alpha, ptr<ExtReal> X,
//            const ExtReal beta,  mut<ExtReal> Y,
//            const Int  rhs_count
//        ) const = 0;
//
//        virtual void SolveMetric(
//            const MeshBase_T & restrict M,
//            const ExtReal alpha, ptr<ExtReal> X,
//            const ExtReal beta,  mut<ExtReal> Y,
//            const Int  rhs_count,
//            const Int  max_iter,
//            const Real tolerance
//        ) const = 0;
        
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
