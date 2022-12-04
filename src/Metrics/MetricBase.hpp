#pragma once

namespace Repulsor
{
    template<typename Real, typename Int, typename SReal, typename ExtReal>
    class MetricBase
    {
        ASSERT_FLOAT(Real   );
        ASSERT_INT  (Int    );
        ASSERT_FLOAT(SReal  );
        ASSERT_FLOAT(ExtReal);
        
    public:
        
        using MeshBase_T         = SimplicialMeshBase<Real,Int,SReal,ExtReal>;
        
        using TangentVector_T    = typename MeshBase_T::TangentVector_T;
        using CotangentVector_T  = typename MeshBase_T::CotangentVector_T;
        
        using BlockClusterTree_T = typename MeshBase_T::BlockClusterTree_T;
        
        using Values_T           = Tensor2<Real,size_t>;
        using ValueContainer_T   = std::unordered_map<std::string,Values_T>;
        
        MetricBase() = default;

        virtual ~MetricBase() = default;
        
    protected:
        
        static constexpr ExtReal ext_one = static_cast<ExtReal>(1);
        
    public:
        
        virtual ValueContainer_T & MetricValues( const MeshBase_T & M ) const = 0;

        virtual void MultiplyMetric(
            const MeshBase_T &             M,
            const ExtReal                  alpha,
            const ExtReal * restrict const X,
            const ExtReal                  beta,
                  ExtReal * restrict const Y,
            const Int                      rhs_count,
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const = 0;
        
        // Return the differential of the energy; use caching.
        void MultiplyMetric(
            const MeshBase_T & M,
            const ExtReal alpha,
            const TangentVector_T & X,
            const ExtReal beta,
                  CotangentVector_T & Y,
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const
        {
            this->MultiplyMetric( M, alpha, X.data(), beta, Y.data(), X.Dimension(1), VF_flag, NF_flag, FF_flag );
        }

        
    public:

        static std::string className()
        {
            return "MetricBase<"+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
        virtual std::string ClassName() const
        {
            return className();
        }
        
    }; // class MetricBase

}// namespace Repulsor
