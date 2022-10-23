#pragma once

#define CLASS MetricBase

namespace Repulsor
{
    template<typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS
    {
        ASSERT_FLOAT(Real   );
        ASSERT_INT  (Int    );
        ASSERT_FLOAT(SReal  );
        ASSERT_FLOAT(ExtReal);
        
    public:
        
        using MeshBase_T        = SimplicialMeshBase<Real,Int,SReal,ExtReal>;
        using TangentVector_T   = typename MeshBase_T::TangentVector_T;
        using CotangentVector_T = typename MeshBase_T::CotangentVector_T;
        
        using Values_T          = Tensor2<Real,Int>;
        using ValueContainer_T  = std::array<Values_T,3>;
        
        CLASS() = default;

        virtual ~CLASS() = default;
        
    protected:
        
        static constexpr ExtReal ext_one = static_cast<ExtReal>(1);
        
    public:
        
        // Return the differential of the energy; use caching.
        ValueContainer_T & MetricValues( const MeshBase_T & M ) const
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
        
        // Actual implementation to be specified by descendants.
        virtual ValueContainer_T compute_metric( const MeshBase_T & M ) const = 0;
        
    public:

        // Return the differential of the energy; use caching.
        void MultiplyMetric(
            const MeshBase_T &             M,
            const ExtReal                  alpha,
            const ExtReal * restrict const X,
            const ExtReal                  beta,
                  ExtReal * restrict const Y,
            const Int                      rhs_count
        ) const
        {
            ptic(ClassName()+"::MultiplyMetric");
            
            auto & S = M.GetBlockClusterTree().GetS();
            auto & T = M.GetBlockClusterTree().GetT();

            T.Pre( X, rhs_count, KernelType::MixedOrder );
            
            S.RequireBuffers( T.BufferDimension() ); // Tell the S-side what it has to expect.
            
            multiply_metric(M,true,true,true);

            S.Post( Y, alpha, beta, KernelType::MixedOrder );
            
            ptoc(ClassName()+"::MultiplyMetric");
        }
        
        // Return the differential of the energy; use caching.
        void MultiplyMetric(
            const MeshBase_T & M,
            const ExtReal alpha,
            const TangentVector_T & X,
            const ExtReal beta,
                  CotangentVector_T & Y
        ) const
        {
            MultiplyMetric( M, alpha, X.data(), beta, Y.data(), X.Dimension(1) );
        }
        
//    protected:
//
//        virtual void multiply_metric( const MeshBase_T & M) = 0
        
    public:
        
        // Actual implementation to be specified by descendants.
        virtual void multiply_metric( const MeshBase_T & M, bool VF_flag, bool NF_flag, bool FF_flag ) const = 0;

        
    public:

        static std::string className()
        {
            return TO_STD_STRING(CLASS)+"<"+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
        virtual std::string ClassName() const
        {
            return className();
        }
        
    };

}// namespace Repulsor

#undef CLASS


