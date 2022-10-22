#pragma once

#define CLASS EnergyBase

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
        
        using MeshBase_T = SimplicialMeshBase<Real,Int,SReal,ExtReal>;

        using TangentVector_T   = Tensor2<ExtReal,Int>;
        using CotangentVector_T = Tensor2<ExtReal,Int>;
        using Values_T          = Tensor2<   Real,Int>;
        using ValueContainer_T  = std::array<Values_T,3>;
        
        explicit CLASS( ExtReal weight_ = static_cast<ExtReal>(1) )
        :
            weight(weight_)
        {}

        virtual ~CLASS() = default;
        
        virtual void Compute( const MeshBase_T & M ) const
        {
            ptic(ClassName()+"::Compute");
            
            // Create some dummies.
            ValueContainer_T metric_values;
            ValueContainer_T prec_values;
            
            std::any energy = weight * compute(M);
            
            M.SetCache( ClassName()+"::Value", energy );
            
            CotangentVector_T diff ( M.VertexCount(), M.AmbDim() );

            M.Assemble_ClusterTree_Derivatives( diff.data(), weight, false );
            
            std::any thing = diff;
            
            M.SetCache( ClassName()+"::Differential", thing );
            
            ptoc(ClassName()+"::Compute");
        }
        
        // Actual implementation to be specified by descendants.
        virtual ExtReal compute( const MeshBase_T & M ) const = 0;
        
        // Return the value of the energy; use caching.
        ExtReal Value( const MeshBase_T & M ) const
        {
            ptic(ClassName()+"::Value");
            
            if( !M.IsCached(ClassName()+"::Value"))
            {
                std::any energy = weight * value(M);
                
                M.SetCache( ClassName()+"::Value", energy );
            }
            
            ptoc(ClassName()+"::Value");
            
            return std::any_cast<ExtReal>( M.GetCache(ClassName()+"::Value") );
        }
        
        // Actual implementation to be specified by descendants.
        virtual ExtReal value( const MeshBase_T & M ) const = 0;
        

        // Return the differential of the energy; use caching.
        const CotangentVector_T & Differential( const MeshBase_T & M ) const
        {
            ptic(ClassName()+"::Differential");
            
            if( !M.IsCached(ClassName()+"::Differential"))
            {
//                M.GetS().CleanseDerivativeBuffers();

                differential(M);
                
                CotangentVector_T diff ( M.VertexCount(), M.AmbDim() );

                M.Assemble_ClusterTree_Derivatives( diff.data(), weight, false );
                
                // TODO: Find out whether this incurs a copy operation.
                std::any thing = diff;
                
                M.SetCache( ClassName()+"::Differential", thing );
            }
            
            ptoc(ClassName()+"::Differential");
            
            return std::any_cast<CotangentVector_T &>(
                M.GetCache(ClassName()+"::Differential")
            );
        }
        
        // Actual implementation to be specified by descendants.
        virtual void differential( const MeshBase_T & M ) const = 0;
        
//        virtual void SimplexEnergies( const MeshBase_T & M, ExtReal * output, bool addTo = false ) const = 0;
//        
//        virtual void Density( const MeshBase_T & M, ExtReal * output, bool addTo = false ) const = 0;
        
        
        // Return the differential of the energy; use caching.
        ValueContainer_T & MetricValues( const MeshBase_T & M ) const
        {
            ptic(ClassName()+"::MetricValues");
            
            if( !M.IsCached(ClassName()+"::MetricValues"))
            {
                std::any thing = compute_metric(M);
                
                M.SetCache( ClassName()+"::MetricValues", thing );
            }
            
            ptoc(ClassName()+"::MetricValues");
            
            return std::any_cast<ValueContainer_T &>(
                M.GetCache(ClassName()+"::MetricValues")
            );
        }
        
        // Actual implementation to be specified by descendants.
        virtual ValueContainer_T compute_metric( const MeshBase_T & M ) const = 0;
        

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
            
            auto & S = M.GetClusterTree();
            auto & T = M.GetClusterTree();
            
            T.Pre( X, rhs_count, KernelType::MixedOrder );
            
            S.RequireBuffers( T.BufferDimension() ); // Tell the S-side what it has to expect.
            
            multiply_metric(M);

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
        
        // Actual implementation to be specified by descendants.
        virtual void multiply_metric( const MeshBase_T & M ) const = 0;
        
        ExtReal GetWeight()  const
        {
            return weight;
        }
        
        void SetWeight( const ExtReal weight_ )
        {
            weight = weight_;
        }

        
    protected:
        
        ExtReal weight = static_cast<ExtReal>(1);
        
    public:
        
//        virtual std::string Stats() const = 0;

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

