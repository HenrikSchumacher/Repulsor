#pragma once

namespace Repulsor
{
    template<typename MeshBase_T> class EnergyBase {};

    template< typename Real, typename Int, typename LInt, typename SReal, typename ExtReal >
    class EnergyBase<SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>>
    {
        
    public:
        
        
        using MeshBase_T             = SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>;
        using BlockClusterTreeBase_T = typename MeshBase_T::BlockClusterTreeBase_T;
        using TangentVector_T        = typename MeshBase_T::TangentVector_T;
        using CotangentVector_T      = typename MeshBase_T::CotangentVector_T;
        using ValueContainer_T       = MetricValueContainer<Real,LInt>;
        
        EnergyBase()
        {}

        virtual ~EnergyBase() = default;
        
    protected:

        // We need some dummy values so that kernels won't complain, in particular when run in parallel.
        mutable ValueContainer_T metric_values;
        
    public:
        
//        virtual void SimplexEnergies( const MeshBase_T & M, ExtReal * output, bool addTo = false ) const = 0;
//
//        virtual void Density( const MeshBase_T & M, ExtReal * output, bool addTo = false ) const = 0;
        
        
    public:
        
        // Return the value of the energy; use caching.
        ExtReal Value( mref<MeshBase_T> M ) const
        {
            ptic(ClassName()+"::Value");
            
            if( !M.InCacheQ(ClassName()+"::Value"))
            {
//                std::any energy = value(M);
                
                M.SetCache( ClassName()+"::Value", value(M) );
            }
            
            ptoc(ClassName()+"::Value");
            
            return std::any_cast<ExtReal>( M.GetCache(ClassName()+"::Value") );
        }
        
    protected:
        
        // Actual implementation to be specified by descendants.
        virtual ExtReal value( cref<MeshBase_T> M ) const = 0;
        
    public:

        // Return the differential of the energy; use caching.
        cref<CotangentVector_T> Differential( mref<MeshBase_T> M ) const
        {
            ptic(ClassName()+"::Differential");
            
            if( !M.InCacheQ(ClassName()+"::Differential"))
            {
//                M.GetS().CleanseDerivativeBuffers();

                const ExtReal en = differential(M);
                
                if( !M.InCacheQ(ClassName()+"::Value"))
                {
                    M.SetCache( ClassName()+"::Value", en );
                }
                
                CotangentVector_T diff ( M.VertexCount(), M.AmbDim() );

                M.Assemble_ClusterTree_Derivatives( diff.data(), ExtReal(1), false );
                
                M.SetCache( ClassName()+"::Differential", std::move(diff) );
            }
            
            ptoc(ClassName()+"::Differential");
            
            return std::any_cast<CotangentVector_T &>(
                M.GetCache(ClassName()+"::Differential")
            );
        }
        
        // Return the differential of the energy to a pointer; don't use any caching.
        ExtReal Differential( mref<MeshBase_T> M, mptr<ExtReal> diff ) const
        {
            ptic(ClassName()+"::Differential (pointer)");
            
            const ExtReal en = differential(M);
            
            M.Assemble_ClusterTree_Derivatives( diff, ExtReal(1), false );
            
            ptoc(ClassName()+"::Differential (pointer)");
            
            return en;
        }
        
    protected:
        
        // Actual implementation to be specified by descendants.
        virtual ExtReal differential( const MeshBase_T & M ) const = 0;

    public:

        static std::string className()
        {
            return std::string("EnergyBase<")+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }
        
        virtual std::string ClassName() const
        {
            return className();
        }
        
    }; // class EnergyBase

}// namespace Repulsor

