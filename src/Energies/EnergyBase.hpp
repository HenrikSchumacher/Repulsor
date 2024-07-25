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
            std::string tag = ClassName()+"::Value";
            
            ptic(tag);
            
            if( !M.InCacheQ(tag) )
            {
                M.SetCache( tag, value(M) );
            }
            
            ptoc(tag);
            
            return M.template GetCache<ExtReal>(tag);
        }
        
    protected:
        
        // Actual implementation to be specified by descendants.
        virtual ExtReal value( cref<MeshBase_T> M ) const = 0;
        
    public:

        // Return the differential of the energy; use caching.
        cref<CotangentVector_T> Differential( mref<MeshBase_T> M ) const
        {
            std::string tag = ClassName()+"::Differential";
                
            ptic(tag);
            
            if( !M.InCacheQ(tag))
            {
                const ExtReal en = differential(M);
                
                if( !M.InCacheQ(ClassName()+"::Value"))
                {
                    M.SetCache( ClassName()+"::Value", en );
                }
                
                CotangentVector_T diff ( M.VertexCount(), M.AmbDim() );

                M.Assemble_ClusterTree_Derivatives( 
                    ExtReal(1), ExtReal(0), diff.data(), M.AmbDim()
                );
                
                M.SetCache( tag, std::move(diff) );
            }
            
            ptoc(tag);
            
            return M.template GetCache<CotangentVector_T>(tag);
        }
        
        // Return the differential of the energy to a pointer; don't use any caching.
        ExtReal Differential( 
            mref<MeshBase_T> M,
            const ExtReal alpha, const ExtReal beta, mptr<ExtReal> Y, const Int ldY
        ) const
        {
            ptic(ClassName()+"::Differential (pointer)");
            
            const ExtReal en = differential(M);
            
            M.Assemble_ClusterTree_Derivatives(
                alpha, beta, Y, ldY
            );
            
            ptoc(ClassName()+"::Differential (pointer)");
            
            return en;
        }
        
        // Return the differential of the energy to a pointer; don't use any caching.
        ExtReal Differential( 
            mref<MeshBase_T> M,
            mptr<ExtReal> Y, bool addto = false
        ) const
        {
            ptic(ClassName()+"::Differential (pointer)");
            
            const ExtReal en = differential(M);
            
            M.Assemble_ClusterTree_Derivatives(
                ExtReal(1), ExtReal(addto), Y, M.AmbDim()
            );
            
            ptoc(ClassName()+"::Differential (pointer)");
            
            return en;
        }
        
    protected:
        
        // Actual implementation to be specified by descendants.
        virtual ExtReal differential( const MeshBase_T & M ) const = 0;
        
        
    public:
        
        // Return the densities of the energy to a pointer; don't use any caching.
        ExtReal SimplexDensities( mref<MeshBase_T> M, mptr<ExtReal> density_ptr ) const
        {
            ptic(ClassName()+"::SimplexDensities (pointer)");
            
            const ExtReal en = density(M);
            
            M.Assemble_ClusterTree_SimplexDensities( density_ptr, ExtReal(1), false );
            
            ptoc(ClassName()+"::SimplexDensities (pointer)");
            
            return en;
        }
        
        // Return the densities of the energy to a pointer; don't use any caching.
        ExtReal VertexDensities( mref<MeshBase_T> M, mptr<ExtReal> density_ptr ) const
        {
            ptic(ClassName()+"::VertexDensities (pointer)");
            
            const ExtReal en = density(M);
            
            M.Assemble_ClusterTree_VertexDensities( density_ptr, ExtReal(1), false );
            
            ptoc(ClassName()+"::VertexDensities (pointer)");
            
            return en;
        }
        
    protected:
        
        // Actual implementation to be specified by descendants.
        virtual ExtReal density( const MeshBase_T & M ) const = 0;

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

