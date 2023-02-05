#pragma once

namespace Repulsor
{
    template<typename Real, typename Int, typename SReal, typename ExtReal>
    class EnergyBase
    {
        ASSERT_FLOAT( Real    );
        ASSERT_FLOAT( SReal   );
        ASSERT_FLOAT( ExtReal );
        ASSERT_INT  ( Int     );
        
    public:
        
        
        using MeshBase_T         = SimplicialMeshBase<Real,Int,SReal,ExtReal>;
        using BlockClusterTree_T = typename MeshBase_T::BlockClusterTree_T;
        using LInt               = size_t;
        using TangentVector_T    = typename MeshBase_T::TangentVector_T;
        using CotangentVector_T  = typename MeshBase_T::CotangentVector_T;
        
        using Values_T           = Tensor2<Real,LInt>;
        using ValueContainer_T   = std::unordered_map<std::string,Values_T>;
        
        EnergyBase() = default;

        virtual ~EnergyBase() = default;
        
    protected:

        
    public:
        
//        virtual void SimplexEnergies( const MeshBase_T & M, ExtReal * output, bool addTo = false ) const = 0;
//
//        virtual void Density( const MeshBase_T & M, ExtReal * output, bool addTo = false ) const = 0;
        
        
    public:
        
        // Return the value of the energy; use caching.
        ExtReal Value( const MeshBase_T & M ) const
        {
            ptic(ClassName()+"::Value");
            
            if( !M.IsCached(ClassName()+"::Value"))
            {
                std::any energy = value(M);
                
                M.SetCache( ClassName()+"::Value", energy );
            }
            
            ptoc(ClassName()+"::Value");
            
            return std::any_cast<ExtReal>( M.GetCache(ClassName()+"::Value") );
        }
        
    protected:
        
        // Actual implementation to be specified by descendants.
        virtual ExtReal value( const MeshBase_T & M ) const = 0;
        
    public:

        // Return the differential of the energy; use caching.
        const CotangentVector_T & Differential( const MeshBase_T & M ) const
        {
            ptic(ClassName()+"::Differential");
            
            if( !M.IsCached(ClassName()+"::Differential"))
            {
//                M.GetS().CleanseDerivativeBuffers();

                differential(M);
                
                CotangentVector_T diff ( M.VertexCount(), M.AmbDim() );

                M.Assemble_ClusterTree_Derivatives( diff.data(), ExtReal(1), false );
                
                std::any thing ( std::move(diff) );
                // TODO: Find out whether this incurs a copy operation.
                M.SetCache( ClassName()+"::Differential", thing );
            }
            
            ptoc(ClassName()+"::Differential");
            
            return std::any_cast<CotangentVector_T &>(
                M.GetCache(ClassName()+"::Differential")
            );
        }
        
    protected:
        
        // Actual implementation to be specified by descendants.
        virtual void differential( const MeshBase_T & M ) const = 0;

    public:

        static std::string className()
        {
            return std::string("EnergyBase<")+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }
        
        virtual std::string ClassName() const
        {
            return className();
        }
        
    }; // class EnergyBase

}// namespace Repulsor

