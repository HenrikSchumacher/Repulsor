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

        using ValueContainer_T = Tensor2<ExtReal,Int>;
        using Differential_T   = Tensor2<ExtReal,Int>;
        using TangentVector_T  = Tensor2<ExtReal,Int>;
        
        explicit CLASS( ExtReal weight_ = static_cast<ExtReal>(1) )
        :
            weight(weight_)
        {}

        virtual ~CLASS() = default;
        
        virtual void Compute( const MeshBase_T & M ) const
        {
            ptic(ClassName()+"::Compute");
            
            // Create some dummies.
            std::array<ValueContainer_T,3> metric_values;
            std::array<ValueContainer_T,3> prec_values;
            
            M.GetClusterTree().CleanseDerivativeBuffers();
            
            std::any energy = compute(M);
            
            M.SetCache(ClassName()+"::Value", energy );
            
            std::any diff = std::make_any<Differential_T>( M.VertexCount(), M.AmbDim() );
            
            M.Assemble_ClusterTree_Derivatives(
                std::any_cast<Differential_T>(diff).data(),
                static_cast<ExtReal>(1),
                false
            );
            
            M.SetCache(ClassName()+"::Differential", diff);
            
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
                std::any energy = value(M);
                
                M.SetCache(ClassName()+"::Value", energy);
                
                return std::any_cast<ExtReal>(energy);
            }
            
            ptoc(ClassName()+"::Value");
            
            return std::any_cast<ExtReal>( M.GetCache(ClassName()+"::Value") );
        }
        
        // Actual implementation to be specified by descendants.
        virtual ExtReal value( const MeshBase_T & M ) const = 0;
        
        // Return the differential of the energy; use caching.
        const Differential_T & Differential( const MeshBase_T & M ) const
        {
            ptic((ClassName()+"::Differential"));
            if( !M.IsCached(ClassName()+"::Differential"))
            {
                M.GetClusterTree().CleanseDerivativeBuffers();
                
                differential(M);
                
                std::any diff = std::make_any<Differential_T>( M.VertexCount(), M.AmbDim() );
                
                M.Assemble_ClusterTree_Derivatives(
                    std::any_cast<Differential_T>(diff).data(),
                    static_cast<ExtReal>(1),
                    false
                );
                
                M.SetCache(ClassName()+"::Differential", diff);
            }
            
            ptoc((ClassName()+"::Differential"));
            
            return std::any_cast<Differential_T &>(
                M.GetCache(ClassName()+"::Differential")
            );
        }
        
        // Actual implementation to be specified by descendants.
        virtual void differential( const MeshBase_T & M ) const = 0;
        
//        virtual ExtReal Differential( const MESH_T & M, ExtReal * output, bool addTo = false ) const = 0;
        
//        virtual Differential_T & Metric( const MeshBase_T & M, const TangentVector_T & U ) const = 0;
        
//        virtual void SimplexEnergies( const MeshBase_T & M, ExtReal * output, bool addTo = false ) const = 0;
//        
//        virtual void Density( const MeshBase_T & M, ExtReal * output, bool addTo = false ) const = 0;
        
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

