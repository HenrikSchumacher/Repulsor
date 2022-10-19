#pragma once

#include "TangentPoint/TP_Kernel_FF.hpp"
#include "TangentPoint/TP_Kernel_NF.hpp"
#include "TangentPoint/TP_Kernel_VF.hpp"
#include "TangentPoint/TP_Traversor.hpp"


#define CLASS TangentPoint

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS
    {
    public:
        
        using Mesh_T                  = SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,SReal,ExtReal>;
        
        using ValueContainer_T        = Tensor2<Real,Int>;
        
        using DifferentialContainer_T = Tensor2<ExtReal,Int>;
        
        using BlockClusterTree_T      = typename Mesh_T::BlockClusterTree_T;
        
        template<typename T_q, typename T_p>
        CLASS( const T_q q_, const T_p p_ )
        :   q( static_cast<Real>(q_) )
        ,   p( static_cast<Real>(p_) )
        {}
        
    protected:
        
        const Real q;
        const Real p;
        
    public:
        
        void Require( const Mesh_T & M )
        {
            // Create some dummies.
            std::array<ValueContainer_T,3> metric_values;
            std::array<ValueContainer_T,3> prec_values;
            
            M.GetClusterTree().CleanseDerivativeBuffers();
            
            TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,false,true,false,false>
                traversor( M.GetBlockClusterTree(), metric_values, prec_values, p, q );
            
            std::any energy = traversor.Compute();
            
            M.SetCache("TangentPoint_Energy", energy );
            
            DifferentialContainer_T diff ( M.VertexCount(), AMB_DIM );
            
            M.Assemble_ClusterTree_Derivatives( diff.data(), static_cast<ExtReal>(1), false );
            
            std::any result = std::move(diff);
            
            M.SetCache("TangentPoint_DEnergy", result );
        }
        
        ExtReal Energy( const Mesh_T & M )
        {
            if( M.IsCached("TangentPoint_Energy"))
            {
                return std::any_cast<ExtReal>( M.GetCache("TangentPoint_Energy") );
            }
            else
            {
                // Create some dummies.
                std::array<ValueContainer_T,3> metric_values;
                std::array<ValueContainer_T,3> prec_values;
                
                TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,true,false,false,false>
                    traversor( M.GetBlockClusterTree(), metric_values, prec_values, p, q );
                
                std::any energy = traversor.Compute();
                
                M.SetCache("TangentPoint_Energy", energy );
                
                return std::any_cast<ExtReal>(energy);
            }
        }
        
        const DifferentialContainer_T & DEnergy( const Mesh_T & M )
        {
            if( !M.IsCached("TangentPoint_DEnergy"))
            {
                // Create some dummies.
                std::array<ValueContainer_T,3> metric_values;
                std::array<ValueContainer_T,3> prec_values;
                
                M.GetClusterTree().CleanseDerivativeBuffers();
                
                TP_Traversor<DOM_DIM,DOM_DIM,BlockClusterTree_T,false,true,false,false>
                    traversor( M.GetBlockClusterTree(), metric_values, prec_values, p, q );
                
                (void)traversor.Compute();
                
                DifferentialContainer_T diff ( M.VertexCount(), AMB_DIM );
                
                M.Assemble_ClusterTree_Derivatives( diff.data(), static_cast<ExtReal>(1), false );
                
                std::any result = std::move(diff);
                
                M.SetCache("TangentPoint_DEnergy", result );
            }
            
            return std::any_cast<DifferentialContainer_T>( M.GetCache("TangentPoint_DEnergy") );
        }
        
        
        std::string ClassName() const
        {
            return className();
        }
        
    private:
        
        static std::string className()
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    };
    
}

#undef CLASS
