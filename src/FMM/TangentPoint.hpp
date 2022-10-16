#pragma once

#define CLASS TangentPoint

namespace Tensors
{
    template< int S_DOM_DIM_, int T_DOM_DIM_, typename BlockClusterTree_T >
    class CLASS
    {
    public:
        
        using Real    = typename BlockClusterTree_T::Real;
        using Int     = typename BlockClusterTree_T::Int;
        using SReal   = typename BlockClusterTree_T::sReal;
        using ExtReal = typename BlockClusterTree_T::ExtReal;
        
        static constexpr Int S_DOM_DIM = S_DOM_DIM_;
        static constexpr Int T_DOM_DIM = T_DOM_DIM_;
        
//        using V_Kernel_T =

                
//        using F_Kernel_T =
        
        using N_Traversor_T = FMM_Tranversor<N_Kernel_T>;
        
        using SparsityPattern_T = SparsityPatternCSR<Int>;
        
    protected:
        
        CLASS() = default;
        
        CLASS( S_Mesh_T & S_M, T_Mesh_T & T_M )
        :   N_kernel( S_M->ClusterTree(), T_M->ClusterTree() )
        ,   N_traversor( )
        {
            
        }

        ~CLASS() = default;
        
    protected:
        
        N_Kernel_T    N_kernel;
        N_Traversor_T N_traversor;
        
    public:
        

//##############################################################################################
//      Matrix multiplication
//##############################################################################################
        
        Real Compute() const
        {
            TP_Kernel_NearField<S_DOM_DIM,T_DOM_DIM,ClusterTree_T,true,true,false,false> N_ker (
                S, T, alpha_half, beta_half
            );
            
            
        }
        
    public:
        
        std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"+kernel.ClassName()+">";
        }
        
    };
    
}// namespace Tensors

#undef CLASS

