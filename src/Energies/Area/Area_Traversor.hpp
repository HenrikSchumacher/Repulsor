#pragma once

#define CLASS Area_Traversor

namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename BlockClusterTree_T,
        bool energy_flag, bool diff_flag, bool hess_flag, bool metric_flag
    >
    class CLASS
    {
    public:
        
        using ClusterTree_T    = typename BlockClusterTree_T::ClusterTree_T;
        using Real             = typename BlockClusterTree_T::Real;
        using Int              = typename BlockClusterTree_T::Int;
        using SReal            = typename BlockClusterTree_T::SReal;
        using ExtReal          = typename BlockClusterTree_T::ExtReal;
        
        static constexpr Int S_DOM_DIM = S_DOM_DIM_;
        static constexpr Int T_DOM_DIM = T_DOM_DIM_;
        static constexpr bool is_symmetric = BlockClusterTree_T::IsSymmetric();
        
        using ValueContainer_T = Tensor2<Real,Int>;
        
        using Configurator_T = FMM_Configurator<ClusterTree_T>;
        
        
    public:
        
        CLASS() = delete;
        
        CLASS(
              const BlockClusterTree_T & bct_,
              std::array<ValueContainer_T,3> & metric_values_,
              std::array<ValueContainer_T,3> & prec_values_,
              const Real q_,
              const Real p_
        )
        :   bct(bct_)
        ,   metric_values ( metric_values_ )
        ,   prec_values   ( prec_values_   )
        {}

        ~CLASS() = default;
        

    protected:

        const BlockClusterTree_T & bct;
        
        std::array<ValueContainer_T,3> & metric_values;
        std::array<ValueContainer_T,3> & prec_values;
        
        Real en = 0;
        
    public:
        

//##############################################################################################
//      Compute
//##############################################################################################
        
        Real Compute()
        {
            en = static_cast<Real>(0);

            VeryNearField();
            NearField    ();
            FarField     ();
            
            return en;
        }
        
    protected:
        
        void VeryNearField()
        {
            ptic(ClassName()+"::VeryNearField");
            
                ptic(ClassName()+"::VeryNearField prepare kernels");
                
                Configurator_T conf (
                    bct.GetS(), bct.GetT(),
                    metric_values[FMM_Type::VF], prec_values[FMM_Type::VF]
                );
                
                using Kernel_T = Area_Kernel_VF<
                    S_DOM_DIM, T_DOM_DIM,
                    ClusterTree_T,
                    BlockClusterTree_T::IsSymmetric(),
                    energy_flag, diff_flag, hess_flag, metric_flag
                >;
                
                Kernel_T ker ( conf, bct.NearFieldSeparationParameter(), 20 );
                
                ptoc(ClassName()+"::VeryNearField prepare kernels");
                
                ptic(ClassName()+"::VeryNearField prepare traversor");
                FMM_Traversor<Kernel_T> traversor ( bct.VeryNear(), ker );
                ptoc(ClassName()+"::VeryNearField prepare traversor");
                
                ptic(ClassName()+"::VeryNearField traversal");
                en += traversor.Compute();
                ptoc(ClassName()+"::VeryNearField traversal");
            
            ptoc(ClassName()+"::VeryNearField");
        }
            
        void NearField()
        {
            ptic(ClassName()+"::NearField");
            
                ptic(ClassName()+"::NearField prepare kernels");
                
                Configurator_T conf (
                    bct.GetS(), bct.GetT(),
                    metric_values[FMM_Type::NF], prec_values[FMM_Type::NF]
                );
        
                using Kernel_T = Area_Kernel_NF<
                    S_DOM_DIM, T_DOM_DIM,
                    ClusterTree_T,
                    BlockClusterTree_T::IsSymmetric(),
                    energy_flag, diff_flag, hess_flag, metric_flag
                >;
                
                Kernel_T ker ( conf );
                
                ptoc(ClassName()+"::NearField prepare kernels");
                
                ptic(ClassName()+"::NearField prepare traversor");
                FMM_Traversor<Kernel_T> traversor ( bct.Near(), ker );
                ptoc(ClassName()+"::NearField prepare traversor");
                
                ptic(ClassName()+"::NearField traversal");
                en += traversor.Compute();
                ptoc(ClassName()+"::NearField traversal");
            
            ptoc(ClassName()+"::NearField");
        }
            
        void FarField()
        {
            ptic(ClassName()+"::FarField");

                ptic(ClassName()+"::FarField prepare kernels");
                
                Configurator_T conf (
                    bct.GetS(), bct.GetT(),
                    metric_values[FMM_Type::FF], prec_values[FMM_Type::FF]
                );
                
                using Kernel_T = Area_Kernel_FF<
                    ClusterTree_T,
                    BlockClusterTree_T::IsSymmetric(),
                    energy_flag, diff_flag, hess_flag, metric_flag
                >;

                Kernel_T ker ( conf );
                
                ptoc(ClassName()+"::FarField prepare kernels");
                
                ptic(ClassName()+"::FarField prepare traversor");
                FMM_Traversor<Kernel_T> traversor ( bct.Far(), ker );
                ptoc(ClassName()+"::FarField prepare traversor");
                
                ptic(ClassName()+"::FarField traversal");
                en += traversor.Compute();
                ptoc(ClassName()+"::FarField traversal");
            
            ptoc(ClassName()+"::FarField");
        }
        
    public:
        
        std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(S_DOM_DIM)+","+ToString(T_DOM_DIM)+","+bct.ClassName()+">";
        }
    };
    
}// namespace Repulsor

#undef CLASS

