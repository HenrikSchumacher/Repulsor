#pragma once

#define CLASS TP_Traversor

namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename BlockClusterTree_T,
        bool energy_flag, bool diff_flag, bool hess_flag, bool metric_flag, bool prec_flag
    >
    class CLASS
    {
    public:
        
        using ClusterTree_T    = typename BlockClusterTree_T::ClusterTree_T;
        using Real             = typename BlockClusterTree_T::Real;
        using Int              = typename BlockClusterTree_T::Int;
        using SReal            = typename BlockClusterTree_T::SReal;
        using ExtReal          = typename BlockClusterTree_T::ExtReal;
        
        static constexpr Int AMB_DIM   = BlockClusterTree_T::AMB_DIM;
        static constexpr Int S_DOM_DIM = S_DOM_DIM_;
        static constexpr Int T_DOM_DIM = T_DOM_DIM_;
        static constexpr bool is_symmetric = BlockClusterTree_T::IsSymmetric();
        
        using Values_T          = Tensor2<Real,Int>;
        using ValueContainer_T  = std::array<Values_T,3>;
        using TangentVector_T   = Tensor2<ExtReal,Int>;
        using CotangentVector_T = Tensor2<ExtReal,Int>;
        
        using Configurator_T   = FMM_Configurator<ClusterTree_T>;
        
        
    public:
        
        CLASS() = delete;
        
        CLASS(
              const BlockClusterTree_T & bct_,
              ValueContainer_T & metric_values_,
              ValueContainer_T & prec_values_,
              const Real q_,
              const Real p_
        )
        :   bct(bct_)
        ,   q             ( q_ )
        ,   p             ( p_ )
        ,   q_half_real   ( q / static_cast<Real>(2) )
        ,   p_half_real   ( p / static_cast<Real>(2) )
        ,   q_half_int    ( q_half_real )
        ,   p_half_int    ( p_half_real )
        ,   q_half_is_int ( q_half_real == static_cast<Real>(q_half_int) )
        ,   p_half_is_int ( p_half_real == static_cast<Real>(p_half_int) )
        ,   metric_values ( metric_values_ )
        ,   prec_values   ( prec_values_   )
        {}

        ~CLASS() = default;
        

    protected:

        const BlockClusterTree_T & bct;
        
        const Real q;
        const Real p;
        const Real q_half_real;
        const Real p_half_real;
        
        const Int q_half_int;
        const Int p_half_int;
        
        const bool q_half_is_int;
        const bool p_half_is_int;
        
        ValueContainer_T & metric_values;
        ValueContainer_T & prec_values;
        
        Real en = 0;
        
    public:
        

//##############################################################################################
//      Compute
//##############################################################################################
        
        Real Compute()
        {
            en = static_cast<Real>(0);

            if( q_half_is_int )
            {
                if( p_half_is_int )
                {
                    Compute_VF<Int,Int>( q_half_int, p_half_int );
                    Compute_NF<Int,Int>( q_half_int, p_half_int );
                    Compute_FF<Int,Int>( q_half_int, p_half_int );
                }
                else
                {
                    Compute_VF<Int,Real>( q_half_int, p_half_real );
                    Compute_NF<Int,Real>( q_half_int, p_half_real );
                    Compute_FF<Int,Real>( q_half_int, p_half_real );
                }
            }
            else
            {
                if( p_half_is_int)
                {
                    Compute_VF<Real,Int>( q_half_real, p_half_int );
                    Compute_NF<Real,Int>( q_half_real, p_half_int );
                    Compute_FF<Real,Int>( q_half_real, p_half_int );
                }
                else
                {
                    Compute_VF<Real,Real>( q_half_real, p_half_real );
                    Compute_NF<Real,Real>( q_half_real, p_half_real );
                    Compute_FF<Real,Real>( q_half_real, p_half_real );
                }
            }
            
            return en;
        }
        
    protected:
        
        template< typename T1, typename T2 >
        void Compute_VF( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::Compute_VF");
            
                ptic(ClassName()+"::Compute_VF prepare kernels");
                
                Configurator_T conf (
                    bct.GetS(), bct.GetT(),
                    metric_values[FMM::Type::VF], prec_values[FMM::Type::VF]
                );
                
                using Kernel_T = TP_Kernel_VF<
                    S_DOM_DIM, T_DOM_DIM,
                    ClusterTree_T, T1, T2,
                    BlockClusterTree_T::IsSymmetric(),
                    energy_flag, diff_flag, hess_flag, metric_flag, prec_flag
                >;
                
                Kernel_T ker ( conf, bct.NearFieldSeparationParameter(), 20, q_half_, p_half_ );
                
                ptoc(ClassName()+"::Compute_VF prepare kernels");
                
                ptic(ClassName()+"::Compute_VF prepare traversor");
                FMM_Traversor<Kernel_T> traversor ( bct.VeryNear(), ker );
                ptoc(ClassName()+"::Compute_VF prepare traversor");
                
                ptic(ClassName()+"::Compute_VF traversal");
                en += traversor.Compute();
                ptoc(ClassName()+"::Compute_VF traversal");
            
                if constexpr ( is_symmetric && metric_flag )
                {
                    wprint(ClassName()+"::Compute_VF: TODO symmetrize metric!");
                }
            
            ptoc(ClassName()+"::Compute_VF");
        }
            
        template< typename T1, typename T2 >
        void Compute_NF( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::Compute_NF");
            
                ptic(ClassName()+"::Compute_NF prepare kernels");
                
                Configurator_T conf (
                    bct.GetS(), bct.GetT(),
                    metric_values[FMM::Type::NF], prec_values[FMM::Type::NF]
                );
        
                using Kernel_T = TP_Kernel_NF<
                    S_DOM_DIM, T_DOM_DIM,
                    ClusterTree_T, T1, T2,
                    BlockClusterTree_T::IsSymmetric(),
                    energy_flag, diff_flag, hess_flag, metric_flag, prec_flag
                >;
                
                Kernel_T ker ( conf, q_half_, p_half_ );
                
                ptoc(ClassName()+"::Compute_NF prepare kernels");
                
                ptic(ClassName()+"::Compute_NF prepare traversor");
                FMM_Traversor<Kernel_T> traversor ( bct.Near(), ker );
                ptoc(ClassName()+"::Compute_NF prepare traversor");
                
                ptic(ClassName()+"::Compute_NF traversal");
                en += traversor.Compute();
                ptoc(ClassName()+"::Compute_NF traversal");
            
                if constexpr ( is_symmetric && metric_flag  )
                {
                    wprint(ClassName()+"::Compute_NF: TODO symmetrize metric!");
                }
            
            ptoc(ClassName()+"::Compute_NF");
        }
            
        template< typename T1, typename T2 >
        void Compute_FF( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::Compute_FF");
            
            ptic(ClassName()+"::Compute_FF prepare kernels");
            
            Configurator_T conf (
                                 bct.GetS(), bct.GetT(),
                                 metric_values[FMM::Type::FF], prec_values[FMM::Type::FF]
                                 );
            
            using Kernel_T = TP_Kernel_FF<
            ClusterTree_T, T1, T2,
            BlockClusterTree_T::IsSymmetric(),
            energy_flag, diff_flag, hess_flag, metric_flag, prec_flag
            >;
            
            Kernel_T ker ( conf, q_half_, p_half_ );
            
            ptoc(ClassName()+"::Compute_FF prepare kernels");
            
            ptic(ClassName()+"::Compute_FF prepare traversor");
            FMM_Traversor<Kernel_T> traversor ( bct.Far(), ker );
            ptoc(ClassName()+"::Compute_FF prepare traversor");
            
            ptic(ClassName()+"::Compute_FF traversal");
            en += traversor.Compute();
            ptoc(ClassName()+"::Compute_FF traversal");
            
            if( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    wprint(ClassName()+"::Compute_FF: TODO symmetrize metric!");
                }
                
                wprint(ClassName()+"::Compute_FF: TODO symmetrize metric!");
            }
            ptoc(ClassName()+"::Compute_FF");
        }
        

        void MultiplyMetric() const
        {
            const Int rhs_count = bct.GetS().BufferDimension()/(1+AMB_DIM);
            if( metric_flag )
            {
                MultiplyMetric_FF();
                MultiplyMetric_NF();
                MultiplyMetric_VF();
            }
        }
        

        void MultiplyMetric_FF( const Int rhs_count ) const
        {
            using Kernel_T = TP_Kernel_FF_MultiplyMetric<
                AMB_DIM,AMB_DIM,
                Real,Int,ExtReal,ExtReal,
                true,true,1,0
            >;
            
            if( metric_flag )
            {
                SparseKernelMatrixCSR_gen<Kernel_T> matrix ( bct.Near() );
                
                matrix.Dot(
                    metric_values[FMM::Type::FF].data()
                    bct.GetT().ClusterInputBuffer().data()
                    bct.GetS().ClusterOutputBuffer().data()
                    rhs_count
               );
            }
        }
        
        void MultiplyMetric_NF( const Int rhs_count ) const
        {
            using Kernel_T = TP_Kernel_NF_MultiplyMetric<
                AMB_DIM,AMB_DIM,
                Real,Int,ExtReal,ExtReal,
                true,true,1,0
            >;
            if( metric_flag )
            {
                
                SparseKernelMatrixCSR_gen<Kernel_T> matrix ( bct.Near() );
                
                matrix.Dot(
                    metric_values[FMM::Type::NF].data()
                    bct.GetT().PrimitiveInputBuffer().data()
                    bct.GetS().PrimitiveOutputBuffer().data()
                    rhs_count
               );
            }
        }
        
        void MultiplyMetric_VF( const Int rhs_count ) const
        {
            using Kernel_T = TP_Kernel_VF_MultiplyMetric<
                AMB_DIM,AMB_DIM,
                Real,Int,ExtReal,ExtReal,
                true,true,1,1   // Caution! We have to add-into and not to overwrite!
            >;
            
            if constexpr ( metric_flag )
            {
                SparseKernelMatrixCSR_gen<Kernel_T> matrix ( bct.VeryNear() );
                
                matrix.Dot(
                    metric_values[FMM::Type::VF].data()
                    bct.GetT().PrimitiveInputBuffer().data()
                    bct.GetS().PrimitiveOutputBuffer().data()
                    rhs_count
                );
            }
        }

        
    public:
        
        std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(S_DOM_DIM)+","+ToString(T_DOM_DIM)+","+bct.ClassName()+">";
        }
    };
    
}// namespace Repulsor

#undef CLASS

