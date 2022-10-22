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
        
        using Kernel_VF_MultiplyMetric_T = TP_Kernel_VF_MultiplyMetric<
            AMB_DIM,AMB_DIM,
            Real,Int,Real,Real,
            true,true,1,1 // Caution! We have to add-into and not to overwrite!
        >;
        
        using Kernel_NF_MultiplyMetric_T = TP_Kernel_NF_MultiplyMetric<
            AMB_DIM,AMB_DIM,
            Real,Int,Real,Real,
            true,true,1,0
        >;
        
        using Kernel_FF_MultiplyMetric_T = TP_Kernel_FF_MultiplyMetric<
            AMB_DIM,AMB_DIM,
            Real,Int,Real,Real,
            true,true,1,0
        >;
        
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
            
            if( diff_flag )
            {
                bct.GetS().CleanseDerivativeBuffers();
                
                if( !is_symmetric )
                {
                    bct.GetT().CleanseDerivativeBuffers();
                }
            }
            
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
                
                FMM_Traversor<Kernel_T> traversor ( bct.VeryNear(), ker );
                en += traversor.Compute();

            
            if( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    SparseKernelMatrixCSR_gen<Kernel_VF_MultiplyMetric_T> matrix ( bct.VeryNear() );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( metric_values[FMM::Type::VF].data() );
                }
                
                wprint(ClassName()+"::Compute_VF: TODO: Correct the diagonal blocks!");
            }
            
            ptoc(ClassName()+"::Compute_VF");
        }
            
        template< typename T1, typename T2 >
        void Compute_NF( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::Compute_NF");
            
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

                FMM_Traversor<Kernel_T> traversor ( bct.Near(), ker );

                en += traversor.Compute();
            
            if( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    SparseKernelMatrixCSR_gen<Kernel_NF_MultiplyMetric_T> matrix ( bct.Near() );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( metric_values[FMM::Type::NF].data() );
                }
                
                wprint(ClassName()+"::Compute_NF: TODO: Correct the diagonal blocks!");
            }
            
            ptoc(ClassName()+"::Compute_NF");
        }
            
        template< typename T1, typename T2 >
        void Compute_FF( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::Compute_FF");
            
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
            
            FMM_Traversor<Kernel_T> traversor ( bct.Far(), ker );
            
            en += traversor.Compute();
            
            if( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    SparseKernelMatrixCSR_gen<Kernel_FF_MultiplyMetric_T> matrix ( bct.Far() );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( metric_values[FMM::Type::FF].data() );
                }
                
                wprint(ClassName()+"::Compute_FF: TODO: Correct the diagonal blocks!");
            }
            ptoc(ClassName()+"::Compute_FF");
        }
        
    public:
        
        void MultiplyMetric() const
        {
            const Int rhs_count = bct.GetS().BufferDimension()/(1+AMB_DIM);
            
            if( metric_flag )
            {
                MultiplyMetric_NF(rhs_count);
                MultiplyMetric_VF(rhs_count);
                MultiplyMetric_FF(rhs_count);
            }
        }
        
    protected:

        void MultiplyMetric_NF( const Int rhs_count ) const
        {
            if( metric_flag )
            {
                SparseKernelMatrixCSR_gen<Kernel_NF_MultiplyMetric_T> matrix ( bct.Near() );
                
                matrix.Dot(
                    metric_values[FMM::Type::NF].data(),
                    static_cast<Real>(1), bct.GetT().PrimitiveInputBuffer().data(),
                    static_cast<Real>(0), bct.GetS().PrimitiveOutputBuffer().data(),
                    rhs_count
               );
            }
        }
        
        void MultiplyMetric_VF( const Int rhs_count ) const
        {
            if constexpr ( metric_flag )
            {
                SparseKernelMatrixCSR_gen<Kernel_VF_MultiplyMetric_T> matrix ( bct.VeryNear() );
                
                matrix.Dot(
                    metric_values[FMM::Type::VF].data(),
                    static_cast<Real>(1), bct.GetT().PrimitiveInputBuffer().data(),
                    static_cast<Real>(1), bct.GetS().PrimitiveOutputBuffer().data(),
                    rhs_count
                );
            }
        }
        
        void MultiplyMetric_FF( const Int rhs_count ) const
        {
            if( metric_flag )
            {
                SparseKernelMatrixCSR_gen<Kernel_FF_MultiplyMetric_T> matrix ( bct.Near() );
                
                matrix.Dot(
                    metric_values[FMM::Type::FF].data(),
                    static_cast<Real>(1), bct.GetT().ClusterInputBuffer().data(),
                    static_cast<Real>(0), bct.GetS().ClusterOutputBuffer().data(),
                    rhs_count
               );
            }
        }
        
        
    public:
        
        std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"
                + ToString(S_DOM_DIM)   + ","
                + ToString(T_DOM_DIM)   + ","
                + bct.ClassName()       + ","
                + ToString(energy_flag) + ","
                + ToString(diff_flag)   + ","
                + ToString(hess_flag)   + ","
                + ToString(metric_flag) + ","
                + ToString(prec_flag)
                + ">";
        }
    };
    
}// namespace Repulsor

#undef CLASS

