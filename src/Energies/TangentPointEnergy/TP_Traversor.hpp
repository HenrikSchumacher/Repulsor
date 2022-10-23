#pragma once

#define CLASS TP_Traversor

namespace Repulsor
{
    using namespace FMM;
    
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename BlockClusterTree_T,
        bool energy_flag, bool diff_flag, bool metric_flag
    >
    class CLASS
    {
    public:
        
        using ClusterTree_T     = typename BlockClusterTree_T::ClusterTree_T;
        using Real              = typename BlockClusterTree_T::Real;
        using Int               = typename BlockClusterTree_T::Int;
        using SReal             = typename BlockClusterTree_T::SReal;
        using ExtReal           = typename BlockClusterTree_T::ExtReal;
        
        static constexpr Int  AMB_DIM      = BlockClusterTree_T::AMB_DIM;
        static constexpr Int  S_DOM_DIM    = S_DOM_DIM_;
        static constexpr Int  T_DOM_DIM    = T_DOM_DIM_;
        static constexpr bool is_symmetric = BlockClusterTree_T::IsSymmetric();
        
        using TangentVector_T   = Tensor2<ExtReal,Int>;
        using CotangentVector_T = Tensor2<ExtReal,Int>;
        
        using Values_T          = Tensor2<Real,Int>;
        using ValueContainer_T  = std::array<Values_T,3>;
        
        using Configurator_T    = FMM_Configurator<ClusterTree_T>;

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
        
    protected:
        
        static constexpr Real zero = 0;
        static constexpr Real one  = 1;
        static constexpr Real two  = 2;
        
    
    public:
        
        CLASS() = delete;
        
        CLASS(
              const BlockClusterTree_T & bct_,
              ValueContainer_T & metric_values,
              const Real q_,
              const Real p_
        )
        :   bct(bct_)
        ,   metric_values ( metric_values )
        ,   q             ( q_ )
        ,   p             ( p_ )
        ,   q_half_real   ( q / two )
        ,   p_half_real   ( p / two )
        ,   q_half_int    ( q_half_real )
        ,   p_half_int    ( p_half_real )
        ,   q_half_is_int ( q_half_real == static_cast<Real>(q_half_int) )
        ,   p_half_is_int ( p_half_real == static_cast<Real>(p_half_int) )
        {}

        ~CLASS() = default;
        

    protected:

        const BlockClusterTree_T & bct;
        
        ValueContainer_T & metric_values;
        
        const Real q;
        const Real p;
        const Real q_half_real;
        const Real p_half_real;
        
        const Int q_half_int;
        const Int p_half_int;
        
        const bool q_half_is_int;
        const bool p_half_is_int;
        
        Real en = zero;
        
    public:
        

//##############################################################################################
//      Compute
//##############################################################################################
        
        Real Compute()
        {
            en = zero;
            
            if constexpr ( diff_flag )
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
                    VF_Compute<Int,Int>( q_half_int, p_half_int );
                    NF_Compute<Int,Int>( q_half_int, p_half_int );
                    FF_Compute<Int,Int>( q_half_int, p_half_int );
                }
                else
                {
                    VF_Compute<Int,Real>( q_half_int, p_half_real );
                    NF_Compute<Int,Real>( q_half_int, p_half_real );
                    FF_Compute<Int,Real>( q_half_int, p_half_real );
                }
            }
            else
            {
                if( p_half_is_int)
                {
                    VF_Compute<Real,Int>( q_half_real, p_half_int );
                    NF_Compute<Real,Int>( q_half_real, p_half_int );
                    FF_Compute<Real,Int>( q_half_real, p_half_int );
                }
                else
                {
                    VF_Compute<Real,Real>( q_half_real, p_half_real );
                    NF_Compute<Real,Real>( q_half_real, p_half_real );
                    FF_Compute<Real,Real>( q_half_real, p_half_real );
                }
            }
            
            return en;
        }
        
//##############################################################################################
//      Compute subroutines
//##############################################################################################
    
    protected:
        
        template< typename T1, typename T2 >
        void VF_Compute( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::VF_Compute");
            
            const Type type = Type::VF;
            const int idx = std::__to_underlying(type);
            
            auto & pattern = bct.template GetMatrix<type>();
        
            Configurator_T conf ( bct.GetS(), bct.GetT(), metric_values[idx] );

            using Kernel_T = TP_Kernel_VF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, T1, T2,
                BlockClusterTree_T::IsSymmetric(),
                energy_flag, diff_flag, metric_flag
            >;

            Kernel_T ker ( conf, bct.NearFieldSeparationParameter(), 20, q_half_, p_half_ );

            FMM_Traversor<Kernel_T> traversor ( pattern, ker );
            en += traversor.Compute();

            if constexpr ( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    SparseKernelMatrixCSR_gen<Kernel_NF_MultiplyMetric_T> matrix ( pattern );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( metric_values[idx].data() );
                }
            }
            
            ptoc(ClassName()+"::VF_Compute");
        }
            
        template< typename T1, typename T2 >
        void NF_Compute( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::NF_Compute");
            
            const Type type = Type::NF;
            const int idx = std::__to_underlying(type);
            
            auto & pattern = bct.template GetMatrix<type>();
            
            Configurator_T conf ( bct.GetS(), bct.GetT(), metric_values[idx] );
    
            using Kernel_T = TP_Kernel_NF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, T1, T2,
                BlockClusterTree_T::IsSymmetric(),
                energy_flag, diff_flag, metric_flag
            >;
            
            Kernel_T ker ( conf, q_half_, p_half_ );

            FMM_Traversor<Kernel_T> traversor ( pattern, ker );

            en += traversor.Compute();

            if constexpr ( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    SparseKernelMatrixCSR_gen<Kernel_NF_MultiplyMetric_T> matrix ( pattern );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( metric_values[idx].data() );
                }
            }
            
            ptoc(ClassName()+"::NF_Compute");
        }
            
        template< typename T1, typename T2 >
        void FF_Compute( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::FF_Compute");
            
            const Type type = Type::FF;
            const int idx = std::__to_underlying(type);
            
            auto & pattern = bct.template GetMatrix<type>();
            
            Configurator_T conf ( bct.GetS(), bct.GetT(), metric_values[idx] );
            
            using Kernel_T = TP_Kernel_FF<
                ClusterTree_T, T1, T2,
                BlockClusterTree_T::IsSymmetric(),
                energy_flag, diff_flag, metric_flag
            >;
            
            Kernel_T ker ( conf, q_half_, p_half_ );
            
            FMM_Traversor<Kernel_T> traversor ( pattern, ker );
            
            en += traversor.Compute();
            
            if constexpr ( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    SparseKernelMatrixCSR_gen<Kernel_NF_MultiplyMetric_T> matrix ( pattern );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( metric_values[idx].data() );
                }
            }
            
            ptoc(ClassName()+"::FF_Compute");
        }
        
//##############################################################################################
//      Multiply metric
//##############################################################################################
        
    public:
        
        void MultiplyMetric(const bool VF_flag, const bool NF_flag, const bool FF_flag ) const
        {
            if constexpr ( metric_flag )
            {
                const Int rhs_count = bct.GetS().BufferDimension() / Kernel_NF_MultiplyMetric_T::ROWS;                
         
                if( NF_flag )
                {
                    NF_MultiplyMetric(rhs_count);
                }
                else
                {
                    bct.GetS().PrimitiveOutputBuffer().SetZero();
                }
                
                if( VF_flag )
                {
                    VF_MultiplyMetric(rhs_count);
                }
                
                if( FF_flag )
                {
                    FF_MultiplyMetric(rhs_count);
                }
                else
                {
                    bct.GetS().ClusterOutputBuffer().SetZero();
                }
            }
        }
        
//##############################################################################################
//      Multiply metric subroutines
//##############################################################################################
        
    protected:

        void NF_MultiplyMetric( const Int rhs_count ) const
        {
            constexpr Type type = Type::NF;
            
            const int idx = std::__to_underlying(type);
            auto & pattern = bct.template GetMatrix<type>();
            
            SparseKernelMatrixCSR_gen<Kernel_NF_MultiplyMetric_T> matrix ( pattern );
            
            matrix.Dot(
                metric_values[idx].data(),
                one,  bct.GetT().template GetBuffer<type,InOut::In >().data(),
                zero, bct.GetS().template GetBuffer<type,InOut::Out>().data(),
                rhs_count
           );
        }
        
        void VF_MultiplyMetric( const Int rhs_count ) const
        {
            if constexpr ( metric_flag )
            {
                constexpr Type type = Type::VF;
                
                const int idx = std::__to_underlying(type);
                auto & pattern = bct.template GetMatrix<type>();
                
                SparseKernelMatrixCSR_gen<Kernel_VF_MultiplyMetric_T> matrix ( pattern );
                
                matrix.Dot(
                    metric_values[idx].data(),
                    one,  bct.GetT().template GetBuffer<type,InOut::In >().data(),
                    one,  bct.GetS().template GetBuffer<type,InOut::Out>().data(),
                    rhs_count
                );
            }
        }
        
        void FF_MultiplyMetric( const Int rhs_count ) const
        {
            if constexpr ( metric_flag )
            {
                constexpr Type type = Type::FF;
                
                const int idx = std::__to_underlying(type);
                auto & pattern = bct.template GetMatrix<type>();
                
                SparseKernelMatrixCSR_gen<Kernel_FF_MultiplyMetric_T> matrix ( pattern );
                
                matrix.Dot(
                    metric_values[idx].data(),
                    one,  bct.GetT().template GetBuffer<type,InOut::In >().data(),
                    zero, bct.GetS().template GetBuffer<type,InOut::Out>().data(),
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
                + ToString(metric_flag) +
                + ">";
        }
    };
    
}// namespace Repulsor

#undef CLASS

