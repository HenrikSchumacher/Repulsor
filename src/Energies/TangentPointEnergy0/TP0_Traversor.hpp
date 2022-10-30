#pragma once

#define CLASS TP0_Traversor

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
        using LInt              = typename BlockClusterTree_T::LInt;
        using SReal             = typename BlockClusterTree_T::SReal;
        using ExtReal           = typename BlockClusterTree_T::ExtReal;
        using Accumulator_T     = typename ClusterTree_T::Accumulator_T;
        
        static constexpr Int  AMB_DIM      = BlockClusterTree_T::AMB_DIM;
        static constexpr Int  S_DOM_DIM    = S_DOM_DIM_;
        static constexpr Int  T_DOM_DIM    = T_DOM_DIM_;
        static constexpr bool is_symmetric = BlockClusterTree_T::IsSymmetric();
        

        
        using TangentVector_T   = Tensor2<ExtReal,Int>;
        using CotangentVector_T = Tensor2<ExtReal,Int>;
        
        using Values_T          = Tensor2<Real,LInt>;
        using ValueContainer_T  = std::unordered_map<std::string,Values_T>;
        
        using Configurator_T    = FMM_Configurator<ClusterTree_T>;

        using Kernel_VF_Multiply_T = TP0_Kernel_VF_MultiplyMetric<
            AMB_DIM,AMB_DIM,
            Real,Real,Real,Int,LInt,
            true,true,1,1 // Caution! We have to add-into and not to overwrite!
        >;
        
        using Kernel_NF_Multiply_T = TP0_Kernel_NF_MultiplyMetric<
            AMB_DIM,AMB_DIM,
            Real,Real,Real,Int,LInt,
            true,true,1,0
        >;
        
        using Kernel_FF_Multiply_T = TP0_Kernel_FF_MultiplyMetric<
            AMB_DIM,AMB_DIM,
            Real,Real,Real,Int,LInt,
            true,true,1,0
        >;
    
        static constexpr Int VF_block_size = Kernel_VF_Multiply_T::ROWS * Kernel_VF_Multiply_T::COLS;
        static constexpr Int NF_block_size = Kernel_NF_Multiply_T::ROWS * Kernel_NF_Multiply_T::COLS;
        static constexpr Int FF_block_size = Kernel_FF_Multiply_T::ROWS * Kernel_FF_Multiply_T::COLS;
        
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
        :   bct( bct_ )
        ,   S( bct.GetS() )
        ,   T( bct.GetT() )
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
        const ClusterTree_T      & S;
        const ClusterTree_T      & T;
        
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
            const Int thread_count = bct.ThreadCount();
            
            en = zero;
            
            DummyAllocators();
            
            if constexpr ( diff_flag )
            {
                S.CleanseDerivativeBuffers();
                
                if( !is_symmetric )
                {
                    T.CleanseDerivativeBuffers();
                }
            }
            if constexpr ( metric_flag )
            {
                ptic("Allocate accumulators");

                metric_values["VF_diag"] = Values_T (S.PrimitiveCount(), VF_block_size );
                metric_values["NF_diag"] = Values_T (S.PrimitiveCount(), NF_block_size );
                metric_values["FF_diag"] = Values_T (S.ClusterCount(),   FF_block_size );

                S.VF_Accumulator() = Accumulator_T( thread_count, S.PrimitiveCount(), VF_block_size );
                S.NF_Accumulator() = Accumulator_T( thread_count, S.PrimitiveCount(), NF_block_size );
                S.FF_Accumulator() = Accumulator_T( thread_count, S.ClusterCount(),   FF_block_size );

                if constexpr ( !is_symmetric )
                {
                    T.VF_Accumulator() = Accumulator_T( thread_count, T.PrimitiveCount(), VF_block_size );
                    T.NF_Accumulator() = Accumulator_T( thread_count, T.PrimitiveCount(), NF_block_size );
                    T.FF_Accumulator() = Accumulator_T( thread_count, T.ClusterCount(),   FF_block_size );
                }
                ptoc("Allocate accumulators");
            }
            else
            {
                DummyAllocators();
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
            
            if constexpr ( metric_flag )
            {
                DummyAllocators();
            }
            
            return en;
        }
     
        void DummyAllocators()
        {
            ptic("Allocate dummy accumulators");

            const Int thread_count = bct.ThreadCount();

            S.VF_Accumulator() = Accumulator_T( thread_count, 1, 1 );
            S.NF_Accumulator() = Accumulator_T( thread_count, 1, 1 );
            S.FF_Accumulator() = Accumulator_T( thread_count, 1, 1 );

            if constexpr ( !is_symmetric )
            {
                T.VF_Accumulator() = Accumulator_T( thread_count, 1, 1 );
                T.NF_Accumulator() = Accumulator_T( thread_count, 1, 1 );
                T.FF_Accumulator() = Accumulator_T( thread_count, 1, 1 );
            }
            ptoc("Allocate dummy accumulators");
        }
        

        
//##############################################################################################
//      Compute subroutines
//##############################################################################################
    
    protected:
        
        template< typename T1, typename T2 >
        void VF_Compute( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::VF_Compute");
            
            auto & pattern = bct.VeryNear();
            
            Configurator_T conf ( S, T, metric_values["VF"] );
            
            using Kernel_T = TP0_Kernel_VF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, T1, T2,
                BlockClusterTree_T::IsSymmetric(),
                energy_flag, diff_flag, metric_flag
            >;

            Kernel_T ker ( conf, bct.NearFieldSeparationParameter(), 20, q_half_, p_half_ );

            FMM_Traversor<Kernel_T,LInt> traversor ( pattern, ker );
            en += traversor.Compute();

            if constexpr ( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    SparseKernelMatrixCSR<Kernel_NF_Multiply_T> matrix ( pattern );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( metric_values["VF"].data() );
                }
                
//                ptic("Reduce VF_Accumulators");
//                // We overwrite diag_values here so that we do not have to zerofy it first.
//                S.VF_Accumulator().AddReduce( metric_values["VF_diag"], false );
//                ptoc("Reduce VF_Accumulators");
            }
            
            ptoc(ClassName()+"::VF_Compute");
        }
            
        
        template< typename T1, typename T2 >
        void NF_Compute( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::NF_Compute");

            auto & pattern = bct.Near();

            Configurator_T conf ( S, T, metric_values["NF"] );

            using Kernel_T = TP0_Kernel_NF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, T1, T2,
                BlockClusterTree_T::IsSymmetric(),
                energy_flag, diff_flag, metric_flag
            >;

            Kernel_T ker ( conf, q_half_, p_half_ );

            FMM_Traversor<Kernel_T,LInt> traversor ( pattern, ker );

            en += traversor.Compute();

            if constexpr ( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    SparseKernelMatrixCSR<Kernel_NF_Multiply_T> matrix ( pattern );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( metric_values["NF"].data() );
                }

                ptic("Reduce NF_Accumulators");
                // We overwrite diag_values here so that we do not have to zerofy it first.
                S.NF_Accumulator().AddReduce( metric_values["NF_diag"], false );
                ptoc("Reduce NF_Accumulators");
            }

            ptoc(ClassName()+"::NF_Compute");
        }
        
            
        template< typename T1, typename T2 >
        void FF_Compute( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::FF_Compute");
            
            auto & pattern = bct.Far();
            
            Configurator_T conf ( S, T, metric_values["FF"] );
            
            using Kernel_T = TP0_Kernel_FF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, T1, T2,
                BlockClusterTree_T::IsSymmetric(),
                energy_flag, diff_flag, metric_flag
            >;
            
            Kernel_T ker ( conf, q_half_, p_half_ );
            
            FMM_Traversor<Kernel_T,LInt> traversor ( pattern, ker );
            
            en += traversor.Compute();
            
            if constexpr ( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    SparseKernelMatrixCSR<Kernel_NF_Multiply_T> matrix ( pattern );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( metric_values["FF"].data() );
                }
                
                ptic("Reduce FF_Accumulators");
                // We overwrite diag_values here so that we do not have to zerofy it first.
                S.FF_Accumulator().AddReduce( metric_values["FF_diag"], false );
                ptoc("Reduce FF_Accumulators");
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
                const Int rhs_count = S.BufferDimension() / Kernel_NF_Multiply_T::ROWS;
         
                if( NF_flag )
                {
                    NF_MultiplyMetric(rhs_count);
                    // TODO: Apply diagonal
                }
                else
                {
                    S.PrimitiveOutputBuffer().SetZero();
                }
                
                if( VF_flag )
                {
                    VF_MultiplyMetric(rhs_count);
                    // TODO: Apply diagonal
                }
                
                
                if( FF_flag )
                {
                    FF_MultiplyMetric(rhs_count);
                    // TODO: Apply diagonal
                }
                else
                {
                    S.ClusterOutputBuffer().SetZero();
                }
            }
        }
        
//##############################################################################################
//      Multiply metric subroutines
//##############################################################################################
        
    protected:

        void NF_MultiplyMetric( const Int rhs_count ) const
        {
            SparseKernelMatrixCSR<Kernel_NF_Multiply_T> matrix ( bct.Near() );
            
            matrix.Dot(
                metric_values["NF"].data(), metric_values["NF_diag"].data(),
                one,  T.PrimitiveInputBuffer().data(),
                zero, S.PrimitiveOutputBuffer().data(),
                rhs_count
           );
        }
        
        
        void VF_MultiplyMetric( const Int rhs_count ) const
        {
            if constexpr ( metric_flag )
            {
                SparseKernelMatrixCSR<Kernel_VF_Multiply_T> matrix ( bct.VeryNear() );
                
                matrix.Dot(
                    metric_values["VF"].data(), metric_values["VF_diag"].data(),
                    one,  T.PrimitiveInputBuffer().data(),
                    one,  S.PrimitiveOutputBuffer().data(),
                    rhs_count
                );
            }
        }
        
        void FF_MultiplyMetric( const Int rhs_count ) const
        {
            if constexpr ( metric_flag )
            {
                SparseKernelMatrixCSR<Kernel_FF_Multiply_T> matrix ( bct.Far() );
                
                matrix.Dot(
                    metric_values["FF"].data(), metric_values["FF_diag"].data(),
                    one,  T.ClusterInputBuffer().data(),
                    zero, S.ClusterOutputBuffer().data(),
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

