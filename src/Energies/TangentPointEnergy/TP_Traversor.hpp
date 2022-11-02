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
        using LInt              = typename BlockClusterTree_T::LInt;
        using SReal             = typename BlockClusterTree_T::SReal;
        using ExtReal           = typename BlockClusterTree_T::ExtReal;
        using Pattern_T         = typename BlockClusterTree_T::Pattern_T;
        using Values_T          = typename BlockClusterTree_T::Values_T;
        using ValueContainer_T  = typename BlockClusterTree_T::ValueContainer_T;
        
        using Accumulator_T     = typename ClusterTree_T::Accumulator_T;
        
        static constexpr Int  AMB_DIM      = BlockClusterTree_T::AMB_DIM;
        static constexpr Int  S_DOM_DIM    = S_DOM_DIM_;
        static constexpr Int  T_DOM_DIM    = T_DOM_DIM_;
        static constexpr bool is_symmetric = BlockClusterTree_T::IsSymmetric();
        
        using TangentVector_T   = Tensor2<ExtReal,Int>;
        using CotangentVector_T = Tensor2<ExtReal,Int>;
        
        using Configurator_T    = FMM_Configurator<BlockClusterTree_T>;

        using Kernel_Block_MulAdd_T = ArrowheadBlockKernel_fixed<
            AMB_DIM+1, AMB_DIM+1, AMB_DIM, true,
            Real, Real, Real, Int, LInt,
            1,    1,                    // CAUTION: We use add-in instead of overwrite!
                         true,
            true, false, true, true,
            true, false,
            true
        >;
        
        using Kernel_Block_Mul_T = ArrowheadBlockKernel_fixed<
            AMB_DIM+1, AMB_DIM+1, AMB_DIM, true,
            Real, Real, Real, Int, LInt,
            1,    0,
                         true,
            true, false, true, true,
            true, false,
            true
        >;
        
        using Kernel_Diag_MulAdd_T = DenseBlockKernel_fixed<
            AMB_DIM+1, AMB_DIM+1, AMB_DIM, true,
            Real, Real, Real, Int, LInt,
            1,    1,
            true, false, true,
            true, false, true, true,
            true, false,
            1   , 2,
            true
        >;
    
        static constexpr Int VF_blk_size = Kernel_Block_MulAdd_T::ROWS * Kernel_Block_MulAdd_T::COLS;
        static constexpr Int NF_blk_size =    Kernel_Block_Mul_T::ROWS *    Kernel_Block_Mul_T::COLS;
        static constexpr Int FF_blk_size =    Kernel_Block_Mul_T::ROWS *    Kernel_Block_Mul_T::COLS;
        
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
        :   bct           ( bct_          )
        ,   S             ( bct.GetS()    )
        ,   T             ( bct.GetT()    )
        ,   metric_values ( metric_values )
        ,   q             ( q_            )
        ,   p             ( p_            )
        ,   q_half_real   ( q / two       )
        ,   p_half_real   ( p / two       )
        ,   q_half_int    ( q_half_real   )
        ,   p_half_int    ( p_half_real   )
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
                bct.GetS().CleanseDerivativeBuffers();
                
                if( !is_symmetric )
                {
                    T.CleanseDerivativeBuffers();
                }
            }
            if constexpr ( metric_flag )
            {
                ptic("Allocate accumulators");
                
                S.VF_Accumulator()
                    = Accumulator_T( thread_count, S.PrimitiveCount(), VF_blk_size, 0 );
                S.NF_Accumulator()
                    = Accumulator_T( thread_count, S.PrimitiveCount(), NF_blk_size, 0 );
                S.FF_Accumulator()
                    = Accumulator_T( thread_count, S.ClusterCount(),   FF_blk_size, 0 );

                if constexpr ( !is_symmetric )
                {
                    T.VF_Accumulator()
                    = Accumulator_T( thread_count, T.PrimitiveCount(), VF_blk_size, 0 );
                    T.NF_Accumulator()
                    = Accumulator_T( thread_count, T.PrimitiveCount(), NF_blk_size, 0 );
                    T.FF_Accumulator()
                    = Accumulator_T( thread_count, T.ClusterCount(),   FF_blk_size, 0 );
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

            S.VF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
            S.NF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
            S.FF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );

            if constexpr ( !is_symmetric )
            {
                T.VF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
                T.NF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
                T.FF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
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
            
            const Pattern_T & pattern = bct.VeryNear();
                  Values_T  & values  = metric_values["VF"];
            
            Configurator_T conf ( bct, values );
            
            using Kernel_T = TP_Kernel_VF<
                S_DOM_DIM, T_DOM_DIM,
                BlockClusterTree_T, T1, T2,
                energy_flag, diff_flag, metric_flag
            >;

            Kernel_T ker ( conf, bct.NearFieldSeparationParameter(), 20, q_half_, p_half_ );

            FMM_Traversor<Pattern_T,Kernel_T> traversor ( pattern, ker );
            en += traversor.Compute();

            if constexpr ( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    SparseKernelMatrixCSR<Kernel_Block_Mul_T> matrix ( pattern );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( values.data() );
                }
                
                ptic("Reduce VF_Accumulators");
                metric_values["VF_diag"] = S.VF_Accumulator().template AddReduce<Real,LInt>();
                ptoc("Reduce VF_Accumulators");
            }
            
            ptoc(ClassName()+"::VF_Compute");
        }
            
        
        template< typename T1, typename T2 >
        void NF_Compute( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::NF_Compute");

            const Pattern_T & pattern = bct.Near();
                  Values_T  & values  = metric_values["NF"];
            
            Configurator_T conf ( bct, values );

            using Kernel_T = TP_Kernel_NF<
                S_DOM_DIM, T_DOM_DIM,
                BlockClusterTree_T, T1, T2,
                energy_flag, diff_flag, metric_flag
            >;

            Kernel_T ker ( conf, q_half_, p_half_ );

            FMM_Traversor<Pattern_T,Kernel_T> traversor ( pattern, ker );

            en += traversor.Compute();

            if constexpr ( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    SparseKernelMatrixCSR<Kernel_Block_Mul_T> matrix ( pattern );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( values.data() );
                }

                ptic("Reduce NF_Accumulators");
                metric_values["NF_diag"] = S.NF_Accumulator().template AddReduce<Real,LInt>();
                ptoc("Reduce NF_Accumulators");
            }

            ptoc(ClassName()+"::NF_Compute");
        }
        
            
        template< typename T1, typename T2 >
        void FF_Compute( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::FF_Compute");
            
            const Pattern_T & pattern = bct.Far();
                  Values_T  & values  = metric_values["FF"];
            
            Configurator_T conf ( bct, values );
            
            using Kernel_T = TP_Kernel_FF<
                BlockClusterTree_T, T1, T2,
                energy_flag, diff_flag, metric_flag
            >;
            
            Kernel_T ker ( conf, q_half_, p_half_ );
            
            FMM_Traversor<Pattern_T,Kernel_T> traversor ( pattern, ker );
            
            en += traversor.Compute();
            
            if constexpr ( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    SparseKernelMatrixCSR<Kernel_Block_Mul_T> matrix ( pattern );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( values.data() );
                }
                
                ptic("Reduce FF_Accumulators");
                metric_values["FF_diag"] = S.FF_Accumulator().template AddReduce<Real,LInt>();
                ptoc("Reduce FF_Accumulators");
            }
            
            ptoc(ClassName()+"::FF_Compute");
        }
        
//##############################################################################################
//      Multiply metric
//##############################################################################################
        
    public:
        
        void MultiplyMetric(
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const
        {
            if constexpr ( metric_flag )
            {
                const Int rhs_count = S.BufferDimension() / Kernel_Block_Mul_T::ROWS;
                
                if( NF_flag )
                {
                    NF_MultiplyMetric(rhs_count);
                }
                else
                {
                    S.PrimitiveOutputBuffer().SetZero();
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
            SparseKernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.Near() );
            
            matrix.Dot(
                metric_values["NF"].data(),
                one,  T.PrimitiveInputBuffer().data(),
                zero, S.PrimitiveOutputBuffer().data(),
                rhs_count
            );
            
            if constexpr ( is_symmetric )
            {
                DiagonalKernelMatrix<Kernel_Diag_MulAdd_T> diag (
                    bct.Near().RowCount(),
                    bct.Near().ThreadCount()
                );

                diag.Dot(
                    metric_values["NF_diag"].data(),
                    one, T.PrimitiveInputBuffer().data(),
                    one, S.PrimitiveOutputBuffer().data(),
                    rhs_count
                );
            }
        }
        
        void VF_MultiplyMetric( const Int rhs_count ) const
        {
            SparseKernelMatrixCSR<Kernel_Block_MulAdd_T> matrix ( bct.VeryNear() );
            
            matrix.Dot(
                metric_values["VF"].data(),
                one, T.PrimitiveInputBuffer().data(),
                one, S.PrimitiveOutputBuffer().data(),
                rhs_count
            );
            
            if constexpr ( is_symmetric )
            {
                DiagonalKernelMatrix<Kernel_Diag_MulAdd_T> diag (
                    bct.VeryNear().RowCount(),
                    bct.VeryNear().ThreadCount()
                );

                diag.Dot(
                    metric_values["VF_diag"].data(),
                    one, T.PrimitiveInputBuffer().data(),
                    one, S.PrimitiveOutputBuffer().data(),
                    rhs_count
                );
            }
        }
        
        void FF_MultiplyMetric( const Int rhs_count ) const
        {
            SparseKernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.Far() );
            
            matrix.Dot(
                metric_values["FF"].data(),
                one,  T.ClusterInputBuffer().data(),
                zero, S.ClusterOutputBuffer().data(),
                rhs_count
           );
            
            if constexpr ( is_symmetric )
            {
                DiagonalKernelMatrix<Kernel_Diag_MulAdd_T> diag (
                    bct.Far().RowCount(),
                    bct.Far().ThreadCount()
                );
                
                diag.Dot(
                    metric_values["FF_diag"].data(),
                    one,  T.ClusterInputBuffer().data(),
                    one, S.ClusterOutputBuffer().data(),
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
