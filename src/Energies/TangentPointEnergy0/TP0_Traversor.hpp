#pragma once

#include "TP0_Kernel_FF.hpp"
#include "TP0_Kernel_NF.hpp"
#include "TP0_Kernel_VF.hpp"
#include "TP0_Kernel_MultiplyMetric.hpp"

namespace Repulsor
{
    using namespace FMM;
    
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename BlockClusterTree_T_,
        bool energy_flag, bool diff_flag, bool metric_flag
    >
    class TP0_Traversor
    {
    public:
        
        using BlockClusterTree_T = BlockClusterTree_T_;
        using ClusterTree_T      = typename BlockClusterTree_T::ClusterTree_T;
        using Real               = typename ClusterTree_T::Real;
        using Int                = typename ClusterTree_T::Int;
        using LInt               = typename BlockClusterTree_T::LInt;
        using SReal              = typename ClusterTree_T::SReal;
        using ExtReal            = typename ClusterTree_T::ExtReal;
        
        using Configurator_T     = FMM_Configurator<ClusterTree_T>;
        using Values_T           = typename Configurator_T::Values_T;
        using ValueContainer_T   = typename Configurator_T::ValueContainer_T;
        
        using Accumulator_T      = typename ClusterTree_T::Accumulator_T;
        
        static constexpr Int  AMB_DIM      = ClusterTree_T::AMB_DIM;
        static constexpr Int  S_DOM_DIM    = S_DOM_DIM_;
        static constexpr Int  T_DOM_DIM    = T_DOM_DIM_;
        static constexpr bool is_symmetric = BlockClusterTree_T::IsSymmetric();
        
        using TangentVector_T   = Tensor2<ExtReal,Int>;
        using CotangentVector_T = Tensor2<ExtReal,Int>;
        
        

        using Kernel_Block_MulAdd_T = TP0_Kernel_MultiplyMetric<
            AMB_DIM,AMB_DIM,
            Real, Real, Real, Int, LInt,
            1, 1 // <-- this 1 means add-in
        >;
        
        using Kernel_Block_Mul_T = TP0_Kernel_MultiplyMetric<
            AMB_DIM,AMB_DIM,
            Real, Real, Real, Int, LInt,
            1, 0 // <-- this 0 means overwrite
        >;

        
        using Kernel_Diag_MulAdd_T = TP0_Kernel_MultiplyMetric<
            AMB_DIM,AMB_DIM,
            Real, Real, Real, Int, LInt,
            1, 1 // <-- this 1 means add-in
        >;
    
        static constexpr Int VF_blk_size = Kernel_Block_MulAdd_T::ROWS * Kernel_Block_MulAdd_T::COLS;
        static constexpr Int NF_blk_size =    Kernel_Block_Mul_T::ROWS *    Kernel_Block_Mul_T::COLS;
        static constexpr Int FF_blk_size =    Kernel_Block_Mul_T::ROWS *    Kernel_Block_Mul_T::COLS;
        
    protected:
        
        static constexpr Real zero = 0;
        static constexpr Real one  = 1;
        static constexpr Real two  = 2;
        
    
    public:
        
        TP0_Traversor() = delete;
        
        TP0_Traversor(
              const BlockClusterTree_T & bct_,
              ValueContainer_T & metric_values_,
              const Real q_,
              const Real p_
        )
        :   bct           ( bct_          )
        ,   conf          ( bct.GetS(), bct.GetT(), metric_values_ )
        ,   metric_values ( metric_values_ )
        ,   q             ( q_            )
        ,   p             ( p_            )
        ,   q_half_real   ( q / two       )
        ,   p_half_real   ( p / two       )
        ,   q_half_int    ( q_half_real   )
        ,   p_half_int    ( p_half_real   )
        ,   q_half_is_int ( q_half_real == static_cast<Real>(q_half_int) )
        ,   p_half_is_int ( p_half_real == static_cast<Real>(p_half_int) )
        {}

        ~TP0_Traversor() = default;
        

    protected:

        const BlockClusterTree_T & bct;
        
        Configurator_T conf;
        
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
                    bct.GetT().CleanseDerivativeBuffers();
                }
            }

            if constexpr ( metric_flag )
            {
                ptic("Allocate accumulators");
                
                bct.GetS().VF_Accumulator()
                = Accumulator_T( thread_count, bct.GetS().PrimitiveCount(), VF_blk_size, 0 );
                bct.GetS().NF_Accumulator()
                = Accumulator_T( thread_count, bct.GetS().PrimitiveCount(), NF_blk_size, 0 );
                bct.GetS().FF_Accumulator()
                = Accumulator_T( thread_count, bct.GetS().ClusterCount(),   FF_blk_size, 0 );

                if constexpr ( !is_symmetric )
                {
                    bct.GetT().VF_Accumulator()
                    = Accumulator_T( thread_count, bct.GetT().PrimitiveCount(), VF_blk_size, 0 );
                    bct.GetT().NF_Accumulator()
                    = Accumulator_T( thread_count, bct.GetT().PrimitiveCount(), NF_blk_size, 0 );
                    bct.GetT().FF_Accumulator()
                    = Accumulator_T( thread_count, bct.GetT().ClusterCount(),   FF_blk_size, 0 );
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

            bct.GetS().VF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
            bct.GetS().NF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
            bct.GetS().FF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );

            if constexpr ( !is_symmetric )
            {
                bct.GetT().VF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
                bct.GetT().NF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
                bct.GetT().FF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
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
            
            using Kernel_T = TP0_Kernel_VF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, T1, T2, BlockClusterTree_T::is_symmetric,
                energy_flag, diff_flag, metric_flag
            >;

            Kernel_T ker ( conf, bct.NearFieldSeparationParameter(), int_cast<Int>(bct.Settings().max_refinement), q_half_, p_half_ );

            FMM_Traversor traversor ( bct.VeryNear(), ker );
            en += traversor.Compute();

            if constexpr ( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    Sparse::KernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.VeryNear() );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( ker.OffDiag().data() );
                }
                
                ptic("Reduce VF_Accumulators");
                ker.Diag() = bct.GetS().VF_Accumulator().template AddReduce<Real,LInt>();
                ptoc("Reduce VF_Accumulators");
            }
            
            ptoc(ClassName()+"::VF_Compute");
        }
            
        
        template< typename T1, typename T2 >
        void NF_Compute( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::NF_Compute");

            using Kernel_T = TP0_Kernel_NF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, T1, T2,
                BlockClusterTree_T::IsSymmetric(),
                energy_flag, diff_flag, metric_flag
            >;

            Kernel_T ker ( conf, q_half_, p_half_ );

            FMM_Traversor traversor ( bct.Near(), ker );

            en += traversor.Compute();

            if constexpr ( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    Sparse::KernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.Near() );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( ker.OffDiag().data() );
                }

                ptic("Reduce NF_Accumulators");
                ker.Diag() = bct.GetS().NF_Accumulator().template AddReduce<Real,LInt>();
                ptoc("Reduce NF_Accumulators");
            }

            ptoc(ClassName()+"::NF_Compute");
        }
        
            
        template< typename T1, typename T2 >
        void FF_Compute( const T1 q_half_, const T1 p_half_ )
        {
            ptic(ClassName()+"::FF_Compute");
            
            using Kernel_T = TP0_Kernel_FF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, T1, T2,
                BlockClusterTree_T::IsSymmetric(),
                energy_flag, diff_flag, metric_flag
            >;
            
            Kernel_T ker ( conf, q_half_, p_half_ );
            
            FMM_Traversor traversor ( bct.Far(), ker );
            
            en += traversor.Compute();
            
            if constexpr ( metric_flag )
            {
                if constexpr ( is_symmetric )
                {
                    Sparse::KernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.Far() );
                    
                    matrix.FillLowerTriangleFromUpperTriangle( ker.OffDiag().data() );
                }
                
                ptic("Reduce FF_Accumulators");
                ker.Diag()= bct.GetS().FF_Accumulator().template AddReduce<Real,LInt>();
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
                const Int rhs_count = bct.GetS().BufferDimension() / Kernel_Block_Mul_T::ROWS;
                
                if( NF_flag && (bct.Near().NonzeroCount() > 0) )
                {
                    NF_MultiplyMetric(rhs_count);
                }
                else
                {
                    bct.GetS().PrimitiveOutputBuffer().SetZero();
                }
                
                if( VF_flag && (bct.VeryNear().NonzeroCount() > 0) )
                {
                    VF_MultiplyMetric(rhs_count);
                }
                
                
                if( FF_flag && (bct.Far().NonzeroCount() > 0) )
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
            Sparse::KernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.Near() );
            
            matrix.Dot(
                metric_values["NF"].data(),
                one,  bct.GetT().PrimitiveInputBuffer().data(),
                zero, bct.GetS().PrimitiveOutputBuffer().data(),
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
                    one, bct.GetT().PrimitiveInputBuffer().data(),
                    one, bct.GetS().PrimitiveOutputBuffer().data(),
                    rhs_count
                );
            }
        }
        
        void VF_MultiplyMetric( const Int rhs_count ) const
        {
            Sparse::KernelMatrixCSR<Kernel_Block_MulAdd_T> matrix ( bct.VeryNear() );
            
            matrix.Dot(
                metric_values["VF"].data(),
                one, bct.GetT().PrimitiveInputBuffer().data(),
                one, bct.GetS().PrimitiveOutputBuffer().data(),
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
                    one, bct.GetT().PrimitiveInputBuffer().data(),
                    one, bct.GetS().PrimitiveOutputBuffer().data(),
                    rhs_count
                );
            }
        }
        
        void FF_MultiplyMetric( const Int rhs_count ) const
        {
            Sparse::KernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.Far() );
            
            matrix.Dot(
                metric_values["FF"].data(),
                one,  bct.GetT().ClusterInputBuffer().data(),
                zero, bct.GetS().ClusterOutputBuffer().data(),
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
                    one, bct.GetT().ClusterInputBuffer().data(),
                    one, bct.GetS().ClusterOutputBuffer().data(),
                    rhs_count
                );
            }
        }
        
        
    public:
        
        std::string ClassName() const
        {
            return "TP0_Traversor<"
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
