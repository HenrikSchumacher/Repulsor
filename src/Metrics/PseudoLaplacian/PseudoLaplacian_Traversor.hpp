#pragma once

#include "PseudoLaplacian_Kernel_FF.hpp"
#include "PseudoLaplacian_Kernel_NF.hpp"
#include "PseudoLaplacian_Kernel_VF.hpp"

namespace Repulsor
{
    using namespace FMM;
    
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename BlockClusterTree_T_,
        bool high_order
    >
    class PseudoLaplacian_Traversor
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

        using Kernel_Block_MulAdd_T = ScalarBlockKernel_fixed<
            high_order ? AMB_DIM : 1, high_order ? AMB_DIM : 1, AMB_DIM, true,
            Real, Real, Real, Int, LInt,
            1,    1,
            true, true, false, true,
            true, true,
            false
        >;
        
        using Kernel_Block_Mul_T = ScalarBlockKernel_fixed<
            high_order ? AMB_DIM : 1, high_order ? AMB_DIM : 1, AMB_DIM, true,
            Real, Real, Real, Int, LInt,
            1,    0,
            true, true, false, true,
            true, true,
            false
        >;
        
        using Kernel_Diag_MulAdd_T = ScalarBlockKernel_fixed<
            high_order ? AMB_DIM : 1, high_order ? AMB_DIM : 1, AMB_DIM, true,
            Real, Real, Real, Int, LInt,
            1,    0,
            true, true, false, true,
            true, true,
            false
        >;

    
        static constexpr Int VF_blk_size = Kernel_Block_MulAdd_T::ROWS * Kernel_Block_MulAdd_T::COLS;
        static constexpr Int NF_blk_size =    Kernel_Block_Mul_T::ROWS *    Kernel_Block_Mul_T::COLS;
        static constexpr Int FF_blk_size =    Kernel_Block_Mul_T::ROWS *    Kernel_Block_Mul_T::COLS;
        
    protected:
        
        static constexpr Real zero = 0;
        static constexpr Real one  = 1;
        static constexpr Real two  = 2;
        
    
    public:
        
        PseudoLaplacian_Traversor() = delete;
        
        PseudoLaplacian_Traversor(
              const BlockClusterTree_T & bct_,
              ValueContainer_T & metric_values_,
              const Real s_
        )
        :   bct           ( bct_          )
        ,   conf          ( bct.GetS(), bct.GetT(), metric_values_ )
        ,   metric_values ( metric_values_ )
        ,   s             ( s_            )
        {}

        ~PseudoLaplacian_Traversor() = default;
        

    protected:

        const BlockClusterTree_T & bct;
        
        Configurator_T conf;
        
        ValueContainer_T & metric_values;

        const Real s;
        
    public:
        

//##############################################################################################
//      Compute
//##############################################################################################
        
        Real Compute()
        {
            const Int thread_count = bct.ThreadCount();
            
            DummyAllocators();

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
            
            VF_Compute();
            NF_Compute();
            FF_Compute();

            DummyAllocators();
            
            return 0;
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
        
        void VF_Compute()
        {
            ptic(ClassName()+"::VF_Compute");
            
            using Kernel_T = PseudoLaplacian_Kernel_VF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, BlockClusterTree_T::is_symmetric,
                high_order
            >;
            
            Kernel_T ker ( conf, bct.NearFieldSeparationParameter(), 20, s );

            FMM_Traversor traversor ( bct.VeryNear(), ker );
            
            (void)traversor.Compute();

            if constexpr ( is_symmetric )
            {
                Sparse::KernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.VeryNear() );
                
                matrix.FillLowerTriangleFromUpperTriangle( ker.OffDiag().data() );
            }
            
            ptic("Reduce VF_Accumulators");
            ker.Diag() = bct.GetS().VF_Accumulator().template AddReduce<Real,LInt>();
            ptoc("Reduce VF_Accumulators");
            
            ptoc(ClassName()+"::VF_Compute");
        }
            
        
        void NF_Compute()
        {
            ptic(ClassName()+"::NF_Compute");

            using Kernel_T = PseudoLaplacian_Kernel_NF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, BlockClusterTree_T::is_symmetric,
                high_order
            >;

            Kernel_T ker ( conf, s );

            FMM_Traversor traversor ( bct.Near(), ker );

            (void)traversor.Compute();

            if constexpr ( is_symmetric )
            {
                Sparse::KernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.Near() );
                
                matrix.FillLowerTriangleFromUpperTriangle( ker.OffDiag().data() );
            }

            ptic("Reduce NF_Accumulators");
            ker.Diag() = bct.GetS().NF_Accumulator().template AddReduce<Real,LInt>();
            ptoc("Reduce NF_Accumulators");

            ptoc(ClassName()+"::NF_Compute");
        }
        
        
        void FF_Compute()
        {
            ptic(ClassName()+"::FF_Compute");
            
            using Kernel_T = PseudoLaplacian_Kernel_FF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, BlockClusterTree_T::is_symmetric,
                high_order
            >;
            
            Kernel_T ker ( conf, s );
            
            FMM_Traversor traversor ( bct.Far(), ker );
            
            (void)traversor.Compute();
            
            if constexpr ( is_symmetric )
            {
                Sparse::KernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.Far() );
                
                matrix.FillLowerTriangleFromUpperTriangle( ker.OffDiag().data() );
            }
            
            ptic("Reduce FF_Accumulators");
            ker.Diag()= bct.GetS().FF_Accumulator().template AddReduce<Real,LInt>();
            ptoc("Reduce FF_Accumulators");
            
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
            return "PseudoLaplacian_Traversor<"
                + ToString(S_DOM_DIM)   + ","
                + ToString(T_DOM_DIM)   + ","
                + bct.ClassName()       + ","
                + ToString(high_order) +
                + ">";
        }

    }; // class PseudoLaplacian_Traversor
    
}// namespace Repulsor
