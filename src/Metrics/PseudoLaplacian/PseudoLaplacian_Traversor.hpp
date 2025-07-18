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
        using ValueContainer_T   = typename Configurator_T::ValueContainer_T;
        using Values_T           = typename ValueContainer_T::Values_T;
        
        using Accumulator_T      = typename ClusterTree_T::Accumulator_T;
        
        static constexpr Int  AMB_DIM      = ClusterTree_T::AMB_DIM;
        static constexpr Int  S_DOM_DIM    = S_DOM_DIM_;
        static constexpr Int  T_DOM_DIM    = T_DOM_DIM_;
        static constexpr bool symmetricQ = BlockClusterTree_T::SymmetricQ();
        
        using TangentVector_T   = Tensor2<ExtReal,Int>;
        using CotangentVector_T = Tensor2<ExtReal,Int>;

        using Kernel_Block_MulAdd_T = ScalarBlockKernel_Tiny<
            high_order ? AMB_DIM : 1, high_order ? AMB_DIM : 1, AMB_DIM,
            Real, Real, Real, Int, LInt,
            Scalar::Flag::Plus, Scalar::Flag::Plus,
            true, true,
            true
        >;
        
        using Kernel_Block_Mul_T = ScalarBlockKernel_Tiny<
            high_order ? AMB_DIM : 1, high_order ? AMB_DIM : 1, AMB_DIM,
            Real, Real, Real, Int, LInt,
            Scalar::Flag::Plus, Scalar::Flag::Zero,
            true, true,
            true
        >;
        
        using Kernel_Diag_MulAdd_T = ScalarBlockKernel_Tiny<
            high_order ? AMB_DIM : 1, high_order ? AMB_DIM : 1, AMB_DIM,
            Real, Real, Real, Int, LInt,
            Scalar::Flag::Plus, Scalar::Flag::Plus,
            true, true,
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
        
        PseudoLaplacian_Traversor(
              cref<BlockClusterTree_T> bct_,
              mref<ValueContainer_T>   metric_values_,
              const Real s_
        )
        :   bct           ( bct_           )
        ,   conf          ( bct.GetS(), bct.GetT(), metric_values_ )
        ,   metric_values ( metric_values_ )
        ,   s             ( s_             )
        {}

        // Default constructor
        PseudoLaplacian_Traversor() = delete;
        // Destructor
        virtual ~PseudoLaplacian_Traversor() = default;
        // Copy constructor
        PseudoLaplacian_Traversor( const PseudoLaplacian_Traversor & other ) = default;
        // Copy assignment operator
        PseudoLaplacian_Traversor & operator=( const PseudoLaplacian_Traversor & other ) = default;
        // Move constructor
        PseudoLaplacian_Traversor( PseudoLaplacian_Traversor && other ) = default;
        // Move assignment operator
        PseudoLaplacian_Traversor & operator=( PseudoLaplacian_Traversor && other ) = default;
        

    protected:

        cref<BlockClusterTree_T> bct;
        
        Configurator_T conf;
        
        mref<ValueContainer_T> metric_values;

        const Real s;
        
    public:
        

//###########################################################
//      Compute
//###########################################################
        
        Real Compute()
        {
            const Int thread_count = bct.ThreadCount();
            
            bct.GetS().VF_Accumulator()
            = Accumulator_T( thread_count, bct.GetS().PrimitiveCount(), VF_blk_size, 0 );
            bct.GetS().NF_Accumulator()
            = Accumulator_T( thread_count, bct.GetS().PrimitiveCount(), NF_blk_size, 0 );
            bct.GetS().FF_Accumulator()
            = Accumulator_T( thread_count, bct.GetS().ClusterCount(),   FF_blk_size, 0 );

            if constexpr ( !symmetricQ )
            {
                bct.GetT().VF_Accumulator()
                = Accumulator_T( thread_count, bct.GetT().PrimitiveCount(), VF_blk_size, 0 );
                bct.GetT().NF_Accumulator()
                = Accumulator_T( thread_count, bct.GetT().PrimitiveCount(), NF_blk_size, 0 );
                bct.GetT().FF_Accumulator()
                = Accumulator_T( thread_count, bct.GetT().ClusterCount(),   FF_blk_size, 0 );
            }
            
            VF_Compute();
            NF_Compute();
            FF_Compute();

            AllocateDummyAccumulator();
            
            return 0;
        }
     
        void AllocateDummyAccumulator()
        {
            const Int thread_count = bct.ThreadCount();

            bct.GetS().VF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
            bct.GetS().NF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
            bct.GetS().FF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );

            if constexpr ( !symmetricQ )
            {
                bct.GetT().VF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
                bct.GetT().NF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
                bct.GetT().FF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
            }
        }
        

        
//###########################################################
//      Compute subroutines
//###########################################################
    
    protected:
        
        void VF_Compute()
        {
            if( bct.VeryNear().NonzeroCount() == 0 )
            {
                return;
            }
            
            TOOLS_PTIC(ClassName()+"::VF_Compute");
            
            using Kernel_T = PseudoLaplacian_Kernel_VF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, BlockClusterTree_T::symmetricQ,
                high_order
            >;
            
            Kernel_T ker ( conf, 0, bct.NearFieldSeparationParameter(), 20, s );

            FMM_Traversor traversor ( bct.VeryNear(), ker );
            
            (void)traversor.Compute();

            if constexpr ( symmetricQ )
            {
                Sparse::KernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.VeryNear() );
                
                matrix.FillLowerTriangleFromUpperTriangle( ker.OffDiag().data() );
            }
            
            ker.Diag() = bct.GetS().VF_Accumulator().template AddReduce<Real,LInt>();
            
            TOOLS_PTOC(ClassName()+"::VF_Compute");
        }
            
        
        void NF_Compute()
        {
            if( bct.Near().NonzeroCount() == 0 )
            {
                return;
            }
            
            TOOLS_PTIC(ClassName()+"::NF_Compute");

            using Kernel_T = PseudoLaplacian_Kernel_NF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, BlockClusterTree_T::symmetricQ,
                high_order
            >;

            Kernel_T ker ( conf, 0, s );

            FMM_Traversor traversor ( bct.Near(), ker );

            (void)traversor.Compute();

            if constexpr ( symmetricQ )
            {
                Sparse::KernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.Near() );
                
                matrix.FillLowerTriangleFromUpperTriangle( ker.OffDiag().data() );
            }

            ker.Diag() = bct.GetS().NF_Accumulator().template AddReduce<Real,LInt>();

            TOOLS_PTOC(ClassName()+"::NF_Compute");
        }
        
        
        void FF_Compute()
        {
            if( bct.Far().NonzeroCount() == 0 )
            {
                return;
            }
            
            TOOLS_PTIC(ClassName()+"::FF_Compute");
            
            using Kernel_T = PseudoLaplacian_Kernel_FF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, BlockClusterTree_T::symmetricQ,
                high_order
            >;
            
            Kernel_T ker ( conf, 0, s );
            
            FMM_Traversor traversor ( bct.Far(), ker );
            
            (void)traversor.Compute();
            
            if constexpr ( symmetricQ )
            {
                Sparse::KernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.Far() );
                
                matrix.FillLowerTriangleFromUpperTriangle( ker.OffDiag().data() );
            }
            
            ker.Diag()= bct.GetS().FF_Accumulator().template AddReduce<Real,LInt>();
            
            TOOLS_PTOC(ClassName()+"::FF_Compute");
        }
        
//###########################################################
//      Multiply metric
//###########################################################
        
    public:
        
        void MultiplyMetric(
            const bool VF_flag = true,
            const bool NF_flag = true,
            const bool FF_flag = true
        ) const
        {
            const Int nrhs = bct.GetS().BufferDim() / Kernel_Block_Mul_T::ROWS;
            
            if( NF_flag && (bct.Near().NonzeroCount() > 0) )
            {
                NF_MultiplyMetric(nrhs);
            }
            else
            {
                bct.GetS().PrimitiveOutputBuffer().SetZero();
            }
            
            if( VF_flag && (bct.VeryNear().NonzeroCount() > 0) )
            {
                VF_MultiplyMetric(nrhs);
            }
            
            
            if( FF_flag && (bct.Far().NonzeroCount() > 0) )
            {
                FF_MultiplyMetric(nrhs);
            }
            else
            {
                bct.GetS().ClusterOutputBuffer().SetZero();
            }
        }
        
//###########################################################
//      Multiply metric subroutines
//###########################################################
        
    protected:

        void NF_MultiplyMetric( const Int nrhs ) const
        {
            Sparse::KernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.Near() );
        
            matrix.Dot(
                metric_values.NF.data(),
                one,  bct.GetT().PrimitiveInputBuffer().data(),
                zero, bct.GetS().PrimitiveOutputBuffer().data(),
                nrhs
            );
            
            if constexpr ( symmetricQ )
            {
                DiagonalKernelMatrix<Kernel_Diag_MulAdd_T> diag (
                    bct.Near().RowCount(),
                    bct.Near().ThreadCount()
                );

                diag.Dot(
                    metric_values.NF_diag.data(),
                    one, bct.GetT().PrimitiveInputBuffer().data(),
                    one, bct.GetS().PrimitiveOutputBuffer().data(),
                    nrhs
                );
            }
        }
        
        void VF_MultiplyMetric( const Int nrhs ) const
        {
            Sparse::KernelMatrixCSR<Kernel_Block_MulAdd_T> matrix ( bct.VeryNear() );
            
            matrix.Dot(
                metric_values.VF.data(),
                one, bct.GetT().PrimitiveInputBuffer().data(),
                one, bct.GetS().PrimitiveOutputBuffer().data(),
                nrhs
            );
            
            if constexpr ( symmetricQ )
            {
                DiagonalKernelMatrix<Kernel_Diag_MulAdd_T> diag (
                    bct.VeryNear().RowCount(),
                    bct.VeryNear().ThreadCount()
                );

                diag.Dot(
                    metric_values.VF_diag.data(),
                    one, bct.GetT().PrimitiveInputBuffer().data(),
                    one, bct.GetS().PrimitiveOutputBuffer().data(),
                    nrhs
                );
            }
        }
        
        void FF_MultiplyMetric( const Int nrhs ) const
        {
            Sparse::KernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.Far() );
            
            matrix.Dot(
                metric_values.FF.data(),
                one,  bct.GetT().ClusterInputBuffer().data(),
                zero, bct.GetS().ClusterOutputBuffer().data(),
                nrhs
           );
            
            if constexpr ( symmetricQ )
            {
                DiagonalKernelMatrix<Kernel_Diag_MulAdd_T> diag (
                    bct.Far().RowCount(),
                    bct.Far().ThreadCount()
                );
                
                diag.Dot(
                    metric_values.FF_diag.data(),
                    one, bct.GetT().ClusterInputBuffer().data(),
                    one, bct.GetS().ClusterOutputBuffer().data(),
                    nrhs
                );
            }
        }
        
        
    public:
        
        std::string ClassName() const
        {
            return "PseudoLaplacian_Traversor<"
                + ToString(S_DOM_DIM)   + ","
                + ToString(T_DOM_DIM)   + ","
//                + bct.ClassName()       + ","
                + "...,"
                + ToString(high_order) +
                + ">";
        }

    }; // class PseudoLaplacian_Traversor
    
}// namespace Repulsor

