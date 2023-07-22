#pragma once

namespace Repulsor
{
    using namespace FMM;
    
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename ClusterTree_T_,
        bool symmetricQ_,
        bool energy_flag, bool diff_flag, bool metric_flag
    >
    class TP_AllPairs_Traversor
    {
    public:
        

        using ClusterTree_T      = ClusterTree_T_;
        
        using Real               = typename ClusterTree_T::Real;
        using Int                = typename ClusterTree_T::Int;
        using SReal              = typename ClusterTree_T::SReal;
        using ExtReal            = typename ClusterTree_T::ExtReal;
        
        using Configurator_T     = FMM_Configurator<ClusterTree_T>;
        using ValueContainer_T   = typename Configurator_T::ValueContainer_T;
        using LInt               = typename Configurator_T::LInt;
        
        using Accumulator_T      = typename ClusterTree_T::Accumulator_T;
        
        static constexpr Int  AMB_DIM      = ClusterTree_T::AMB_DIM;
        static constexpr Int  S_DOM_DIM    = S_DOM_DIM_;
        static constexpr Int  T_DOM_DIM    = T_DOM_DIM_;
        static constexpr bool symmetricQ = symmetricQ_;
        
        using TangentVector_T   = Tensor2<ExtReal,Int>;
        using CotangentVector_T = Tensor2<ExtReal,Int>;
        
        
        using Kernel_Block_MulAdd_T = ArrowheadBlockKernel_fixed<
            AMB_DIM+1, AMB_DIM+1, AMB_DIM,
            Real, Real, Real, Int, LInt,
            Scalar::Flag::Plus, Scalar::Flag::Plus, // CAUTION: We use add-in instead of overwrite!
            true,
            true, true,
            true
        >;
        
        using Kernel_Block_Mul_T = ArrowheadBlockKernel_fixed<
            AMB_DIM+1, AMB_DIM+1, AMB_DIM,
            Real, Real, Real, Int, LInt,
            Scalar::Flag::Plus, Scalar::Flag::Zero,
            true,
            true, true,
            true
        >;
  
        using Kernel_Diag_MulAdd_T = DenseBlockKernel_fixed<
            AMB_DIM+1, AMB_DIM+1, AMB_DIM,
            Real, Real, Real, Int, LInt,
            Scalar::Flag::Plus, Scalar::Flag::Plus,
            true, false, true,
            true, true,
            true
        >;
        
//        using Kernel_Diag_MulAdd_T = DenseBlockKernel_fixed<
//            AMB_DIM+1, AMB_DIM+1, AMB_DIM, true,
//            Real, Real, Real, Int, LInt,
//            Scalar::Flag::Plus, Scalar::Flag::Plus,
//            true, false, true,
//            true, false, true, true,
//            true, false,
//            1   , 2,
//            true
//        >;

        static constexpr Int NF_blk_size =    Kernel_Block_Mul_T::ROWS *    Kernel_Block_Mul_T::COLS;
        
    protected:
        
        static constexpr Real zero = 0;
        static constexpr Real one  = 1;
        static constexpr Real two  = 2;
        
    
    public:
        
        TP_AllPairs_Traversor() = delete;
        
        TP_AllPairs_Traversor(
              cref<ClusterTree_T> S_,
              cref<ClusterTree_T> T_,
              mref<ValueContainer_T> metric_values_,
              const Real q_,
              const Real p_
        )
        :   S             ( S_          )
        ,   T             ( T_          )
        ,   conf          ( S, T, metric_values_ )
        ,   metric_values ( metric_values_ )
        ,   q             ( q_            )
        ,   p             ( p_            )
        ,   q_half_real   ( q / two       )
        ,   p_half_real   ( p / two       )
        ,   q_half_int    ( static_cast<Int>(q_half_real) )
        ,   p_half_int    ( static_cast<Int>(p_half_real) )
        ,   q_half_is_int ( q_half_real == static_cast<Real>(q_half_int) )
        ,   p_half_is_int ( p_half_real == static_cast<Real>(p_half_int) )
        {}

        ~TP_AllPairs_Traversor() = default;
        

    protected:

        cref<ClusterTree_T> S;
        cref<ClusterTree_T> T;
        
        Configurator_T conf;
        
        mref<ValueContainer_T> metric_values;
        
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
            const Int thread_count = std::min( S.ThreadCount(), T.ThreadCount() );
            
            en = zero;
            
            DummyAllocators();

            if constexpr ( diff_flag )
            {
                S.CleanseDerivativeBuffers();
                
                if( !symmetricQ )
                {
                    T.CleanseDerivativeBuffers();
                }
            }

            if constexpr ( metric_flag )
            {
                ptic("Allocate accumulators");
                
                S.NF_Accumulator()
                = Accumulator_T( thread_count, S.PrimitiveCount(), NF_blk_size, 0 );

                if constexpr ( !symmetricQ )
                {
                    T.NF_Accumulator()
                    = Accumulator_T( thread_count, T.GetT().PrimitiveCount(), NF_blk_size, 0 );
                }
                ptoc("Allocate accumulators");
            }
            else
            {
                DummyAllocators();
            }
            
            if( q_half_real >= one )
            {
                if( q_half_is_int )
                {
                    if( p_half_is_int )
                    {
                        NF_Compute<Int,Int,2>( q_half_int, p_half_int );
                    }
                    else
                    {
                        NF_Compute<Int,Real,2>( q_half_int, p_half_real );
                    }
                }
                else
                {
                    if( p_half_is_int)
                    {
                        NF_Compute<Real,Int,2>( q_half_real, p_half_int );
                    }
                    else
                    {
                        NF_Compute<Real,Real,2>( q_half_real, p_half_real );
                    }
                }
            }
            else if( q_half_real == zero )
            {
                if( p_half_is_int )
                {
                    NF_Compute<Int,Int,0>( q_half_int, p_half_int );
                }
                else
                {
                    NF_Compute<Int,Real,0>( q_half_int, p_half_real );
                }
            }
            else
            {
                if( q_half_real < zero )
                {
                    eprint("Tangent-point kernels not well-defined for q < 0.");
                }
                else
                {
                    eprint("Tangent-point kernels not implement for 0 < q < 2.");
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

            const Int thread_count = std::min( S.ThreadCount(), T.ThreadCount() );
            S.NF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
            if constexpr ( !symmetricQ )
            {
                T.NF_Accumulator() = Accumulator_T( thread_count, 1, 1, 0 );
            }
            ptoc("Allocate dummy accumulators");
        }
        

        
//##############################################################################################
//      Compute subroutines
//##############################################################################################
    
    protected:
        
        template< typename T1, typename T2, int q_flag>
        void NF_Compute( const T1 q_half_, const T2 p_half_ )
        {
            ptic(ClassName()+"::NF_Compute");

            using Kernel_T = TP_Kernel_NF<
                S_DOM_DIM, T_DOM_DIM,
                ClusterTree_T, T1, T2, q_flag,
                symmetricQ,
                energy_flag, diff_flag, metric_flag
            >;
            
            Kernel_T ker ( conf, 0, q_half_, p_half_ );

            AllPairs_Traversor traversor ( S.PrimitiveCount(), T.PrimitiveCount(), ker );
            
            en += traversor.Compute();
            
            if constexpr ( metric_flag )
            {
                if constexpr ( symmetricQ )
                {
                    // TODO: Need dense matrix here.
//                    SparseKernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.Near() );
//
//                    matrix.FillLowerTriangleFromUpperTriangle( ker.OffDiag().data() );
                }

                ptic("Reduce NF_Accumulators");
                ker.Diag() = S.NF_Accumulator().template AddReduce<Real,LInt>();
                ptoc("Reduce NF_Accumulators");
            }

            ptoc(ClassName()+"::NF_Compute");
        }
        
//##############################################################################################
//      Multiply metric
//##############################################################################################
        
    public:
        
        void MultiplyMetric(
            const bool NF_flag = true
        ) const
        {
            if constexpr ( metric_flag )
            {
                const Int rhs_count = S.BufferDimension() / Kernel_Block_Mul_T::ROWS;
                
                if( NF_flag && ( S.PrimitiveCount()> 0) )
                {
                    NF_MultiplyMetric(rhs_count);
                }
                else
                {
                    S.PrimitiveOutputBuffer().SetZero();
                }
            }
        }
        
//##############################################################################################
//      Multiply metric subroutines
//##############################################################################################
        
    protected:

        void NF_MultiplyMetric( const Int rhs_count ) const
        {
            // TODO: DenseKernelMatrix
//            SparseKernelMatrixCSR<Kernel_Block_Mul_T> matrix ( bct.Near() );
//
//            matrix.Dot(
//                metric_values.NF.data(),
//                one,  T.PrimitiveInputBuffer().data(),
//                zero, S.PrimitiveOutputBuffer().data(),
//                rhs_count
//            );
            
            if constexpr ( symmetricQ )
            {
                DiagonalKernelMatrix<Kernel_Diag_MulAdd_T> diag (
                    S.PrimitiveCount(),
                    S.ThreadCount()
                );

                diag.Dot(
                    metric_values.NF_diag.data(),
                    one, T.PrimitiveInputBuffer().data(),
                    one, S.PrimitiveOutputBuffer().data(),
                    rhs_count
                );
            }
        }
        
    public:
        
        std::string ClassName() const
        {
            return "TP_AllPairs_Traversor<"
                + ToString(S_DOM_DIM)   + ","
                + ToString(T_DOM_DIM)   + ","
                + S.ClassName()       + ","
                + ToString(energy_flag) + ","
                + ToString(diff_flag)   + ","
                + ToString(metric_flag) +
                + ">";
        }
        
    }; // class TP_AllPairs_Traversor
    
}// namespace Repulsor
