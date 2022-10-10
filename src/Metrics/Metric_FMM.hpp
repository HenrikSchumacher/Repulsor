#pragma once

#define CLASS Metric_FMM
#define BASE  MetricBase<Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, int DEGREE, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:
        
        using BlockClusterTree_T = BlockClusterTree <AMB_DIM,Real,Int,SReal,ExtReal>;
        using A_Kernel_T         = Metric__NFK_Adaptive <DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>;
        using N_Kernel_T         = Metric__NFK          <DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>;
        using F_Kernel_T         = Metric__FFK_FMM  <AMB_DIM,DEGREE,Real,Int>;
        
        CLASS() : BASE() {}
        
        CLASS(
            const BlockClusterTree_T & bct_,
            const A_Kernel_T & A,
            const N_Kernel_T & N,
            const F_Kernel_T & F,
            const ExtReal weight_ = static_cast<ExtReal>(1)
        )
        :   BASE  ( weight_   )
        ,   bct   ( &bct_     )
        ,   A_ker ( A.Clone() )
        ,   N_ker ( N.Clone() )
        ,   F_ker ( F.Clone() )
        {}
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   BASE  ( other.weight         )
        ,   bct   ( other.bct            )
        ,   A_ker ( other.A_ker->Clone() )
        ,   N_ker ( other.N_ker->Clone() )
        ,   F_ker ( other.F_ker->Clone() )
        {}

        virtual ~CLASS() override = default;
        
    protected:
        
        using BASE::weight;
        
        const BlockClusterTree_T * bct = nullptr;
        
        std::unique_ptr<A_Kernel_T> A_ker;
        std::unique_ptr<N_Kernel_T> N_ker;
        std::unique_ptr<F_Kernel_T> F_ker;

        Tensor1<Real,Int> verynear_values;
        Tensor1<Real,Int>     near_values;
        Tensor1<Real,Int>      far_values;
        
        mutable bool metrics_initialized = false;
        
    public:

//        virtual Int VeryNearValueSize() const = 0;
//
//        virtual Int NearValueSize() const = 0;
//
//        virtual Int FarValueSize() const = 0;
        
        virtual const Tensor1<Real,Int> & VeryNearFieldValues() const override
        {
            return verynear_values;
        }
        
        virtual const Tensor1<Real,Int> & NearFieldValues() const override
        {
            return near_values;
        }

        virtual const Tensor1<Real,Int> & FarFieldValues() const override
        {
            return far_values;
        }

        void RequireMetrics() const override
        {
            if( metrics_initialized )
            {
                return;
            }
            
            ptic(ClassName()+"::RequireMetrics");

            VeryNearField();
            
            NearField();
            
            FarField ();
            
            if( bct->IsSymmetric() )
            {
                Far_FillLowerTriangleFromUpperTriangle();
            }
            
            ComputeDiagonals();
            
            metrics_initialized = true;
                
            ptoc(ClassName()+"::RequireMetrics");
            
        } // RequireMetrics
        
    protected:
        

//##############################################################################################
//      Compute diagonals
//##############################################################################################
        
    protected:
        
        virtual void ComputeDiagonals() const
        {
//            ptic(ClassName()+"::ComputeDiagonals");
//
//            Int cols = 1;
//
//            auto & T = bct->GetT();
//            auto & S = bct->GetS();
//
//            T.RequireBuffers(cols);
//            S.RequireBuffers(cols);
//
////            S.CleanseBuffers();
////            T.CleanseBuffers();
//
////            S.PrimitiveOutputBuffer().Fill( static_cast<Real>(0) );
//
//            {
//                      Real * restrict const in_buffer = T.PrimitiveInputBuffer().data();
//                const Real * restrict const near_data = T.PrimitiveNearFieldData().data();
//                const Int near_dim = T.NearDim();
//                const Int primitive_count = T.PrimitiveCount();
//
//                for( Int i = 0; i < primitive_count; ++i )
//                {
//                    in_buffer[i] = near_data[ near_dim * i ];
//                }
//            }
//
//            T.PrimitivesToClusters(false);
//
//            T.PercolateUp();
//
//            for( auto const & f : near_values)
//            {
//                ApplyNearFieldKernel( static_cast<Real>(1), f.first );
//
//                ApplyFarFieldKernel ( static_cast<Real>(1), f.first );
//
//                S.PercolateDown();
//
//                S.ClustersToPrimitives( true );
//
//                Real const * restrict const u = S.PrimitiveOutputBuffer().data();
//                Real const * restrict const a = S.PrimitiveNearFieldData().data();
//                Real       * restrict const values = near_values.find(f.first)->second.data();
//                Int  const * restrict const diag_ptr = bct->Near().Diag().data();
//
//                const Int m        = S.PrimitiveCount();
//                const Int near_dim = S.NearDim();
//
//                #pragma omp parallel for
//                for( Int i = 0; i < m; ++i )
//                {
//                    values[diag_ptr[i]] -= u[i]/a[ near_dim * i ];
//                }
//            }
            
            ptoc(ClassName()+"::ComputeDiagonals");
            
        } // ComputeDiagonals
        
        
        
        
    public:
        

        
        virtual void ApplyFarFieldKernel( const Real alpha ) const override
        {
            ptic(ClassName()+"::ApplyFarFieldKernel");

            const Int  * restrict const rp = bct->Far().Outer().data();
            const Int  * restrict const ci = bct->Far().Inner().data();
            const Real * restrict const a  = far_values.data();
            const Real * restrict const X  = bct->GetT().ClusterInputBuffer().data();
                  Real * restrict const Y  = bct->GetT().ClusterOutputBuffer().data();

            
            const auto & job_ptr = bct->Far().JobPtr();
            
            #pragma omp parallel for num_threads( job_ptr.Size()-1 )
            for( Int thread = 0; thread < job_ptr.Size()-1; ++thread )
            {
                std::unique_ptr<F_Kernel_T> F = F_ker->Clone();
                
                F->LoadInputBuffer(X);
                F->LoadOutputBuffer(Y);
                
                const Int i_begin = job_ptr[thread  ];
                const Int i_end   = job_ptr[thread+1];
                
                for( Int i = i_begin; i < i_end; ++i )
                {
                    const Int l_begin = rp[i  ];
                    const Int l_end   = rp[i+1];
                    
                    __builtin_prefetch( &ci[l_end] );
                    
                    __builtin_prefetch( &a[l_end] );
                    
                    F->ClearVector();

                    for( Int l = l_begin; l < l_end-1; ++l )
                    {
                        F->ApplyBlock( alpha, l, ci[l] );
                    }
                    
                    // incorporate the local updates into Y-buffer
                    F->WriteVector(i);
                }
            }
            
            ptoc(ClassName()+"::ApplyFarFieldKernel");
        }
        
        void Multiply_DenseMatrix(
            const Real    alpha, const ExtReal * X,
            const ExtReal  beta,       ExtReal * Y,
            const Int cols
        ) const override
        {
            ptic(ClassName()+"::Multiply_DenseMatrix");
            
            if( bct == nullptr )
            {
                eprint(ClassName()+"::Multiply_DenseMatrix: BlockClusterTree not initialized.");
                ptoc(ClassName()+"::Multiply_DenseMatrix");
                return;
            }
                        
            auto & S = bct->GetS();
            auto & T = bct->GetT();
            
            T.Pre( X, cols, KernelType::MixedOrder );
            
            const Int buffer_dim = T.BufferDimension();

            S.RequireBuffers(buffer_dim); // Tell the S-side what it has to expect.
            
            ApplyVeryNearFieldKernel( weight * alpha );
            
            ApplyNearFieldKernel    ( weight * alpha );
            
            ApplyFarFieldKernel     ( weight * alpha );

            S.Post( Y, KernelType::MixedOrder, beta );
            
            ptoc(ClassName()+"::Multiply_DenseMatrix");
        }

        
//        void TestDiagonal() const override
//        {
//            ptic(ClassName()+"::TestDiagonal");
//
//            Int cols = 1;
//
//            auto & T = bct->GetT();
//            auto & S = bct->GetS();
//
//            T.RequireBuffers(cols);
//            S.RequireBuffers(cols);
//
//            S.CleanseBuffers();
//            T.CleanseBuffers();
//
//            S.PrimitiveOutputBuffer().Fill( static_cast<Real>(0) );
//
//            for( Int i = 0; i < T.PrimitiveCount(); ++i )
//            {
//                T.PrimitiveInputBuffer()[i] = T.PrimitiveNearFieldData()(i,0);
//            }
//
//            T.PrimitivesToClusters(false);
//
//            T.PercolateUp();
//
//            for( auto const & f : near_values)
//            {
//                ApplyNearFieldKernel( static_cast<Real>(1), f.first );
//
//                ApplyFarFieldKernel ( static_cast<Real>(1), f.first );
//
//                S.PercolateDown();
//
//                S.ClustersToPrimitives(true);
//
//                Real maxnorm = static_cast<Real>(0);
//
//                for( Int i = 0 ; i < S.PrimitiveCount(); ++i)
//                {
//                    maxnorm = std::max( maxnorm, S.PrimitiveOutputBuffer()[i] );
//                }
//
//                valprint("maxnorm", maxnorm);
//                DUMP(maxnorm);
//
//            }
//
//            ptoc(ClassName()+"::TestDiagonal");
//        }
        
    

        
        
        
        void Far_FillLowerTriangleFromUpperTriangle() const 
        {
            ptic(ClassName()+"::Far_FillLowerTriangleFromUpperTriangle");
            
            // FarField
            if( bct->Far().WellFormed() )
            {
                const Int * restrict const diag   = bct->Far().Diag().data();
                const Int * restrict const outer  = bct->Far().Outer().data();
                const Int * restrict const inner  = bct->Far().Inner().data();
                
                const auto & job_ptr = bct->Far().LowerTriangularJobPtr();
                
                const Int thread_count = job_ptr.Size()-1;
                
                #pragma omp parallel for num_threads( thread_count )
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    auto F = F_ker->Clone();
                    
                    F->LoadValueBuffers( far_values );
                    
                    const Int i_begin = job_ptr[thread];
                    const Int i_end   = job_ptr[thread+1];
                    
                    for( Int i = i_begin; i < i_end; ++i )
                    {
                        const Int k_begin = outer[i];
                        const Int k_end   =  diag[i];
                        
                        for( Int k = k_begin; k < k_end; ++k )
                        {
                            const Int j = inner[k];
                            
                            Int L =  diag[j];
                            Int R = outer[j+1]-1;
                            
                            while( L < R )
                            {
                                const Int M   = R - (R-L)/static_cast<Int>(2);
                                const Int col = inner[M];

                                if( col > i )
                                {
                                    R = M-1;
                                }
                                else
                                {
                                    L = M;
                                }
                            }
                            
                            F->TransposeBlock(L, k);
                            
                        } // for( Int k = k_begin; k < k_end; ++k )

                    } // for( Int i = i_begin; i < i_end; ++i )
                    
                } // #pragma omp parallel
            }
            
            ptoc(ClassName()+"::Far_FillLowerTriangleFromUpperTriangle");
        }
        
//##############################################################################################
//      Matrix multiplication
//##############################################################################################

        
    virtual void ApplyKernel(
        const SparseBinaryMatrixCSR<Int> matrix,
        const Metric_Kernel<Real,Int> & kernel,
        const Real * restrict const a,
        const Real alpha,
        const Real * restrict const X,
              Real * restrict const Y
    ) const
    {
        const auto & matrix = bct->Near();
        
        const Int  * restrict const rp = matrix.Outer().data();
        const Int  * restrict const ci = matrix.Inner().data();
        
        const auto & job_ptr = matrix.JobPtr();
        
        const Int thread_count = job_ptr.Size()-1;
        
        #pragma omp parallel for num_threads( thread_count )
        for( Int thread = 0; thread < thread_count; ++thread )
        {
            std::unique_ptr<N_Kernel_T> N = N_ker->Clone();
            
            N->LoadInputBuffer(X);
            N->LoadOutputBuffer(Y);
            
            const Int i_begin = job_ptr[thread  ];
            const Int i_end   = job_ptr[thread+1];
            
            for( Int i = i_begin; i < i_end; ++i )
            {
                const Int l_begin = rp[i  ];
                const Int l_end   = rp[i+1];
                
                __builtin_prefetch( &ci[l_end] );
                
                __builtin_prefetch( &a[l_end] );
                
                N->ClearVector();

                for( Int l = l_begin; l < l_end-1; ++l )
                {
                    N->ApplyBlock( alpha, l, ci[l] );
                }
                
                // incorporate the local updates into Y-buffer
                N->WriteVector(i);
            }
        }
        
        ptoc(ClassName()+"::ApplyNearFieldKernel" );
    }
        
    public:
        
        virtual std::string Stats() const override
        {
            std::stringstream s;
            
            s << ClassName() << ": \n\n";
            if( N_ker != nullptr )
            {
                s << "          NearFieldKernel = " << N_ker->Stats() << "\n\n";
            }
            
            if( F_ker != nullptr )
            {
                s << "          FarFieldKernel = " << F_ker->Stats() << "\n";
            }
            return s.str();
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+ToString(DEGREE)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    };
    
}// namespace Repulsor

#undef BASE
#undef CLASS
