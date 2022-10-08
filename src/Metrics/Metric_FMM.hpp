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
        using N_Kernel_T         = Metric__NFK      <DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>;
        using F_Kernel_T         = Metric__FFK_FMM  <AMB_DIM,DEGREE,Real,Int>;
        
        CLASS() : BASE() {}
        
        CLASS(
            const BlockClusterTree_T & bct_,
            const N_Kernel_T & N,
            const F_Kernel_T & F,
            const ExtReal weight_ = static_cast<ExtReal>(1)
        )
        :   BASE  ( weight_   )
        ,   bct   ( &bct_     )
        ,   N_ker ( N.Clone() )
        ,   F_ker ( F.Clone() )
        {}
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   BASE  ( other.weight         )
        ,   bct   ( other.bct            )
        ,   N_ker ( other.N_ker->Clone() )
        ,   F_ker ( other.F_ker->Clone() )
        {}

        virtual ~CLASS() override = default;
        
    protected:
        
        using BASE::weight;
        
        const BlockClusterTree_T * bct = nullptr;
        
        std::unique_ptr<N_Kernel_T> N_ker;
        std::unique_ptr<F_Kernel_T> F_ker;
        
        mutable std::map<KernelType, Tensor1<Real,Int>> near_values;
        mutable std::map<KernelType, Tensor1<Real,Int>>  far_values;
        
        mutable bool metrics_initialized = false;
        
    public:
        
        virtual const std::map<KernelType, Tensor1<Real,Int>> & NearFieldValues() const override
        {
            return near_values;
        }

        virtual const std::map<KernelType, Tensor1<Real,Int>> & FarFieldValues() const override
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
            
            NearField();
            
            if( bct->IsSymmetric() && bct->NearUpperTriangular() )
            {
                Near_FillLowerTriangleFromUpperTriangle();
            }
            
            FarField ();
            
            if( bct->IsSymmetric() && bct->FarUpperTriangular() )
            {
                Far_FillLowerTriangleFromUpperTriangle();
            }
            
            ComputeDiagonals();
            
            metrics_initialized = true;
                
            ptoc(ClassName()+"::RequireMetrics");
            
        } // RequireMetrics
        
    protected:
        
        virtual void NearField() const
        {
            // Loops through all entries of the sparse matrix bct->Near() and computes the interaction data for the blocks.
            // Delegates the computation of the interaction data to the the near field kernel  N_ker (and its clones).
            // If bct->IsSymmetric() && bct->NearUpperTriangular() evaluates to 'true', then
            //      1. only the upper triangular part is traversed upon the computation pass, and
            //      2. the interaction data gets symmetrized afterwards via Near_FillLowerTriangleFromUpperTriangle, which delegates the transposition of blocks to N_ker (because only that knows what really has to happen).
            
            ptic(ClassName()+"::NearField");
            
            auto & near = bct->Near();
            
            N_ker->AllocateValueBuffers(near_values, near.NonzeroCount());
            
            if( near.NonzeroCount() <= 0 )
            {
                wprint(ClassName()+"::NearField no nonseparated blocks detected. Skipping.");
                ptoc(ClassName()+"::NearField");
                return;
            }
            
            const bool uppertriangular = bct->IsSymmetric() && bct->NearUpperTriangular();
            
            const auto & job_ptr = uppertriangular
                ? near.UpperTriangularJobPtr()
                : near.JobPtr();
            
            Int const * restrict const row_begin = uppertriangular
                ? near.Diag().data()
                : near.Outer().data();

            Int const * restrict const row_end = near.Outer().data()+1;
            Int const * restrict const inner   = near.Inner().data();
            
            const Int thread_count = job_ptr.Size()-1;
            
            #pragma omp parallel for num_threads( thread_count )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                std::unique_ptr<N_Kernel_T> N = N_ker->Clone();
                
                N->LoadValueBuffers(near_values);

                // This looks massive, but is only exchange of pointers.
                (void)N->LoadNearField(
                    bct->GetS().PrimitiveNearFieldData(),
                    bct->GetT().PrimitiveNearFieldData()
                );
                (void)N->LoadPrimitiveSerializedData(
                    bct->GetS().PrimitiveSerializedData(),
                    bct->GetT().PrimitiveSerializedData()
                );
                
                const Int i_begin = job_ptr[thread  ];
                const Int i_end   = job_ptr[thread+1];
                
                for( Int i = i_begin; i < i_end; ++i )
                {
//                        N->StartRow();
                    N->LoadS(i);

                    const Int k_begin = row_begin[i];
                    const Int k_end   = row_end  [i];
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        const Int j = inner[k];
                        
                        N->LoadT(j);
                    
                        N->Metric(k);

                    } // for (Int k = k_begin; k < k_end; ++k)
                    
//                        N->FinishRow( diag[i] );
                    
                } // for( Int i = i_begin; i < i_end; ++i )
                
            } // #pragma omp parallel
            
            ptoc(ClassName()+"::NearField");
        }
        
        virtual void FarField() const
        {
            // Loops through all entries of the sparse matrix bct->Far() and computes the interaction data for the blocks.
            // Delegates the computation of the interaction data to the the far field kernel  F_ker (and its clones).
            // If bct->IsSymmetric() && bct->FarUpperTriangular() evaluates to 'true', then
            //      1. only the upper triangular part is traversed upon the computation pass, and
            //      2. the interaction data gets symmetrized afterwards via Far_FillLowerTriangleFromUpperTriangle, which delegates the transposition of blocks to F_ker (because only that knows what really has to happen).
            
            ptic(ClassName()+"::FarField");

            auto & far = bct->Far();

            F_ker->AllocateValueBuffers(far_values, far.NonzeroCount());
            
            if( far.NonzeroCount() <= 0 )
            {
//                wprint(ClassName()+"::NearField no separated blocks detected. Skipping.");
                ptoc(ClassName()+"::FarField");
                return;
            }
            
            const bool uppertriangular = bct->IsSymmetric() && bct->FarUpperTriangular();
            
            {
                const auto & job_ptr = uppertriangular
                    ? far.UpperTriangularJobPtr()
                    : far.JobPtr();
                
                
                Int const * restrict const row_begin = uppertriangular
                    ? far.Diag().data()
                    : far.Outer().data();
                
                Int const * restrict const row_end   = far.Outer().data()+1;
                Int const * restrict const inner     = far.Inner().data();
                
                const Int thread_count = job_ptr.Size()-1;
                
                #pragma omp parallel for num_threads( thread_count )
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    std::unique_ptr<F_Kernel_T> F = F_ker->Clone();
                    
                    F->LoadValueBuffers(far_values);
                    
                    (void)F->LoadFarField(
                        bct->GetS().ClusterFarFieldData(),
                        bct->GetT().ClusterFarFieldData()
                    );

                    const Int i_begin = job_ptr[thread];
                    const Int i_end   = job_ptr[thread+1];
                    
                    for( Int i = i_begin; i < i_end; ++i )
                    {
                        F->LoadS(i);
    
                        const Int k_begin = row_begin[i];
                        const Int k_end   = row_end  [i];
                        
                        for( Int k = k_begin; k < k_end; ++k )
                        {
                            const Int j = inner[k];

                            F->LoadT(j);
                        
                            F->Metric(k);

                        } // for( Int k = k_begin; k < k_end; ++k )

                    } // for( Int i = i_begin; i < i_end; ++i )
                    
                } // #pragma omp parallel
            }
            
            ptoc(ClassName()+"::FarField");
        }
        

//##############################################################################################
//      Compute diagonals
//##############################################################################################
        
    protected:
        
        virtual void ComputeDiagonals() const
        {
            ptic(ClassName()+"::ComputeDiagonals");

            Int cols = 1;
            
            auto & T = bct->GetT();
            auto & S = bct->GetS();
            
            T.RequireBuffers(cols);
            S.RequireBuffers(cols);
            
//            S.CleanseBuffers();
//            T.CleanseBuffers();
            
//            S.PrimitiveOutputBuffer().Fill( static_cast<Real>(0) );
            
            {
                      Real * restrict const in_buffer = T.PrimitiveInputBuffer().data();
                const Real * restrict const near_data = T.PrimitiveNearFieldData().data();
                const Int near_dim = T.NearDim();
                const Int primitive_count = T.PrimitiveCount();
                
                for( Int i = 0; i < primitive_count; ++i )
                {
                    in_buffer[i] = near_data[ near_dim * i ];
                }
            }
            
            T.PrimitivesToClusters(false);
            
            T.PercolateUp();
            
            for( auto const & f : near_values)
            {
                ApplyNearFieldKernel( static_cast<Real>(1), f.first );
                
                ApplyFarFieldKernel ( static_cast<Real>(1), f.first );
                
                S.PercolateDown();
                
                S.ClustersToPrimitives( true );
            
                Real const * restrict const u = S.PrimitiveOutputBuffer().data();
                Real const * restrict const a = S.PrimitiveNearFieldData().data();
                Real       * restrict const values = near_values.find(f.first)->second.data();
                Int  const * restrict const diag_ptr = bct->Near().Diag().data();
                
                const Int m        = S.PrimitiveCount();
                const Int near_dim = S.NearDim();
                
                #pragma omp parallel for
                for( Int i = 0; i < m; ++i )
                {
                    values[diag_ptr[i]] -= u[i]/a[ near_dim * i ];
                }
            }
            
            ptoc(ClassName()+"::ComputeDiagonals");
            
        } // ComputeDiagonals
        
        
        
        
    public:
        
//##############################################################################################
//      Matrix multiplication
//##############################################################################################

        virtual void ApplyNearFieldKernel( const Real factor, const KernelType type ) const override
        {
            ptic(ClassName()+"::ApplyNearFieldKernel ("+KernelTypeName[type]+")" );
            
            bct->Near().Dot(
                near_values.find(type)->second.data(),
                factor,               bct->GetT().PrimitiveInputBuffer().data(),
                static_cast<Real>(0), bct->GetS().PrimitiveOutputBuffer().data(),
                bct->GetT().BufferDimension()
            );
            
            ptoc(ClassName()+"::ApplyNearFieldKernel ("+KernelTypeName[type]+")" );
        }
        
        virtual void ApplyFarFieldKernel( const Real factor, const KernelType type ) const override
        {
            ptic(ClassName()+"::ApplyFarFieldKernel ("+KernelTypeName[type]+")");
            
            bct->Far().Dot(
                far_values.find(type)->second.data(),
                factor,               bct->GetT().ClusterInputBuffer().data(),
                static_cast<Real>(0), bct->GetS().ClusterOutputBuffer().data(),
                bct->GetT().BufferDimension()
            );
            
            ptoc(ClassName()+"::ApplyFarFieldKernel ("+KernelTypeName[type]+")");
        }
        
        void Multiply_DenseMatrix(
            const Real    alpha, const ExtReal * X,
            const ExtReal  beta,       ExtReal * Y,
            const Int cols,
            const KernelType type
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
            
            T.Pre( X, cols, type );
            
            const Int buffer_dim = T.BufferDimension();

            S.RequireBuffers(buffer_dim); // Tell the S-side what it has to expect.
            
            ApplyNearFieldKernel( weight * alpha, type );
            
            ApplyFarFieldKernel ( weight * alpha, type );

            S.Post( Y, type, beta );
            
            ptoc(ClassName()+"::Multiply_DenseMatrix");
        }

        void Multiply_DenseMatrix(
            const Real    alpha, const Tensor1<ExtReal,Int> & X,
            const ExtReal beta,        Tensor1<ExtReal,Int> & Y, KernelType type
        ) const override
        {
            Multiply_DenseMatrix( alpha, X.data(), beta, Y.data(), X.Dimension(1), type );
        }
        
        void Multiply_DenseMatrix(
            const Real    alpha, const Tensor2<ExtReal,Int> & X,
            const ExtReal beta,        Tensor2<ExtReal,Int> & Y, KernelType type
        ) const override
        {
            if( X.Dimension(1) == Y.Dimension(1))
            {
                Multiply_DenseMatrix( alpha, X.data(), beta, Y.data(), X.Dimension(1), type );
            }
            else
            {
                eprint(ClassName()+"::Multiply_DenseMatrix : number of columns of input and output do not match.");
            }
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
        
    
        void Near_FillLowerTriangleFromUpperTriangle() const
        {
            ptic(ClassName()+"::Near_FillLowerTriangleFromUpperTriangle");
            
            // NearField
            if( bct->Near().WellFormed() )
            {
                const Int * restrict const diag   = bct->Near().Diag().data();
                const Int * restrict const outer  = bct->Near().Outer().data();
                const Int * restrict const inner  = bct->Near().Inner().data();
                
                const auto & job_ptr = bct->Near().LowerTriangularJobPtr();
                
                const Int thread_count = job_ptr.Size()-1;
                
                #pragma omp parallel for num_threads( thread_count )
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    auto N = N_ker->Clone();
                    
                    N->LoadValueBuffers( near_values );
                    
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
                            
                            N->TransposeBlock(L, k);
                            
                        } // for( Int k = k_begin; k < k_end; ++k )

                    } // for( Int i = i_begin; i < i_end; ++i )
                    
                } // #pragma omp parallel
            }

            ptoc(ClassName()+"::Near_FillLowerTriangleFromUpperTriangle");
        }
        
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
