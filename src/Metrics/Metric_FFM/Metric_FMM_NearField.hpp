protected:

    virtual void NearField() const
    {
        // Loops through all entries of the sparse matrix bct->Near() and computes the interaction data for the blocks.
        // Delegates the computation of the interaction data to the the near field kernel  N_ker (and its clones).
        // If bct->IsSymmetric() && bct->NearUpperTriangular() evaluates to 'true', then
        //      1. only the upper triangular part is traversed upon the computation pass, and
        //      2. the interaction data gets symmetrized afterwards via Near_FillLowerTriangleFromUpperTriangle, which delegates the transposition of blocks to N_ker (because only that knows what really has to happen).
        
        ptic(ClassName()+"::NearField");
        
        const auto & matrix = bct->Near();
        
        N_ker->AllocateValueBuffers(near_values, matrix.NonzeroCount());
        
        if( matrix.NonzeroCount() <= 0 )
        {
            wprint(ClassName()+"::NearField no nonseparated blocks detected. Skipping.");
            ptoc(ClassName()+"::NearField");
            return;
        }
        
        const bool uppertriangular = bct->IsSymmetric();
        
        const auto & job_ptr = uppertriangular
            ? matrix.UpperTriangularJobPtr()
            : matrix.JobPtr();
        
        Int const * restrict const row_begin = uppertriangular
            ? matrix.Diag().data()
            : matrix.Outer().data();

        Int const * restrict const row_end = matrix.Outer().data()+1;
        Int const * restrict const inner   = matrix.Inner().data();
        
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
                N->LoadS(i);

                const Int k_begin = row_begin[i];
                const Int k_end   = row_end  [i];
                
                for( Int k = k_begin; k < k_end; ++k )
                {
                    const Int j = inner[k];
                    
                    N->LoadT(j);
                
                    N->ComputeBlock(k);

                } // for (Int k = k_begin; k < k_end; ++k)
                
            } // for( Int i = i_begin; i < i_end; ++i )
            
        } // #pragma omp parallel
        
        
        if( bct->IsSymmetric() )
        {
            Near_FillLowerTriangleFromUpperTriangle();
        }
        
        ptoc(ClassName()+"::NearField");
    }


    void Near_FillLowerTriangleFromUpperTriangle() const
    {
        ptic(ClassName()+"::Near_FillLowerTriangleFromUpperTriangle");
        
        const auto & matrix = bct->Near();
        
        if( matrix.WellFormed() )
        {
            const Int * restrict const diag   = matrix.Diag().data();
            const Int * restrict const outer  = matrix.Outer().data();
            const Int * restrict const inner  = matrix.Inner().data();
            
            const auto & job_ptr = matrix.LowerTriangularJobPtr();
            
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


//##############################################################################################
//      Matrix multiplication
//##############################################################################################

        
    virtual void ApplyNearFieldKernel( const Real alpha ) const override
    {
        ptic(ClassName()+"::ApplyNearFieldKernel" );
        
        const auto & matrix = bct->Near();
        
        const Int  * restrict const rp = matrix.Outer().data();
        const Int  * restrict const ci = matrix.Inner().data();
        const Real * restrict const a  = near_values.data();
        const Real * restrict const X  = bct->GetT().ClusterInputBuffer().data();
              Real * restrict const Y  = bct->GetT().ClusterOutputBuffer().data();
        
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
