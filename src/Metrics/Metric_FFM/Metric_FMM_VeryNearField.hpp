protected:

    virtual void VeryNearField() const
    {
        // Loops through all entries of the sparse matrix bct->VeryNear() and computes the interaction data for the blocks.
        // Delegates the computation of the interaction data to the the very near field kernel  A_ker (and its clones).
        // If bct->IsSymmetric() evaluates to 'true', then
        //      1. only the upper triangular part is traversed upon the computation pass, and
        //      2. the interaction data gets symmetrized afterwards via VeryNear_FillLowerTriangleFromUpperTriangle, which delegates the transposition of blocks to A_ker (because only that knows what really has to happen).
        
        ptic(ClassName()+"::VeryNearField");
        
        const auto & matrix = bct->VeryNear();
        
        A_ker->AllocateValueBuffers(
            verynear_values,
            matrix.NonzeroCount() * A_ker->ValueSize()
        );
        
        if( verynear.NonzeroCount() <= 0 )
        {
            wprint(ClassName()+"::VeryNearField no nonseparated blocks detected. Skipping.");
            ptoc(ClassName()+"::VeryNearField");
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
            std::unique_ptr<A_Kernel_T> A = A_ker->Clone();
            
            A->LoadValueBuffers(verynear_values);

            // This looks massive, but is only exchange of pointers.
            (void)A->LoadNearField(
                bct->GetS().PrimitiveNearFieldData(),
                bct->GetT().PrimitiveNearFieldData()
            );
            (void)A->LoadPrimitiveSerialized(
                bct->GetS().PrimitiveSerialized(),
                bct->GetT().PrimitiveSerialized()
            );
            
            const Int i_begin = job_ptr[thread  ];
            const Int i_end   = job_ptr[thread+1];
            
            for( Int i = i_begin; i < i_end; ++i )
            {
                A->LoadS(i);

                const Int k_begin = row_begin[i];
                const Int k_end   = row_end  [i];
                
                for( Int k = k_begin; k < k_end; ++k )
                {
                    const Int j = inner[k];
                    
                    A->LoadT(j);
                
                    A->ComputeBlock(k);

                } // for (Int k = k_begin; k < k_end; ++k)
                
            } // for( Int i = i_begin; i < i_end; ++i )
            
        } // #pragma omp parallel
        
        
        if( bct->IsSymmetric() )
        {
            VeryNear_FillLowerTriangleFromUpperTriangle();
        }
        
        ptoc(ClassName()+"::VeryNearField");
    }


    void VeryNear_FillLowerTriangleFromUpperTriangle() const
    {
        ptic(ClassName()+"::VeryNear_FillLowerTriangleFromUpperTriangle");
        
        const auto & matrix = bct->VeryNear();
        
        // NearField
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
                auto A = A_ker->Clone();
                
                A->LoadValueBuffers( near_values );
                
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
                        
                        A->TransposeBlock(L, k);
                        
                    } // for( Int k = k_begin; k < k_end; ++k )

                } // for( Int i = i_begin; i < i_end; ++i )
                
            } // #pragma omp parallel
        }

        ptoc(ClassName()+"::VeryNear_FillLowerTriangleFromUpperTriangle");
    }

    virtual void ApplyVeryNearFieldKernel( const Real alpha ) const override
    {
        ptic(ClassName()+"::ApplyVeryNearFieldKernel" );
        
        const auto & matrix = bct->VeryNear();
        
        const Int  * restrict const rp = matrix.Outer().data();
        const Int  * restrict const ci = matrix.Inner().data();
        const Real * restrict const a  = verynear_values.data();
        const Real * restrict const X  = bct->GetT().ClusterInputBuffer().data();
              Real * restrict const Y  = bct->GetT().ClusterOutputBuffer().data();

        const auto & job_ptr = matrix.JobPtr();
        
        const Int thread_count = job_ptr.Size()-1;
        
        #pragma omp parallel for num_threads( thread_count )
        for( Int thread = 0; thread < thread_count; ++thread )
        {
            std::unique_ptr<A_Kernel_T> A = A_ker->Clone();
            
            A->LoadInputBuffer(X);
            A->LoadOutputBuffer(Y);
            
            const Int i_begin = job_ptr[thread  ];
            const Int i_end   = job_ptr[thread+1];
            
            for( Int i = i_begin; i < i_end; ++i )
            {
                const Int l_begin = rp[i  ];
                const Int l_end   = rp[i+1];
                
                __builtin_prefetch( &ci[l_end] );
                
                __builtin_prefetch( &a[l_end] );
                
                A->ClearVector();

                for( Int l = l_begin; l < l_end-1; ++l )
                {
                    A->ApplyBlock( alpha, l, ci[l] );
                }
                
                // incorporate the local updates into Y-buffer
                A->WriteVector(i);
            }
        }
        
        ptoc(ClassName()+"::ApplyVeryNearFieldKernel" );
    }
