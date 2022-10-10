protected:

    virtual void FarField() const
    {
        // Loops through all entries of the sparse matrix bct->Far() and computes the interaction data for the blocks.
        // Delegates the computation of the interaction data to the the far field kernel  F_ker (and its clones).
        // If bct->IsSymmetric() evaluates to 'true', then
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
        
        const bool uppertriangular = bct->IsSymmetric();
        
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
                    
                        F->ComputeBlock(k);

                    } // for( Int k = k_begin; k < k_end; ++k )

                } // for( Int i = i_begin; i < i_end; ++i )
                
            } // #pragma omp parallel
        }
        
        
        if( bct->IsSymmetric() )
        {
            Far_FillLowerTriangleFromUpperTriangle();
        }
        
        ptoc(ClassName()+"::FarField");
    }
