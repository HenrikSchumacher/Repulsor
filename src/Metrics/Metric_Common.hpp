//##############################################################################################
//      Computation kernels
//##############################################################################################

public:

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
        
        N_ker.AllocateValueBuffers(near_values, near.NonzeroCount());
        
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
            N_Kernel_T N (N_ker);
            
            N.LoadValueBuffers(near_values);

            // This looks massive, but is only exchange of pointers.
            (void)N.LoadNearField(
                bct->GetS().PrimitiveNearFieldData(),
                bct->GetT().PrimitiveNearFieldData()
            );
            (void)N.LoadPrimitiveSerializedData(
                bct->GetS().PrimitiveSerializedData(),
                bct->GetT().PrimitiveSerializedData()
            );
            
            const Int i_begin = job_ptr[thread  ];
            const Int i_end   = job_ptr[thread+1];
            
            for( Int i = i_begin; i < i_end; ++i )
            {
                N.LoadS(i);

                const Int k_begin = row_begin[i];
                const Int k_end   = row_end  [i];
                
                for( Int k = k_begin; k < k_end; ++k )
                {
                    const Int j = inner[k];
                    
                    N.LoadT(j);
                
                    N.Metric(k);

                } // for (Int k = k_begin; k < k_end; ++k)
                
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
                F_Kernel_T F (F_ker);
                
                F.LoadValueBuffers(far_values);
                
                (void)F.LoadFarField(
                    bct->GetS().ClusterFarFieldData(),
                    bct->GetT().ClusterFarFieldData()
                );

                const Int i_begin = job_ptr[thread];
                const Int i_end   = job_ptr[thread+1];
                
                for( Int i = i_begin; i < i_end; ++i )
                {
                    F.LoadS(i);

                    const Int k_begin = row_begin[i];
                    const Int k_end   = row_end  [i];
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        const Int j = inner[k];

                        F.LoadT(j);
                    
                        F.Metric(k);

                    } // for( Int k = k_begin; k < k_end; ++k )

                } // for( Int i = i_begin; i < i_end; ++i )
                
            } // #pragma omp parallel
        }
        
        ptoc(ClassName()+"::FarField");
    }


public:
    
//##############################################################################################
//      Matrix multiplication
//##############################################################################################

    virtual void ApplyNearFieldKernel( const Real factor, const KernelType type ) const override
    {
        ptic(ClassName()+"::ApplyNearFieldKernel ("+KernelTypeName[type]+")" );
        
        eprint(ClassName()+"::ApplyNearFieldKernel is yet to be implemented.");
    //            bct->Near().Multiply_DenseMatrix(
    //                factor,
    //                near_values.find(type)->second.data(),
    //                bct->GetT().PrimitiveInputBuffer().data(),
    //                static_cast<Real>(0),
    //                bct->GetS().PrimitiveOutputBuffer().data(),
    //                bct->GetT().BufferDimension()
    //            );
        
        ptoc(ClassName()+"::ApplyNearFieldKernel ("+KernelTypeName[type]+")" );
    }

    virtual void ApplyFarFieldKernel( const Real factor, const KernelType type ) const override
    {
        ptic(ClassName()+"::ApplyFarFieldKernel ("+KernelTypeName[type]+")");
        
        eprint(ClassName()+"::ApplyFarFieldKernel is yet to be implemented.");
        
    //            bct->Far().Multiply_DenseMatrix(
    //                factor,
    //                far_values.find(type)->second.data(),
    //                bct->GetT().ClusterInputBuffer().data(),
    //                static_cast<Real>(0),
    //                bct->GetS().ClusterOutputBuffer().data(),
    //                bct->GetT().BufferDimension()
    //            );
        
        ptoc(ClassName()+"::ApplyFarFieldKernel ("+KernelTypeName[type]+")");
    }
