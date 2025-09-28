private:

    void Split( Cluster_T * root )
    {
        TOOLS_PTIMER(timer,ClassName()+"::Split");
        
        P_score_buffer = iota<SReal,Int>( PrimitiveCount() );
        
        P_perm_buffer = Tensor1<Int,Int>( PrimitiveCount() );
        
        const Int top_level_count = int_cast<Int>(settings.parallel_perc_depth);
        
        logprint("Breadth-first scan to split the first few levels.");
        
        tree_rows_ptr = std::vector<std::vector<Cluster_T *>>(top_level_count+1);
        
        tree_rows_ptr[0].push_back(root);

        // TODO: The distribution of threads is still not optimal.
        
        std::mutex tree_rows_ptr_mutex;
        
        Split_Cluster<false>( 0, tree_rows_ptr[0][0], tree_rows_ptr_mutex, ThreadCount() );
        
        for( Int level = 1; level < top_level_count; ++level )
        {
            mref<std::vector<Cluster_T *>> row = tree_rows_ptr[level];

            const Int row_size  = static_cast<Int>(row.size());

            ParallelDo_Dynamic(
                [=,this,&row,&tree_rows_ptr_mutex]( const Int thread, const Int i )
                {
                    Split_Cluster<false>(
                        thread,
                        row[i],
                        tree_rows_ptr_mutex,
                        Max( one, ThreadCount() / row_size )
                    );
                },
                null, row_size, one,
                Min(row_size, ThreadCount() )
           );
        }
        
        // Now `tree_rows_ptr` contains all nodes of level at most `top_level_count`.
        // In particular `tree_rows_ptr[top_level_count]` contains all the nodes at level `top_level_count`.
        
        
        logprint("Parallel depth-first scans to create the subtrees.");
        
        mref<std::vector<Cluster_T *>> tree_row = tree_rows_ptr[top_level_count];
        
        ParallelDo_Dynamic(
            [=,this,&tree_row,&tree_rows_ptr_mutex]( const Int thread, const Int i )
            {
                Split_Cluster<true>( thread, tree_row[i], tree_rows_ptr_mutex, one );
            },
            null, static_cast<Int>(tree_row.size()), one,
            ThreadCount()
        );
        
        logprint("Reverse breadth-first scan to calculate statistics.");
        for( Int level = top_level_count; level --> 0 ; )
        {
            for( Cluster_T * C : tree_rows_ptr[level] )
            {
                if( C->left != nullptr  )
                {
                    Split_Cluster_Post(C);
                }
            }
        }
        
        P_perm_buffer  = Tensor1<Int,Int>();
        P_score_buffer = Tensor1<SReal,Int>();
    }


    template<bool recursive>
    void Split_Cluster( const Int thread, Cluster_T * C, std::mutex & tree_rows_ptr_mutex, const Int thread_count )
    {
        
        const Int begin = C->begin;
        const Int end   = C->end;
        
        if( end - begin > SplitThreshold() )
        {
            // Split finds a nice split of the cluster and reorders the primitives begin,...,end-1 so that
            // primtives begin,...,split_index-1 belong to left  new cluster
            // primtives split_index-1,...,end-1 belong to right new cluster
            // Split has to return a number split_index <= begin if it is not successful and a value begin < split_index < end otherwise.
            // Split is also responsible for computing the bounding volumes of the children, if successful.
            // Remark: Some bounding volume types, e.g., AABBs can use some information from the Split pass to compute the children's bounding volumes. This is why we merge the splitting pass with the computation of the children's bounding columes.
            
            // Remark: Make sure that bounding volumes are already computed for the child clusters. Moreover, we want that the serialized data is stored in the thread's storage that _created_ the new clusters. This is why we do NOT compute the bounding volumes at the beginning of Split; C is possibly created by another thread and we _must not_ write to that thread's memory.
            
            const Int  left_ID = thread_cluster_counter(thread,0)+1;
            const Int right_ID = thread_cluster_counter(thread,0)+2;
            
            mref<Primitive_T> P = *P_proto[thread];
            
            const Int split_index = C_proto[thread]->Split(
                P,                                                  // prototype for primitives
                P_serialized.data(), begin, end,                    // which primitives are in question
                P_ordering.data(),                                  // which primitives are in question
                C_thread_serialized.data(C->thread), C->l_ID,       // where to get   the bounding volume info for current cluster
                C_thread_serialized.data(   thread),  left_ID,      // where to store the bounding volume info for left  child (if successful!)
                C_thread_serialized.data(   thread), right_ID,      // where to store the bounding volume info for right child (if successful!)
                P_score_buffer.data(),                              // some scratch space for storing local scores
                P_perm_buffer.data(),                               // scratch space for storing local permutations
                P_inverse_ordering.data(),                          // abusing P_inverse_ordering as scratch space for storing inverses of local permutations
                //                 ( recursive ? static_cast<Int>(1) : ThreadCount() )
                thread_count
            );
            
            // Remark: We implicitly assume that the split was successful.
            
            // create new nodes...
            thread_cluster_counter(thread,0) += 2;
            
            const Int next_depth = C->depth+1;
            
            // We use raw pointers for performance reasons because we want to delete Cluster instances in the parallel serialization pass.
            C->left  = new Cluster_T ( thread,  left_ID, begin,       split_index, next_depth );
            C->right = new Cluster_T ( thread, right_ID, split_index, end,         next_depth );
            
            if constexpr( recursive )
            {
                Split_Cluster<recursive>( thread, C->left,  tree_rows_ptr_mutex, one );
                Split_Cluster<recursive>( thread, C->right, tree_rows_ptr_mutex, one );
                
                Split_Cluster_Post(C);
            }
            else
            {
                const std::lock_guard<std::mutex> lock( tree_rows_ptr_mutex );
                
                tree_rows_ptr[next_depth].push_back( C->left  );
                tree_rows_ptr[next_depth].push_back( C->right );
            }
        }
        else
        {
            // count cluster as leaf cluster
            // counting ourselves as descendant, too!
            C->descendant_count      = one;
            C->descendant_leaf_count = one;
        }

    }

    void Split_Cluster_Post( mptr<Cluster_T> C )
    {
        cptr<Cluster_T> L = C->left;
        cptr<Cluster_T> R = C->right;
        
        // collecting statistics for the later serialization
        // counting ourselves as descendant, too!
        C->descendant_count = one + L->descendant_count + R->descendant_count;
        C->descendant_leaf_count = L->descendant_leaf_count + R->descendant_leaf_count;
        C->max_depth = Max( L->max_depth, R->max_depth );
    }
