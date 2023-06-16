private:

    void Split( Cluster_T * root )
    {
        ptic(className()+"::Split (OpenMP)");
        
        P_score_buffer = iota<SReal,Int>( PrimitiveCount() );
        
        P_perm_buffer = Tensor1<Int,Int>( PrimitiveCount() );
        
        #pragma omp parallel num_threads( ThreadCount() ) shared( root )
        {
            #pragma omp single nowait
            {
                Split_Cluster( root, ThreadCount() );
            }
        }
        
        P_perm_buffer  = Tensor1<Int,Int>();
        P_score_buffer = Tensor1<SReal,Int>();
        
        ptoc(className()+"::Split (OpenMP)");
    } // Split

    void Split_Cluster( Cluster_T * C, const Int free_thread_count )
    {
        const Int thread = omp_get_thread_num();
        
        const Int begin = C->begin;
        const Int end   = C->end;
        
        if( end - begin > SplitThreshold() )
        {
            // TODO: Many things to do here:
            // Split finds a nice split of the cluster and reorders the primitives begin,...,end-1 so that
            // primtives begin,...,split_index-1 belong to left  new cluster
            // primtives split_index-1,...,end-1 belong to right new cluster
            // Split has to return a number split_index <= begin if it is not successful and a value begin < split_index < end otherwise.
            // Split is also responsible for computing the bounding volumes of the children, if successful.
            // Remark: Some bounding volume types, e.g., AABBs can use some information from the Split pass to compute the children's bounding volumes. This is why we merge the splitting pass with the computation of the children's bounding columes.
            
            // Remark: Make sure that bounding volumes are already computed for the child clusters. Moreover, we want that the serialized data is stored in the thread's storage that _created_ the new clusters. This is why we do NOT compute the bounding volumes at the beginning of Split; C is possibly created by another thread and we _must not_ write to that thread's memory.
            
            const Int  left_ID = thread_cluster_counter(thread,0)+1;
            const Int right_ID = thread_cluster_counter(thread,0)+2;

            const Int begin = C->begin;
            const Int end   = C->end;
            
            Primitive_T & P = *P_proto[thread];
            
            const Int split_index = C_proto[thread]->Split(
                 P,                                                  // prototype for primitves
                 P_serialized.data(), begin, end,                    // which primitives are in question
                 P_ordering.data(),                                  // which primitives are in question
                 C_thread_serialized.data(C->thread), C->l_ID,       // where to get   the bounding volume info for current cluster
                 C_thread_serialized.data(   thread),  left_ID,      // where to store the bounding volume info for left  child (if successful!)
                 C_thread_serialized.data(   thread), right_ID,      // where to store the bounding volume info for right child (if successful!)
                 P_score_buffer.data(),                              // some scratch space for storing local scores
                 P_perm_buffer.data(),                               // scratch space for storing local permutations
                 P_inverse_ordering.data(),                          // abusing P_inverse_ordering as scratch space for storing inverses of local permutations
                 thread_count
            );
            
            // Remark: We implicitly assume that the split was successful.
            
            // create new nodes...
            thread_cluster_counter(thread,0) += 2;
            
            // We use raw pointers for performance reasons because we want to delete Cluster instances in the parallel serialization pass.
            C->left  = new Cluster_T ( thread,  left_ID, begin,       split_index, C->depth+1 );
            C->right = new Cluster_T ( thread, right_ID, split_index, end,         C->depth+1 );
            
            // ... and split them in parallel
            #pragma omp task final(free_thread_count<1)
            {
                Split_Cluster( C->left, free_thread_count/2 );
            }
            #pragma omp task final(free_thread_count<1)
            {
                Split_Cluster( C->right, free_thread_count - free_thread_count/2 );
            }
            #pragma omp taskwait
            
            // collecting statistics for the later serialization
            // counting ourselves as descendant, too!
            C->descendant_count = 1 + C->left->descendant_count + C->right->descendant_count;
            C->descendant_leaf_count = C->left->descendant_leaf_count + C->right->descendant_leaf_count;
            C->max_depth = std::max( C->left->max_depth, C->right->max_depth );
        }
        else
        {
            // count cluster as leaf cluster
            // counting ourselves as descendant, too!
            C->descendant_count = 1;
            C->descendant_leaf_count = 1;
            return;
        }
    } // Split_Cluster
