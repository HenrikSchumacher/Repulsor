private:

    void Serialize( Cluster_T * const root )
    {
        ptic(className()+"::Serialize (OpenMP)");
        
        //            tree_max_depth = root->max_depth;
        //
        // We have to allocate these two arrays first, so that ClusterCount() and LeafClusterCount() return correct results.
        C_serialized  = Tensor2<SReal,Int>( root->descendant_count, C_proto[0]->Size() );
        leaf_clusters = Tensor1<Int,Int>( root->descendant_leaf_count );
        
        std::vector<std::thread> threads;
        threads.reserve(7);
        
        threads.emplace_back( [=](){ C_left  = Tensor1<Int,Int>( ClusterCount() ); } );
        threads.emplace_back( [=](){ C_right = Tensor1<Int,Int>( ClusterCount() ); } );
        threads.emplace_back( [=](){ C_begin = Tensor1<Int,Int>( ClusterCount() ); } );
        threads.emplace_back( [=](){ C_end   = Tensor1<Int,Int>( ClusterCount() ); } );
        threads.emplace_back( [=](){ C_depth = Tensor1<Int,Int>( ClusterCount() ); } );
        threads.emplace_back( [=](){ C_next  = Tensor1<Int,Int>( ClusterCount() ); } );
        threads.emplace_back( [=](){ leaf_cluster_lookup = Tensor1<Int,Int>( ClusterCount(), -1 ); } );
        
        for( auto & thread : threads )
        {
            thread.join();
        }
        
        #pragma omp parallel num_threads( ThreadCount() )
        {
            #pragma omp single nowait
            {
                serialize( root, ThreadCount() );
            }
        }
        
        depth = root->max_depth;
        
        {
            const Int last = PrimitiveCount();
            cptr<Int> ord     = P_ordering.data();
            mptr<Int> inv_ord = P_inverse_ordering.data();
            
            for( Int i = 0; i < last; ++i )
            {
                inv_ord[ord[i]] = i;
            }
        }
        
        leaf_cluster_ptr = Tensor1<Int,Int> ( LeafClusterCount() + 1 );
        leaf_cluster_ptr[0] = 0;
        
        for( Int i = 0; i < LeafClusterCount(); ++i )
        {
            leaf_cluster_ptr[ i + 1 ] = C_end[leaf_clusters[i]];
        }
        
        
        ptoc(className()+"::Serialize");
        
    } // Serialize
    
    void serialize( Cluster_T * C, const Int free_thread_count )
    {
        const Int thread = omp_get_thread_num();
        
        const Int g_ID = C->g_ID;
        
        // enumeration in depth-first order
        C_begin[g_ID] = C->begin;
        C_end  [g_ID] = C->end;
        C_depth[g_ID] = C->depth;
        C_next [g_ID] = g_ID + C->descendant_count;
        
        // TODO: Potentially, some false sharing could occur during write here. But that should be seldom as we use depth-first order.
        C_proto[thread]->SetPointer( C_thread_serialized.data(C->thread), C->l_ID );
        C_proto[thread]->Write( C_serialized.data(), g_ID );
        
        Cluster_T * L = C->left;
        Cluster_T * R = C->right;
        
        if( ( L != nullptr ) && ( R != nullptr ) )
        {
            L->g_ID  = C_left [g_ID] = g_ID + 1;
            L->leafs_before_count = C->leafs_before_count;
            
            R->g_ID = C_right[g_ID] = g_ID + 1 + L->descendant_count;
            R->leafs_before_count = C->leafs_before_count + L->descendant_leaf_count;
            //
            #pragma omp task final( free_thread_count < 1 )
            {
                serialize( L, free_thread_count/2 );
            }
            #pragma omp task final( free_thread_count < 1 )
            {
                serialize( R, free_thread_count - free_thread_count/2 );
            }
            #pragma omp taskwait
            
            // Cleaning up after ourselves to prevent a destructor cascade.
            delete L;
            delete R;
        }
        else
        {
            C_left [g_ID] = -1;
            C_right[g_ID] = -1;
            
            leaf_clusters[C->leafs_before_count] = g_ID;
            leaf_cluster_lookup[g_ID] = C->leafs_before_count;
        }
        
        ptoc(className()+"::Serialize (OpenMP)");
    } //serialize
