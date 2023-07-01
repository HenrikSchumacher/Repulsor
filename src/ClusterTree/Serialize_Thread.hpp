private:

    void Serialize( Cluster_T * const root )
    {
        ptic(className()+"::Serialize");
        
        //            tree_max_depth = root->max_depth;
        //
        // We have to allocated these two arrays first, so that ClusterCount() and LeafClusterCount() return correct results.
        C_serialized  = Tensor2<SReal,Int>( root->descendant_count, C_proto[0]->Size() );
        leaf_clusters = Tensor1<Int,Int>( root->descendant_leaf_count );
        
        std::vector<std::thread> threads;
        threads.reserve(7);
        threads.emplace_back( [=](){ C_left  = Tensor1<Int,Int>( ClusterCount() ); } );
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
        
        
        logprint("Breadth-first scan to serialize the first few levels.");
        
        const Int top_level_count = int_cast<Int>(settings.parallel_perc_depth);
        
        for( Int level = 0; level < top_level_count; ++level )
        {
            for( Cluster_T * C : tree_rows_ptr[level] )
            {
                Serialize_Cluster<false>( null, C );
            }
        }
        
        logprint("Parallel depth-first scans for each subtree.");
        
        std::vector<Cluster_T *> & restrict tree_row = tree_rows_ptr[top_level_count];
        
        ParallelDo_Dynamic(
            [=,&tree_row]( const Int thread, const Int i )
            {
                Serialize_Cluster<true>( thread, tree_row[i] );
            },
            null, static_cast<Int>(tree_row.size()), one,
            ThreadCount()
        );
        

        logprint("Reverse breadth-first scans for clean-up.");
        for( Int level = top_level_count; level --> 0 ; )
        {
            for( Cluster_T * C : tree_rows_ptr[level] )
            {
                if( C->left )
                {
                    Serialize_Cluster_Post( C );
                }
            }
        }
        
        depth = root->max_depth;
        
        {
            const Int last = PrimitiveCount();
            ptr<Int> ord     = P_ordering.data();
            mut<Int> inv_ord = P_inverse_ordering.data();
            
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
        
        
        tree_rows_ptr = std::vector<std::vector<Cluster_T *>>();
        
        ptoc(className()+"::Serialize");
        
    } // Serialize
    
    template<bool recursive>
    void Serialize_Cluster( const Int thread, Cluster_T * C )
    {
        const Int g_ID = C->g_ID;
        
        // enumeration in depth-first order
        C_begin[g_ID] = C->begin;
        C_end  [g_ID] = C->end;
        C_depth[g_ID] = C->depth;
        C_next [g_ID] = g_ID + C->descendant_count;
        
        C_proto[thread]->SetPointer( C_thread_serialized.data(C->thread), C->l_ID );
        C_proto[thread]->Write( C_serialized.data(), g_ID );
        
        Cluster_T * L = C->left;
        Cluster_T * R = C->right;
        
        if( L != nullptr )
        {
            L->g_ID  = C_left [g_ID] = g_ID + 1;
            L->leafs_before_count = C->leafs_before_count;
            
            R->g_ID = C_right[g_ID] = g_ID + 1 + L->descendant_count;
            R->leafs_before_count = C->leafs_before_count + L->descendant_leaf_count;
            
            if constexpr( recursive )
            {
                Serialize_Cluster<recursive>( thread, L );
                Serialize_Cluster<recursive>( thread, R );

                Serialize_Cluster_Post(C);
            }
        }
        else
        {
            C_left [g_ID] = -1;
            C_right[g_ID] = -1;
            
            leaf_clusters[C->leafs_before_count] = g_ID;
            leaf_cluster_lookup[g_ID] = C->leafs_before_count;
        }
    } // Serialize_Cluster


    void Serialize_Cluster_Post( Cluster_T * C )
    {
        
        // Cleaning up after ourselves to prevent a destructor cascade.
        delete C->left;
        delete C->right;
        
    } // Serialize_Cluster_Post
