public:

    virtual cref<Tensor1<Int,Int>> NestedDissectionOrdering( const Int local_thread_count = 1 ) const override
    {
        std::string tag ("NestedDissectionOrdering");
        
        if( !this->InPersistentCacheQ(tag))
        {
            ptic(ClassName()+"::"+tag);
            
            auto & tree = GetClusterTree();
            
            // Single-threaded is actually faster here.
            // TODO: Maybe this is just because of the mutex?
            
            const Int m = VertexCount();
            const Int n = SimplexCount();
            
            Tensor1<Int,Int> perm_0 = iota<Int>( m );
            Tensor1<Int,Int> perm_1 ( m );
            
            cref<SparseMatrix_T> Adj = tree.LowOrderPostProcessor();
            
            cref<Tensor1<Int,Int>> C_begin = tree.ClusterBegin();
            cref<Tensor1<Int,Int>> C_end   = tree.ClusterEnd();
            cref<Tensor1<Int,Int>> C_left  = tree.ClusterLeft();
            cref<Tensor1<Int,Int>> C_right = tree.ClusterRight();
            
            Tensor1<Int,Int> V_begin ( tree.ClusterCount() );
            Tensor1<Int,Int> V_end   ( tree.ClusterCount() );
            
            Tensor2<Real,Int> v ( n, 2 );
            Tensor2<Real,Int> w ( m, 2 );
            
            constexpr Int16 s_zero = 0;
            constexpr Int16 s_one  = 1;
            constexpr Int16 s_two  = 2;
            
            Tensor1<Int16,Int> type ( m );
            
            V_begin[0] = 0;
            V_end  [0] = m;
            
            std::vector<Int> queue_0;
            std::vector<Int> queue_1;
            
            std::mutex mutex;
            
            queue_0.reserve( tree.ClusterCount() );
            queue_1.reserve( tree.ClusterCount() );
            
            // Push root cluster.
            queue_0.push_back( Scalar::Zero<Int> );
            
            Int level = 0;
            
            while( queue_0.size() > 0 )
            {
                debug_print("Starting level " + ToString(level) + "." );
                
                queue_1.resize(0);
                
                debug_print("queue_0 = "+ToString(queue_0));
                
                perm_1.Read( perm_0.data() );
                
                v.SetZero();
                
                debug_print("Compute indicators.");
                
//                ptic("Compute indicators");
                ParallelDo(
                    [&]( const Int k )
                    {
                        const Int C = queue_0[k];
                        const Int L = C_left [C];
                        const Int R = C_right[C];
                        
                        if( L >= 0 )
                        {
                            
                            const Int L_begin = C_begin[L];
                            const Int L_end   = C_end  [L];
                            
                            const Int R_begin = C_begin[R];
                            const Int R_end   = C_end  [R];
                            
                            for( Int i = L_begin; i < L_end; ++i )
                            {
                                v[i][0] = Scalar::One<Real>;
                            }
                            
                            for( Int i = R_begin; i < R_end; ++i )
                            {
                                v[i][1] = Scalar::One<Real>;
                            }
                        }
                    },
                    static_cast<Int>(queue_0.size()), local_thread_count
                );
//                ptoc("Compute indicators");
                
                Adj.template Dot_<2>(
                    Scalar::One <Real>, v.data(), 2,
                    Scalar::Zero<Real>, w.data(), 2,
                    2
                );

                debug_print("Compute types.");
//                ptic("Compute types");
                ParallelDo(
                    [&]( const Int i )
                    {
                        const Int j = perm_0[i];
                        
                        type[i] = - s_one
                            +
                            ( (w[j][0] > Scalar::Zero<Real>) ? s_one : s_zero )
                            +
                            ( (w[j][1] > Scalar::Zero<Real>) ? s_two : s_zero );
                    },
                    m, ThreadCount()
                );
//                ptoc("Compute types");
                
                debug_print("Modify permutation.");
//                ptic("Modify permutation");
                ParallelDo(
                    [&]( const Int k )
                    {
                        const Int C = queue_0[k];
                        const Int L = C_left [C];
                        const Int R = C_right[C];

                        if( L >= 0 )
                        {
                            const Int i_begin = V_begin[C];
                            const Int i_end   = V_end  [C];
                            
                            Tiny::Vector<3,Int,Int> ctr ( Scalar::Zero<Int> );
                            
                            // Compute counters;
                            for( Int i = i_begin; i < i_end; ++i )
                            {
                                const Int16 t = type[i];
                                
                                if( t >= s_zero )
                                {
                                    ++ctr[t];
                                }
                            }
                            
                            Tiny::Vector<4,Int,Int> pos;
                            
                            pos[0] = i_begin;
                            pos[1] = pos[0] + ctr[0];
                            pos[2] = pos[1] + ctr[1];
                            pos[3] = pos[2] + ctr[2];
                            
                            if( (ctr[0] > 0) && (ctr[1] > 0) && ( pos[3] == i_end ) )
                            {
                                V_begin(L) = pos[0];
                                V_end  (L) = pos[1];
                                V_begin(R) = pos[1];
                                V_end  (R) = pos[2];
                                
                                {
                                    const std::lock_guard<std::mutex> lock( mutex );
                                    queue_1.push_back( L );
                                    queue_1.push_back( R );
                                }
                                
                                // Modify permutation;
                                for( Int i = i_begin; i < i_end; ++i )
                                {
                                    const Int16 t = type[i];
                                    
                                    if( t >= s_zero )
                                    {
                                        const Int p = pos[t]++;
                                        
                                        perm_1[p] = perm_0[i];
                                    }
                                }
                            }
                        }
                    },
                    static_cast<Int>(queue_0.size()), local_thread_count
                );
//                ptoc("Modify permutation");
                
                debug_print("Swapping.");
                
                std::swap( queue_0, queue_1 );
                
                swap( perm_0, perm_1 );
                
                debug_print("Done with level " + ToString(level) + "." );
                
                ++level;
            }
            
            this->SetPersistentCache( tag, std::any( std::move(perm_0) ) );
            
            ptoc(ClassName()+"::"+tag);
        }
        
        return std::any_cast<Tensor1<Int,Int> &>( this->GetPersistentCache(tag) );
    }
