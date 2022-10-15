protected:
            
    void FindCollidingClusters_DFS( const Int i_0, const Int j_0 ) const
    {
        const Int thread = omp_get_thread_num();
        
        Int i_stack[max_depth] = {};
        Int j_stack[max_depth] = {};
        
        Int stack_ptr = null;
        i_stack[0] = i_0;
        j_stack[0] = j_0;
        
        std::vector<Int>   & collision_idx  = collision_i[thread];
        std::vector<Int>   & collision_jdx  = collision_j[thread];
        std::vector<SReal> & collision_time = collision_t[thread];
        
        const Int * restrict const S_left   = S.ClusterLeft().data();
        const Int * restrict const S_right  = S.ClusterRight().data();
        
        const Int * restrict const T_left   = T.ClusterLeft().data();
        const Int * restrict const T_right  = T.ClusterRight().data();
        
        const Tensor2<SReal,Int> & S_C_ser    = S.ClusterSerializedData();
        const Tensor2<SReal,Int> & S_C_up_ser = S.ClusterUpdatedSerializedData();
        const Tensor2<SReal,Int> & T_C_ser    = T.ClusterSerializedData();
        const Tensor2<SReal,Int> & T_C_up_ser = T.ClusterUpdatedSerializedData();
        //
        SReal a = static_cast<SReal>(0);
        SReal b = static_cast<SReal>(1);
        
        while( (null <= stack_ptr) && (stack_ptr < max_depth - 4) )
        {
            const Int i = i_stack[stack_ptr];
            const Int j = j_stack[stack_ptr];
            stack_ptr--;
            
            AABB_CollisionTimeInterval<AMB_DIM, SReal, Int> (
                S_C_ser.data(i), S_C_up_ser.data(i),
                T_C_ser.data(j), T_C_up_ser.data(j),
                a, b
            );
            
            if( a < static_cast<SReal>(1) )
            {
                const Int left_i = S_left[i];
                const Int left_j = T_left[j];
                
                // Warning: This assumes that both children in a cluster tree are either defined or empty.
                if( left_i >= null || left_j >= null )
                {
                    
                    const Int right_i = S_right[i];
                    const Int right_j = T_right[j];
                    
                    // TODO: Improve score.
                    
                    const SReal score_i = static_cast<SReal>(left_i >= null) * S_C_ser(i,0);
                    const SReal score_j = static_cast<SReal>(left_j >= null) * T_C_ser(j,0);
                    
                    if( (score_i == score_j) /*&& (score_i > static_cast<SReal>(0)) && score_j > static_cast<SReal>(0)*/ )
                    {
                        if constexpr ( is_symmetric )
                        {
                            if( i == j )
                            {
                                //  Creating 3 blockcluster children, since there is one block is just the mirror of another one.
                                
                                ++stack_ptr;
                                i_stack[stack_ptr] = left_i;
                                j_stack[stack_ptr] = right_j;
                                
                                ++stack_ptr;
                                i_stack[stack_ptr] = right_i;
                                j_stack[stack_ptr] = right_j;
                                
                                ++stack_ptr;
                                i_stack[stack_ptr] = left_i;
                                j_stack[stack_ptr] = left_j;
                            }
                            else
                            {
                                // tie breaker: split both clusters
                                
                                // This is a very seldom case; still required to preserve symmetry.
                                // This happens only if i and j represent_diffent clusters with same radii.
                                
                                ++stack_ptr;
                                i_stack[stack_ptr] = right_i;
                                j_stack[stack_ptr] = right_j;
                                
                                ++stack_ptr;
                                i_stack[stack_ptr] = left_i;
                                j_stack[stack_ptr] = right_j;
                                
                                ++stack_ptr;
                                i_stack[stack_ptr] = right_i;
                                j_stack[stack_ptr] = left_j;
                                
                                ++stack_ptr;
                                i_stack[stack_ptr] = left_i;
                                j_stack[stack_ptr] = left_j;
                                
                            }
                        }
                        else
                        {
                            ++stack_ptr;
                            i_stack[stack_ptr] = right_i;
                            j_stack[stack_ptr] = right_j;
                            
                            ++stack_ptr;
                            i_stack[stack_ptr] = left_i;
                            j_stack[stack_ptr] = right_j;
                            
                            ++stack_ptr;
                            i_stack[stack_ptr] = right_i;
                            j_stack[stack_ptr] = left_j;
                            
                            ++stack_ptr;
                            i_stack[stack_ptr] = left_i;
                            j_stack[stack_ptr] = left_j;
                        }
                    }
                    else
                    {
                        // split only larger cluster
                        if (score_i > score_j)
                        {
                            ++stack_ptr;
                            i_stack[stack_ptr] = right_i;
                            j_stack[stack_ptr] = j;
                            
                            //split cluster i
                            ++stack_ptr;
                            i_stack[stack_ptr] = left_i;
                            j_stack[stack_ptr] = j;
                        }
                        else //score_i < score_j
                        {
                            //split cluster j
                            ++stack_ptr;
                            i_stack[stack_ptr] = i;
                            j_stack[stack_ptr] = right_j;
                            
                            ++stack_ptr;
                            i_stack[stack_ptr] = i;
                            j_stack[stack_ptr] = left_j;
                        }
                    }
                }
                else
                {
                    collision_time.push_back(a * t_init);
                    
                    if constexpr ( is_symmetric )
                    {
                        // For creating an upper triangle matrix we store the pair  { min(i,j), max(i,j) }.
                        if (i <= j)
                        {
                            collision_idx.push_back(i);
                            collision_jdx.push_back(j);
                        }
                        else
                        {
                            collision_idx.push_back(j);
                            collision_jdx.push_back(i);
                        }
                    }
                    else
                    {
                        // No symmetry exploited.
                        collision_idx.push_back(i);
                        collision_jdx.push_back(j);
                    }
                }
            }
        }
    } // FindCollidingClusters_DFS
