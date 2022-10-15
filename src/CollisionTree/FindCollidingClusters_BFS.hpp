protected:

    void FindCollidingClusters_BFS( const Int max_leaves ) const
    {
        i_queue.push_back(null);
        j_queue.push_back(null);
        
        std::vector<Int>   & collision_idx  = collision_i[0];
        std::vector<Int>   & collision_jdx  = collision_j[0];
        std::vector<SReal> & collision_time = collision_t[0];

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

        while( !i_queue.empty() && ( static_cast<Int>(i_queue.size()) < max_leaves ) )
        {
            const Int i = i_queue.front();
            const Int j = j_queue.front();
            i_queue.pop_front();
            j_queue.pop_front();
            
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
                                            
                    const SReal score_i = (left_i>=null) * S_C_ser(i,0);
                    const SReal score_j = (left_j>=null) * T_C_ser(j,0);
                    
                    if( score_i == score_j )
                    {
                        // tie breaker: split both clusters

                        if constexpr ( is_symmetric )
                        {
                            if(i == j)                            {
                                //  Creating 3 blockcluster children, since there is one block is just the mirror of another one.
                                
                                i_queue.push_back(left_i);
                                j_queue.push_back(right_j);
                                
                                i_queue.push_back(right_i);
                                j_queue.push_back(right_j);
                                
                                i_queue.push_back(left_i);
                                j_queue.push_back(left_j);
                            }
                            else
                            {
                                // In case of exploit_symmetry !=0, this is a very seldom case; still requird to preserve symmetry.
                                // This happens only if i and j represent_diffent clusters with same radii.
                                
                                i_queue.push_back(right_i);
                                j_queue.push_back(right_j);
                                
                                i_queue.push_back(left_i);
                                j_queue.push_back(right_j);
                                
                                i_queue.push_back(right_i);
                                j_queue.push_back(left_j);
                                
                                i_queue.push_back(left_i);
                                j_queue.push_back(left_j);
                            }
                        }
                        else
                        {
                            i_queue.push_back(right_i);
                            j_queue.push_back(right_j);
                            
                            i_queue.push_back(left_i);
                            j_queue.push_back(right_j);
                            
                            i_queue.push_back(right_i);
                            j_queue.push_back(left_j);
                            
                            i_queue.push_back(left_i);
                            j_queue.push_back(left_j);
                        }
                    }
                    else
                    {
                        // split only the larger cluster
                        if (score_i > score_j)
                        {
                            i_queue.push_back(right_i);
                            j_queue.push_back(j);
                            
                            //split cluster i
                            i_queue.push_back(left_i);
                            j_queue.push_back(j);
                        }
                        else //score_i < score_j
                        {
                            //split cluster j
                            i_queue.push_back(i);
                            j_queue.push_back(right_j);
                            
                            i_queue.push_back(i);
                            j_queue.push_back(left_j);
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
    } // FindCollidingClusters_BFS
