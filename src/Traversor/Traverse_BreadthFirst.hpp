public:

    void Traverse_BreadthFirst( const Int i0, const Int j0, const Int max_leaves )
    {
        i_queue = std::deque<Int>();
        j_queue = std::deque<Int>();
        i_queue.push_back(static_cast<Int>(i0));
        j_queue.push_back(static_cast<Int>(j0));
        
        const Int thread = omp_get_thread_num();
        
        Kernel_T K = kernels[thread];
        
        while( !i_queue.empty() && ( static_cast<Int>(i_queue.size()) < max_leaves ) )
        {
            const Int i = i_queue.front();
            const Int j = j_queue.front();
            i_queue.pop_front();
            j_queue.pop_front();
            
            K.LoadClusterS(i);
            K.LoadClusterT(j);
            
            const bool admissable = (!( is_symmetric && (i == j) )) && K.IsAdmissable();

            if( !admissable )
            {
                const Int left_i = S_C_left[i];
                const Int left_j = T_C_left[j];

                // Warning: This assumes that either both children are defined or empty.
                if( left_i >= zero || left_j >= zero )
                {
                    const Int right_i = S_C_right[i];
                    const Int right_j = T_C_right[j];
                    
                    const SReal score_i = (left_i>=zero) * K.ScoreS();
                    const SReal score_j = (left_j>=zero) * K.ScoreT();

                    if( score_i == score_j )
                    {
                        // tie breaker: split both clusters

                        if constexpr ( is_symmetric )
                        {
                            if( i == j )
                            {
                                //  Creating 3 blockcluster children, since there is one block that is just the mirror of another one.
                                
                                i_queue.push_back(left_i);
                                j_queue.push_back(right_j);
                                
                                i_queue.push_back(right_i);
                                j_queue.push_back(right_j);
                                
                                i_queue.push_back(left_i);
                                j_queue.push_back(left_j);
                            }
                            else
                            {
                                // This is a very seldom case; still required to preserve symmetry.
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
                            // Split both clusters
                            
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
                        // split only larger cluster
                        if( score_i > score_j )
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
                else // left_i < zero && left_j < zero
                {
                    // We know that i and j are leaf clusters and that they belong either to the near field, the very near field or contain intersecting primitives.
                    
                    // We have to go through all the primitive pairs to classify them.
                    
                    if constexpr ( is_symmetric )
                    {
                        if( i == j )
                        {
                            // This is a diagonal block.
                            // Only traverse uppper triangle part.
                            const Int ii_begin = S_C_begin[i];
                            const Int ii_end   = S_C_end[i];
                            
                            for( Int ii = ii_begin; ii < ii_end; ++ii )
                            {
                                K.LoadPrimitiveS(ii);
                                
                                K.ComputeLeafDiagonal();
                                
                                for( Int jj = ii+1; jj < ii_end; ++jj )
                                {
                                    K.LoadPrimitiveT(jj);
                                    
                                    K.ComputeLeaf();
                                }
                            }
                        }
                        else // i != j
                        {
                            // Traverse whole block, but sort primitive indices.
                            const Int ii_begin = S_C_begin[i];
                            const Int ii_end   = S_C_end[i];
                            const Int jj_begin = T_C_begin[j];
                            const Int jj_end   = T_C_end[j];
                            
                            for( Int ii = ii_begin; ii < ii_end; ++ii )
                            {
                                K.LoadPrimitiveS(ii);
                                
                                for( Int jj = jj_begin; jj < jj_end; ++jj )
                                {
                                    
                                    K.LoadPrimitiveT(jj);
                                    
                                    // No primitive can be member in more than one leaf cluster.
                                    assert( ii != jj );

                                    if( ii < jj )
                                    {
                                        K.ComputeLeaf();
                                    }
                                    else
                                    {
                                        K.ComputeLeafSwapped();
                                    }
                                }
                            }
                        }
                    }
                    else // !is_symmetric
                    {
                        // Just traverse the whole block.
                        const Int ii_begin = S_C_begin[i];
                        const Int ii_end   = S_C_end[i];
                        const Int jj_begin = T_C_begin[j];
                        const Int jj_end   = T_C_end[j];
                        
                        for( Int ii = ii_begin; ii < ii_end; ++ii )
                        {
                            K.LoadPrimitiveS(ii);
                            
                            for( Int jj = jj_begin; jj < jj_end; ++jj )
                            {
                                K.LoadPrimitiveT(jj);
                                
                                K.ComputeLeaf();
                            }
                        }
                    }
                }
            }
            else // admissable
            {
                if constexpr ( is_symmetric )
                {
                    if( i <= j )
                    {
                        K->ComputeAdmissable();
                    }
                    else
                    {
                        K->ComputeAdmissableSwapped();
                    }
                }
                else
                {
                    // No symmetry exploited.
                    K->ComputeAdmissable();
                }
            }
        }
        
    } // Traverse_BreadthFirst
