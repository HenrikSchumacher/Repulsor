public:

    void Traverse_DepthFirst( const Int i0, const Int j0 )
    {
        const Int thread = omp_get_thread_num();
        
        Kernel_T ker = kernels[thread];
        
        Int i_stack[max_depth] = {};
        Int j_stack[max_depth] = {};

        Int stack_ptr = 0;
        i_stack[0] = i0;
        j_stack[0] = j0;
        
        while( (zero <= stack_ptr) && (stack_ptr < max_depth) )
        {
            const Int i = i_stack[stack_ptr];
            const Int j = j_stack[stack_ptr];
            
            stack_ptr--;
            
            ker.LoadClusterS(i);
            ker.LoadClusterT(j);
            
            const bool admissable = (!( is_symmetric && (i == j) )) && ker.IsAdmissable();
            
            if( !admissable )
            {
                const Int left_i = S_C_left[i];
                const Int left_j = T_C_left[j];

                // Warning: This assumes that either both children are defined or empty.
                if( left_i >= zero || left_j >= zero )
                {
                    const Int right_i = S_C_right[i];
                    const Int right_j = T_C_right[j];
                    
                    const SReal score_i = (left_i>=zero) * ker.ClusterScoreS();
                    const SReal score_j = (left_j>=zero) * ker.ClusterScoreT();

                    if( score_i == score_j )
                    {
                        // tie breaker: split both clusters

                        if constexpr ( is_symmetric )
                        {
                            if( i == j )
                            {
                                //  Creating 3 blockcluster children, since there is one block that is just the mirror of another one.
                                
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
                            // Split both clusters
                            
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
                        if( score_i > score_j )
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
                else // left_i < zero && left_j < zero
                {
                    // We know that i and j are leaf clusters and that they belong either to the near field, the very near field or contain intersecting primitives.
                    
                    // We have to go through all the primitive pairs to classify them.
                    
                    if constexpr ( is_symmetric )
                    {
                        if( i == j )
                        {
                            if constexpr ( leafs_are_singletons )
                            {
                                ker.LoadPrimitiveS( S_C_begin[i] );
                                ker.ComputeLeafDiagonal();
                            }
                            else
                            {
                                // This is a diagonal block.
                                // Only traverse uppper triangle part.
                                const Int ii_begin = S_C_begin[i];
                                const Int ii_end   = S_C_end[i];
                                
                                for( Int ii = ii_begin; ii < ii_end; ++ii )
                                {
                                    ker.LoadPrimitiveS(ii);
                                    
                                    ker.ComputeLeafDiagonal();
                                    
                                    for( Int jj = ii+1; jj < ii_end; ++jj )
                                    {
                                        ker.LoadPrimitiveT(jj);
                                        
                                        ker.ComputeLeaf();
                                    }
                                }
                            }
                        }
                        else // i != j
                        {
                            // Traverse whole block, but sort primitive indices.
                            if constexpr ( leafs_are_singletons )
                            {
                                ker.LoadPrimitiveS( S_C_begin[i] );
                                ker.LoadPrimitiveS( T_C_begin[j] );
                                
                                if( ii < jj )
                                {
                                    ker.ComputeLeaf();
                                }
                                else
                                {
                                    ker.ComputeLeafSwapped();
                                }
                            }
                            else
                            {
                                const Int ii_begin = S_C_begin[i];
                                const Int ii_end   = S_C_end[i];
                                const Int jj_begin = T_C_begin[j];
                                const Int jj_end   = T_C_end[j];
                                
                                for( Int ii = ii_begin; ii < ii_end; ++ii )
                                {
                                    ker.LoadPrimitiveS(ii);
                                    
                                    for( Int jj = jj_begin; jj < jj_end; ++jj )
                                    {
                                        ker.LoadPrimitiveT(jj);
                                        
                                        if( ii < jj )
                                        {
                                            ker.ComputeLeaf();
                                        }
                                        else
                                        {
                                            ker.ComputeLeafSwapped();
                                        }
                                    }
                                }
                            }
                        }
                    }
                    else // !is_symmetric
                    {
                        if constexpr ( leafs_are_singletons )
                        {
                            ker.LoadPrimitiveS( S_C_begin[i] );
                            ker.LoadPrimitiveS( T_C_begin[j] );
                            ker.ComputeLeaf();
                        }
                        else
                        {
                            // Just traverse the whole block.
                            const Int ii_begin = S_C_begin[i];
                            const Int ii_end   = S_C_end[i];
                            const Int jj_begin = T_C_begin[j];
                            const Int jj_end   = T_C_end[j];
                            
                            for( Int ii = ii_begin; ii < ii_end; ++ii )
                            {
                                ker.LoadPrimitiveS(ii);
                                
                                for( Int jj = jj_begin; jj < jj_end; ++jj )
                                {
                                    ker.LoadPrimitiveT(jj);
                                    
                                    ker.ComputeLeaf();
                                }
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
                        ker->ComputeAdmissable();
                    }
                    else
                    {
                        ker->ComputeAdmissableSwapped();
                    }
                }
                else
                {
                    // No symmetry exploited.
                    ker->ComputeAdmissable();
                }
            }
        }
    } // Traverse_DepthFirst
    
    
