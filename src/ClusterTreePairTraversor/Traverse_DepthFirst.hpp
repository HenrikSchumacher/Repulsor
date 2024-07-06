public:

    force_flattening void Traverse_DepthFirst(
        const Int thread, const Int i0, const Int j0
    )
    {
        debug_print(className()+"::Traverse_DepthFirst("+ToString(thread)+")...");
        
        mref<Kernel_T> K = kernels[thread];
        
        const Int stack_size = 4 * max_depth;
        
        Int i_stack[stack_size] = {};
        Int j_stack[stack_size] = {};

        Int stack_ptr = null;
        i_stack[0] = i0;
        j_stack[0] = j0;
        
        while( (null <= stack_ptr) && (stack_ptr < stack_size-4) )
        {
            const Int i = i_stack[stack_ptr];
            const Int j = j_stack[stack_ptr];
            --stack_ptr;
            
            K.LoadClusterS(i);
            K.LoadClusterT(j);
            
            const bool admissableQ = (!( symmetricQ && (i == j) )) && K.AdmissableQ();
            
            if( !admissableQ )
            {
                const Int left_i = S_C_left[i];
                
                // Using that children are either both interiors or both leaves.
                if( (left_i >= null) /*|| (left_j >= null)*/ )
                {
                    const Int left_j = T_C_left[j];
                    const Int right_i = S_C_right[i];
                    const Int right_j = T_C_right[j];
                    
                    const SReal score_i = (left_i>=null) * K.ClusterScoreS();
                    const SReal score_j = (left_j>=null) * K.ClusterScoreT();
                    
                    if( score_i == score_j && score_i > zero )
                    {
                        if constexpr ( symmetricQ )
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
                else // (left_i < null) && (left_j < null)
                {
                    // We know that i and j are leaf clusters and that they belong either to the near field, the very near field or contain intersecting primitives.
                    
                    // We have to go through all the primitive pairs to classify them.
                    
                    if constexpr ( symmetricQ )
                    {
                        if( i == j )
                        {
                            if constexpr ( leaves_are_singletons )
                            {
                                K.LoadPrimitiveS( S_C_begin[i] );
                                K.ComputeLeafDiagonal();
                            }
                            else
                            {
                                // This is a diagonal block.
                                // Only traverse uppper triangle part.
                                const Int ii_begin = S_C_begin[i];
                                const Int ii_end   = S_C_end  [i];
                                
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
                        }
                        else // i != j
                        {
                            if constexpr ( leaves_are_singletons )
                            {
                                const Int ii = S_C_begin[i];
                                const Int jj = T_C_begin[j];
                                
                                K.LoadPrimitiveS( ii );
                                K.LoadPrimitiveT( jj );
                                
                                if( ii < jj )
                                {
                                    K.ComputeLeaf();
                                }
                                else
                                {
                                    K.ComputeLeafSwapped();
                                }
                            }
                            else
                            {
                                // Traverse whole block, but sort primitive indices.
                                
                                const Int ii_begin = S_C_begin[i];
                                const Int ii_end   = S_C_end  [i];
                                const Int jj_begin = T_C_begin[j];
                                const Int jj_end   = T_C_end  [j];
                                
                                for( Int ii = ii_begin; ii < ii_end; ++ii )
                                {
                                    K.LoadPrimitiveS(ii);
                                    
                                    for( Int jj = jj_begin; jj < jj_end; ++jj )
                                    {
                                        K.LoadPrimitiveT(jj);
                                        
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
                    }
                    else // !symmetricQ
                    {
                        if constexpr ( leaves_are_singletons )
                        {
                            K.LoadPrimitiveS( S_C_begin[i] );
                            K.LoadPrimitiveT( T_C_begin[j] );
                            K.ComputeLeaf();
                        }
                        else
                        {
                            // Just traverse the whole block.
                            const Int ii_begin = S_C_begin[i];
                            const Int ii_end   = S_C_end  [i];
                            const Int jj_begin = T_C_begin[j];
                            const Int jj_end   = T_C_end  [j];
                            
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
            }
            else // admissableQ
            {
                if constexpr ( symmetricQ )
                {
                    if( i <= j )
                    {
                        K.ComputeAdmissable();
                    }
                    else
                    {
                        K.ComputeAdmissableSwapped();
                    }
                }
                else
                {
                    // No symmetry exploited.
                    K.ComputeAdmissable();
                }
            }
            
        } // while
        
        if( stack_ptr >= stack_size-4 )
        {
            eprint( ClassName() + "::Traverse_DepthFirst: Stack overflow." );
        }
        
        debug_print(className()+"::Traverse_DepthFirst("+ToString(thread)+") done.");
        
    } // Traverse_DepthFirst
    
    
