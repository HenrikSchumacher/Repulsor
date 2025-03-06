public:

    TOOLS_FORCE_FLATTENING void Traverse_DepthFirst(
        const Int thread, const Int i0, const Int j0
    )
    {
        TOOLS_DEBUG_PRINT(className()+"::Traverse_DepthFirst("+ToString(thread)+")...");
        
        mref<Kernel_T> K = kernels[static_cast<Size_T>(thread)];
        
        constexpr Int stack_size = 4 * max_depth;
        
        Int stack [stack_size][2] = {};
        Int stack_ptr = -1;
        
        // Helper routine to manage the pair_stack.
        auto push = [&stack,&stack_ptr]( const Int i, const Int j )
        {
            ++stack_ptr;
            stack[stack_ptr][0] = i;
            stack[stack_ptr][1] = j;
        };
        
        // Helper routine to manage the pair_stack.
        auto pop = [&stack,&stack_ptr]()
        {
            auto result = std::pair( stack[stack_ptr][0], stack[stack_ptr][1] );
            stack_ptr--;
            return result;
        };
        
        auto continueQ = [&stack_ptr,this,thread]()
        {
            const bool overflowQ = (stack_ptr >= stack_size - 4);
            
            if( (0 <= stack_ptr) && (!overflowQ) ) [[likely]]
            {
                return true;
            }
            else
            {
                if ( overflowQ ) [[unlikely]]
                {
                    eprint(this->className()+"::Traverse_DepthFirst("+ToString(thread)+": Stack overflow.");
                }
                return false;
            }
        };
        
        push(i0,j0);
        
        while( continueQ() )
        {
            auto [i,j] = pop();
            
            K.LoadClusterS(i);
            K.LoadClusterT(j);
            
            const bool admissableQ = (!( symmetricQ && (i == j) )) && K.AdmissableQ();
            
            if( !admissableQ )
            {
                const Int L_i = S_C_L[i];
                const Int L_j = T_C_L[j];
                
                // Using that children are either both interiors or both leaves.
                if( (L_i >= null) || (L_j >= null) )
                {
                    const Int R_i = S_C_R[i];
                    const Int R_j = T_C_R[j];
                    
                    const SReal score_i = (L_i>=null) * K.ClusterScoreS();
                    const SReal score_j = (L_j>=null) * K.ClusterScoreT();
                    
                    if( score_i == score_j && score_i > zero )
                    {
                        if constexpr ( symmetricQ )
                        {
                            if( i == j )
                            {
                                //  Creating 3 blockcluster children, since there is one block that is just the mirror of another one.
                                push( L_i, R_j );
                                push( R_i, R_j );
                                push( L_i, L_j );
                            }
                            else
                            {
                                // This is a very seldom case; still required to preserve symmetry.
                                // This happens only if i and j represent_diffent clusters with same radii.
                                push( R_i, R_j );
                                push( L_i, R_j );
                                push( R_i, L_j );
                                push( L_i, L_j );
                            }
                        }
                        else
                        {
                            // Split both clusters
                            push( R_i, R_j );
                            push( L_i, R_j );
                            push( R_i, L_j );
                            push( L_i, L_j );
                        }
                    }
                    else
                    {
                        // split only larger cluster
                        if( score_i > score_j )
                        {
                            push( R_i, j );
                            push( L_i, j );
                        }
                        else //score_i < score_j
                        {
                            //split cluster j
                            push( i, R_j );
                            push( i, L_j );
                        }
                    }
                }
                else // (L_i < null) && (L_j < null)
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
        
        TOOLS_DEBUG_PRINT(className()+"::Traverse_DepthFirst("+ToString(thread)+") done.");
        
    } // Traverse_DepthFirst    
    
