public:

    TOOLS_FORCE_FLATTENING void Traverse_BreadthFirst(
        const Int thread, const Int i0, const Int j0, const Int max_leaves
    )
    {
        TOOLS_PTIC(className()+"::Traverse_BreadthFirst");
        
        mref<Kernel_T> K = kernels[static_cast<Size_T>(thread)];
        
        queue = std::deque<std::pair<Int,Int>>();
        
        auto push = [this]( const Int i, const Int j )
        {
            queue.push_back( std::pair(i,j) );
        };
        
        push(i0,j0);
        
        while( !queue.empty() && ( static_cast<Int>(queue.size()) < max_leaves ) )
        {
            auto [i,j] = queue.front();
            queue.pop_front();
            
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
                        // tie breaker: split both clusters

                        if constexpr ( symmetricQ )
                        {
                            if( i == j )
                            {
                                //  Creating 3 blockcluster children, since there is one block that is just the mirror of another one.
                                push(L_i,R_j);
                                push(R_i,R_j);
                                push(L_i,L_j);
                            }
                            else
                            {
                                // This is a very seldom case; still required to preserve symmetry.
                                // This happens only if i and j represent_diffent clusters with same radii.
                                push(R_i,R_j);
                                push(L_i,R_j);
                                push(R_i,L_j);
                                push(L_i,L_j);
                            }
                        }
                        else
                        {
                            // Split both clusters
                            push(R_i,R_j);
                            push(L_i,R_j);
                            push(R_i,L_j);
                            push(L_i,L_j);
                        }
                    }
                    else
                    {
                        // split only larger cluster
                        if( score_i > score_j )
                        {
                            push(R_i,j);
                            push(L_i,j);
                        }
                        else //score_i < score_j
                        {
                            //split cluster j
                            push(i,R_j);
                            push(i,L_j);
                        }
                    }
                }
                else // L_i < null && L_j < null
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
                        }
                        else // i != j
                        {
                            // Traverse whole block, but sort primitive indices.
                            // Traverse whole block, but sort primitive indices.
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
                                        
                                        // No primitive can be member in more than one leaf cluster.
                                        
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
        }
        
        TOOLS_PTOC(className()+"::Traverse_BreadthFirst");
        
    } // Traverse_BreadthFirst
