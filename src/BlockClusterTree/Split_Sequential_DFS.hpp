protected:

    void Split_Sequential_DFS( const Int i0, const Int j0 )
    {
//            ptic(className()+"::Split_Sequential_DFS");
        
        const Int thread = omp_get_thread_num();
        
        const Int max_depth = 128;
        Int i_stack[128] = {};
        Int j_stack[128] = {};

        Int stack_ptr = 0;
        i_stack[0] = i0;
        j_stack[0] = j0;

        const SparseBinaryMatrixCSR<Int> & A = S.PrimitiveAdjacencyMatrix();
        const Real intersection_theta2 = NearFieldSeparationParameter() * NearFieldSeparationParameter();
        
        std::vector<Int> &    inter_idx =    inter_i[thread];
        std::vector<Int> &    inter_jdx =    inter_j[thread];
        std::vector<Int> & verynear_idx = verynear_i[thread];
        std::vector<Int> & verynear_jdx = verynear_j[thread];
        std::vector<Int> &     near_idx =     near_i[thread];
        std::vector<Int> &     near_jdx =     near_j[thread];
        std::vector<Int> &      far_idx =      far_i[thread];
        std::vector<Int> &      far_jdx =      far_j[thread];

        
        std::shared_ptr<BoundingVolume_T> S_C_proto_ptr = S.ClusterPrototype().Clone();
        std::shared_ptr<BoundingVolume_T> T_C_proto_ptr = T.ClusterPrototype().Clone();
        std::shared_ptr<Primitive_T>      S_P_proto_ptr = S.PrimitivePrototype().Clone();
        std::shared_ptr<Primitive_T>      T_P_proto_ptr = T.PrimitivePrototype().Clone();
        
        BoundingVolume_T & S_C_proto = *S_C_proto_ptr;
        BoundingVolume_T & T_C_proto = *T_C_proto_ptr;
        Primitive_T      & S_P_proto = *S_P_proto_ptr;
        Primitive_T      & T_P_proto = *T_P_proto_ptr;
        GJK_T G;
        
        const Int   * restrict const S_C_left       = S.ClusterLeft().data();
        const Int   * restrict const S_C_right      = S.ClusterRight().data();
        const Int   * restrict const S_C_begin      = S.ClusterBegin().data();
        const Int   * restrict const S_C_end        = S.ClusterEnd().data();
        
        const Int   * restrict const T_C_left       = T.ClusterLeft().data();
        const Int   * restrict const T_C_right      = T.ClusterRight().data();
        const Int   * restrict const T_C_begin      = T.ClusterBegin().data();
        const Int   * restrict const T_C_end        = T.ClusterEnd().data();
        
              SReal * restrict const S_C_serialized = S.ClusterSerializedData().data();
              SReal * restrict const T_C_serialized = T.ClusterSerializedData().data();
    
              SReal * restrict const S_P_serialized = S.PrimitiveSerializedData().data();
              SReal * restrict const T_P_serialized = T.PrimitiveSerializedData().data();
        
        while( (zero <= stack_ptr) && (stack_ptr < max_depth) )
        {
            const Int i = i_stack[stack_ptr];
            const Int j = j_stack[stack_ptr];
            stack_ptr--;
            
            S_C_proto.SetPointer(S_C_serialized,i);
            T_C_proto.SetPointer(T_C_serialized,j);
            
            const bool separatedQ = (!( is_symmetric && (i == j) )) && G.MultipoleAcceptanceCriterion( S_C_proto, T_C_proto, far_theta2 );
            
            
            if( !separatedQ )
            {
                const Int left_i = S_C_left[i];
                const Int left_j = T_C_left[j];

                // Warning: This assumes that either both children are defined or empty.
                if( left_i >= zero || left_j >= zero )
                {

                    const Int right_i = S_C_right[i];
                    const Int right_j = T_C_right[j];
                    
                    const SReal score_i = (left_i>=zero) * S_C_proto.SquaredRadius();
                    const SReal score_j = (left_j>=zero) * T_C_proto.SquaredRadius();

                    if( score_i == score_j && score_i > static_cast<SReal>(0) /*&& score_j > static_cast<SReal>(0)*/ )
                    {
                        // tie breaker: split both clusters

                        if ( is_symmetric )
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
                    
                    if( is_symmetric )
                    {
                        if( i == j )
                        {
                            // This is a diagonal block.
                            // Only traverse uppper triangle part.
                            const Int ii_begin = S_C_begin[i];
                            const Int ii_end   = S_C_end[i];
                            
                            for( Int ii = ii_begin; ii < ii_end; ++ii )
                            {
                                // We need the diagonals of the near field matrix to be part of the sparsity structure.
                                inter_idx.push_back(ii);
                                inter_jdx.push_back(ii);
                                
                                for( Int jj = ii+1; jj < ii_end; ++jj )
                                {
                                    const bool neighbor_found = A.FindNonzeroPosition(ii,jj) >= 0;
                                    
                                    S_P_proto.SetPointer(S_P_serialized,ii);
                                    T_P_proto.SetPointer(T_P_serialized,jj);

                                    const bool admissable = neighbor_found || G.MultipoleAcceptanceCriterion( S_P_proto, T_P_proto, near_theta2);
                                    
                                    if( admissable )
                                    {
                                        near_idx.push_back(ii);
                                        near_jdx.push_back(jj);
                                    }
                                    else
                                    {
                                        const bool intersecting = (!G.SeparatedQ()) && G.IntersectingQ( S_P_proto, T_P_proto, intersection_theta2 );
                                        
                                        if( intersecting )
                                        {
                                            inter_idx.push_back(ii);
                                            inter_jdx.push_back(jj);
                                        }
                                        else
                                        {
                                            verynear_idx.push_back(ii);
                                            verynear_jdx.push_back(jj);
                                        }
                                    }
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
                                for( Int jj = jj_begin; jj < jj_end; ++jj )
                                {
                                    // No primitive can be member in more than one leaf cluster.
                                    assert( ii != jj );
                                    
                                    Int ii_;
                                    Int jj_;
                                    
                                    if( ii < jj )
                                    {
                                        ii_ = ii;
                                        jj_ = jj;
                                    }
                                    else
                                    {
                                        ii_ = jj;
                                        jj_ = ii;
                                    }
                                    
                                    const bool neighbor_found = A.FindNonzeroPosition(ii,jj) >= 0;
                                    
                                    S_P_proto.SetPointer(S_P_serialized,ii);
                                    T_P_proto.SetPointer(T_P_serialized,jj);

                                    const bool admissable = neighbor_found || G.MultipoleAcceptanceCriterion( S_P_proto, T_P_proto, near_theta2);
                                    
                                    if( admissable )
                                    {
                                        near_idx.push_back(ii_);
                                        near_jdx.push_back(jj_);
                                    }
                                    else
                                    {
                                        const bool intersecting = (!G.SeparatedQ()) && G.IntersectingQ( S_P_proto, T_P_proto, intersection_theta2 );
                                        
                                        if( intersecting )
                                        {
                                            inter_idx.push_back(ii_);
                                            inter_jdx.push_back(jj_);
                                        }
                                        else
                                        {
                                            verynear_idx.push_back(ii_);
                                            verynear_jdx.push_back(jj_);
                                        }
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
                            for( Int jj = jj_begin; jj < jj_end; ++jj )
                            {
                                S_P_proto.SetPointer(S_P_serialized,ii);
                                T_P_proto.SetPointer(T_P_serialized,jj);
                                
                                const bool admissable = G.MultipoleAcceptanceCriterion( S_P_proto, T_P_proto, near_theta2);
                                
                                if( admissable )
                                {
                                    near_idx.push_back(ii);
                                    near_jdx.push_back(jj);
                                }
                                else
                                {
                                    const bool intersecting = (!G.SeparatedQ()) && G.IntersectingQ( S_P_proto, T_P_proto, intersection_theta2 );
                                    
                                    if( intersecting )
                                    {
                                        inter_idx.push_back(ii);
                                        inter_jdx.push_back(jj);
                                    }
                                    else
                                    {
                                        verynear_idx.push_back(ii);
                                        verynear_jdx.push_back(jj);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else // separatedQ
            {
                //create far field leaf blockcluster
                if( is_symmetric )
                {
                    // For creating an upper triangle matrix we store the pair  { min(i,j), max(i,j) }.
                    if (i <= j)
                    {
                        far_idx.push_back(i);
                        far_jdx.push_back(j);
                    }
                    else
                    {
                        far_idx.push_back(j);
                        far_jdx.push_back(i);
                    }
                }
                else
                {
                    // No symmetry exploited.
                    far_idx.push_back(i);
                    far_jdx.push_back(j);
                }
            }
        }
//            ptoc(className()+"::Split_Sequential_DFS");
        
    } // Split_Sequential_DFS
