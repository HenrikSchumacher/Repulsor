protected:

    SReal MaximumSafeStepSize_DFS( const Int i_0, const Int j_0, const SReal t_init_ ) const
    {
        Int i_stack[max_depth] = {};
        Int j_stack[max_depth] = {};
        
        Int stack_ptr = null;
        i_stack[0] = i_0;
        j_stack[0] = j_0;
        
        const Int * restrict const S_left     = S.ClusterLeft().data();
        const Int * restrict const S_right    = S.ClusterRight().data();
        const Int * restrict const T_left     = T.ClusterLeft().data();
        const Int * restrict const T_right    = T.ClusterRight().data();
        
        const Int * restrict const S_begin    = S.ClusterBegin().data();
        const Int * restrict const S_end      = S.ClusterEnd().data();
        const Int * restrict const T_begin    = T.ClusterBegin().data();
        const Int * restrict const T_end      = T.ClusterEnd().data();
        
        const Tensor2<SReal,Int> & S_C_ser    = S.ClusterSerializedData();
        const Tensor2<SReal,Int> & S_C_up_ser = S.ClusterUpdatedSerializedData();
        const Tensor2<SReal,Int> & T_C_ser    = T.ClusterSerializedData();
        const Tensor2<SReal,Int> & T_C_up_ser = T.ClusterUpdatedSerializedData();
        
        const Tensor2<SReal,Int> & S_P_ser    = S.PrimitiveSerializedData();
        const Tensor2<SReal,Int> & S_P_v_ser  = S.PrimitiveVelocitiesSerializedData();
        const Tensor2<SReal,Int> & T_P_ser    = T.PrimitiveSerializedData();
        const Tensor2<SReal,Int> & T_P_v_ser  = T.PrimitiveVelocitiesSerializedData();

        const SparseBinaryMatrixCSR<Int> & A = S.PrimitiveAdjacencyMatrix();
        
        PrimitiveCollisionFinder_T C ( S.MovingPrimitivePrototype(), T.MovingPrimitivePrototype() );
        
        SReal a = static_cast<SReal>(0);
        SReal b = static_cast<SReal>(1);
        
        SReal t_max = t_init_;
        
        while( (null <= stack_ptr) && (stack_ptr < max_depth - 4) )
        {
            const Int i = i_stack[stack_ptr];
            const Int j = j_stack[stack_ptr];
            stack_ptr--;
            
            Compute_AABB_CollisionTimeInterval<AMB_DIM, SReal, Int> (
                S_C_ser.data(i), S_C_up_ser.data(i),
                T_C_ser.data(j), T_C_up_ser.data(j),
                a, b
            );
            
            if( a * t_init_ < t_max )
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
                    
                    if( (score_i == score_j) && (score_i > static_cast<SReal>(0)) /*&& score_j > static_cast<SReal>(0)*/ )
                    {
                        // tie breaker: split both clusters
                        
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
                                // In case of is_symmetric !=0, this is a very seldom case; still requird to preserve symmetry.
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
                        else // !symmetric
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
                    else // score_i != score_j
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
                else
                {
                    if constexpr ( is_symmetric )
                    {
                        if( i == j )
                        {
                            const Int ii_begin = S_begin[i];
                            const Int ii_end   = S_end  [i];
                            
                            for( Int ii = ii_begin; ii < ii_end; ++ii )
                            {
                                for( Int jj = ii+1; jj < ii_end; ++jj )
                                {
                                    if( A.FindNonzeroPosition(ii,jj) < null )
                                    {
                                        const SReal t = C.FindMaximumSafeStepSize(
                                            S_P_ser.data(ii), S_P_v_ser.data(ii),
                                            T_P_ser.data(jj), T_P_v_ser.data(jj),
                                            t_max, false
                                        );
                                        
                                        t_max = std::min(t, t_max);
                                    }
                                    
                                } // for( Int jj = jj_begin; jj < jj_end; ++jj )
                                
                            } // for( Int ii = ii_begin; ii < ii_end; ++ii )
                        }
                        else
                        {
                            const Int ii_begin = S_begin[i];
                            const Int ii_end   = S_end[i];
                            
                            const Int jj_begin = T_begin[j];
                            const Int jj_end   = T_end[j];
                            
                            for( Int ii = ii_begin; ii < ii_end; ++ii )
                            {
                                for( Int jj = jj_begin; jj < jj_end; ++jj )
                                {
                                    if( A.FindNonzeroPosition(ii,jj) < null )
                                    {
                                        const SReal t = C.FindMaximumSafeStepSize(
                                            S_P_ser.data(ii), S_P_v_ser.data(ii),
                                            T_P_ser.data(jj), T_P_v_ser.data(jj),
                                            t_max, false
                                        );
                                        
                                        t_max = std::min(t, t_max);
                                    }
                                    
                                } // for( Int jj = jj_begin; jj < jj_end; ++jj )
                                
                            } // for( Int ii = ii_begin; ii < ii_end; ++ii )
                        }
                    }
                    else
                    {
                        const Int ii_begin = S_begin[i];
                        const Int ii_end   = S_end  [i];
                        
                        const Int jj_begin = T_begin[j];
                        const Int jj_end   = T_end[j];
                        
                        for( Int ii = ii_begin; ii < ii_end; ++ii )
                        {
                            for( Int jj = jj_begin; jj < jj_end; ++jj )
                            {
                                const SReal t = C.FindMaximumSafeStepSize(
                                    S_P_ser.data(ii), S_P_v_ser.data(ii),
                                    T_P_ser.data(jj), T_P_v_ser.data(jj),
                                    t_max, false
                                );
                                
                                t_max = std::min(t, t_max);
                                
                            } // for( Int jj = jj_begin; jj < jj_end; ++jj )
                            
                        } // for( Int ii = ii_begin; ii < ii_end; ++ii )
                    }
                
                } // if( left_i >= null || left_j >= null )
            }
        
        } // while
        
        
        return t_max;
        
    } // MaximumSafeStepSize_DFS
    
