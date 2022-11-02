public:

    force_inline Real Compute( const LInt k_global )
    {
        Real sum = static_cast<Real>(0);
        
//        ++primitive_count;
        
        bool from_above = true;
        bool shall_continue = true;
        
        S_Tree.ToChild(0);
        T_Tree.ToChild(0);
    
        while( shall_continue )
        {
            if( from_above )
            {
                if( S_Tree.Level() >= this->max_level )
                {
                    // If at lowest level and inadmissable then we just compute the energy and move up.
//                    max_level_reached = max_level;
//                    block_count++;
//                    bottom_count++;
                    sum += compute();
                    S_Tree.ToParent();
                    T_Tree.ToParent();
                    from_above = false;
                }
                else
                {
                    // If not at lowest level, then we have to check for admissability.
                    auto & P = S_Tree.SimplexPrototype();
                    auto & Q = T_Tree.SimplexPrototype();
                    
                    const bool admissable = this->gjk.MultipoleAcceptanceCriterion(P, Q, this->theta2);
                    
                    if( admissable )
                    {
                        // We compute energy, go to parent, and prepare the next child of the parent.
//                        max_level_reached = std::max( max_level_reached, S_Tree.Level() );
//                        block_count++;
                        sum += compute();
                        S_Tree.ToParent();
                        T_Tree.ToParent();
                        from_above = false;
                    }
                    else
                    {
                        // If inadmissable, we go a level deeper.

                        S_Tree.ToChild(0);
                        T_Tree.ToChild(0);
                        from_above = true;
                    }
                }
            }
            else
            {
                // If we come from below, we have to find the next pair of simplices to visit.
                
                Int S_k = S_Tree.FormerChildID();
                Int T_k = T_Tree.FormerChildID();
                
//                    print("Coming from "+ToString(S_k)+"-th child of S and "+ToString(T_k)+"-th child of T.");
                
                if( T_k < T_Tree.ChildCount()-1 )
                {
                    S_Tree.ToChild(S_k);
                    T_Tree.ToChild(T_k+1);
                    from_above = true;
                }
                else
                {
                    if( S_k < S_Tree.ChildCount()-1 )
                    {
                        S_Tree.ToChild(S_k+1);
                        T_Tree.ToChild(0);
                        from_above = true;
                    }
                    else
                    {
                        // No further unvisited children. Either move up or break.
                        if( S_Tree.Level() == 0 )
                        {
                            shall_continue = false;
                            from_above = true;
                        }
                        else
                        {
                            S_Tree.ToParent();
                            T_Tree.ToParent();
                            from_above = false;
                        }
                    }
                }
                
            } // if( from_above )
            
        } // while( shall_continue )

        this->total_sum += this->symmetry_factor * sum;
        
        if constexpr ( metric_flag )
        {
            copy_buffer( &ij_block[0], &metric_data[BLOCK_NNZ * k_global], BLOCK_NNZ );
        }
        
        return this->symmetry_factor * sum;
    }
