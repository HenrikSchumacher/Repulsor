public:

    force_inline Real Compute( const LInt k_global )
    {
        Real sum = static_cast<Real>(0);
        
        bool from_above = true;
        bool shall_continue = true;
        
        S_Tree.ToChild(0);
        T_Tree.ToChild(0);
    
        while( shall_continue )
        {
//            logprint("level { "+ToString(S_Tree.Level())+","+ToString(T_Tree.Level())+" }");
            
            if( from_above )
            {
//                logprint("from_above");
                if( S_Tree.Level() >= this->max_refinement )
                {
                    // If at lowest level and inadmissable then we just compute the energy and move up.
                    this->max_level_reached = std::max( this->max_level_reached,this->max_refinement );
                    ++this->evaluations;
                    sum += compute();
                    S_Tree.ToParent();
                    T_Tree.ToParent();
                    from_above = false;
                }
                else
                {
//                    logprint("not max_level");
                    // If not at lowest level, then we have to check for admissability.
                    const bool admissableQ = this->gjk.MultipoleAcceptanceCriterion(
                        S_Tree.SimplexPrototype(),
                        T_Tree.SimplexPrototype(),
                        this->theta2
                    );
                    
                    if( admissableQ )
                    {
//                        logprint("admissable");
                        // We compute energy, go to parent, and prepare the next child of the parent.
                        this->max_level_reached = std::max( this->max_level_reached, S_Tree.Level() );
                        ++this->evaluations;
                        sum += compute();
                        S_Tree.ToParent();
                        T_Tree.ToParent();
                        from_above = false;
                    }
                    else
                    {
//                        logprint("inadmisable");
                        // If inadmissable, we go a level deeper.

                        S_Tree.ToChild(0);
                        T_Tree.ToChild(0);
                        from_above = true;
                    }
                }
            }
            else
            {
//                logprint("from_below");
                // If we come from below, we have to find the next pair of simplices to visit.
                
                const typename Base_T::S_Tree_T::Child_T S_k = S_Tree.FormerChildID();
                const typename Base_T::T_Tree_T::Child_T T_k = T_Tree.FormerChildID();
                
//                logprint("Coming from "+ToString(S_k)+"-th child of S and "+ToString(T_k)+"-th child of T.");
                
                if( T_k + 1 < T_Tree.ChildCount() )
                {
                    S_Tree.ToChild(S_k);
                    T_Tree.ToChild(T_k+1);
                    from_above = true;
                }
                else
                {
                    if( S_k + 1 < S_Tree.ChildCount() )
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
            copy_buffer<BLOCK_NNZ>( &ij_block[0], &metric_data[BLOCK_NNZ * k_global] );
        }
        
        return this->symmetry_factor * sum;
    }


    void PrintReport( Int thread, float time ) const
    {
        logprint("Thread "+ToString(thread)+": "
            + "\n\t kernel:        \t" + this->ClassName()
            + "\n\t time elapsed:  \t" + ToString(time)
            + "\n\t evaluations:   \t" + ToString(this->evaluations)
            + "\n\t deepest level: \t" + ToString(this->max_level_reached)
            + "\n"
        );
        
        if( this->max_level_reached >= this->max_refinement )
        {
            wprint(this->ClassName()+" on thread "+ToString(thread)+" reached the maximal allowed refinement level "+ToString(this->max_refinement)+".");
        }
    }
