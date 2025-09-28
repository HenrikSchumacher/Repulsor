//############################################################
//##    Unify edge lengths
//############################################################

public:
        
    virtual bool UnifyEdgeLengths(
        const Real lower_bound,
        const Real upper_bound,
        const Int  max_iter = 100
    ) override
    {
        TOOLS_PTIMER(timer,className()+"::UnifyEdgeLengths");

        if( lower_bound > upper_bound  )
        {
            eprint(className()+"::UnifyEdgeLengths: lower_bound > upper_bound. Aborting");
            return 0;
        }

        const Real split_threshold_2    = upper_bound * upper_bound;
        const Real collapse_threshold_2 = lower_bound * lower_bound;

        PairAggregator<Int,Real,Int> splits    (edge_count);
        PairAggregator<Int,Real,Int> collapses (edge_count);
        
//        TwoArrayQuickSort<Real,Int,Int> quick_sort;
        
        
        TwoArraySort<Real,Int,Int,VarSize,std::less   <Real>> sort;
        TwoArraySort<Real,Int,Int,VarSize,std::greater<Real>> reverse_sort;
        
        
        Int total_split_count    = 0;
        Int       split_count    = 1;
        Int total_collapse_count = 0;
        Int       collapse_count = 1;
        Int iter                 = 0;
        
        while( (split_count > zero || collapse_count > zero) && (iter < max_iter) )
        {
            splits.Clear();
            collapses.Clear();
            
            for( Int e = 0; e < edge_count; ++e )
            {
                // Unblock all vertices.
                std::fill( &V_modified[0], &V_modified[vertex_count], false );
                
                if( !E_active[e] )
                {
    #ifdef REMESHER_VERBATIM
                    wprint(className()+"::UnifyEdgeLengths: Skipping edge "+ToString(e)+" because it is inactive.");
    #endif
                    continue;
                }

                const Real L2 = SquaredEdgeLength(e);

                if( L2 > split_threshold_2 )
                {
                    splits.Push(e,L2);
                }
                else if( L2 < collapse_threshold_2 )
                {
                    collapses.Push(e,L2);
                }

            }
            
            // Order such that shortest edges are collapsed first.
            sort( collapses.data_1(), collapses.data_0(), collapses.Size() );

//            print(className()+"::UnifyEdgeLengths: iteration "+ToString(iter)+":");
            
            collapse_count = CollapseEdges( collapses.data_0(), collapses.Size() );
            total_collapse_count += collapse_count;
//            valprint("  collapse_count",collapse_count);
            
            // Order such that longest edges are split first.
            reverse_sort( splits.data_1(), splits.data_0(), splits.Size() );

            split_count = SplitEdges( splits.data_0(), splits.Size() );
            total_split_count += split_count;
//            valprint("  split_count   ",split_count);
            
            ++iter;
        }
        
        if( collapses.Size() > Int(0) )
        {
            wprint(className()+"::UnifyEdgeLengths: "+ToString(collapses.Size())+" short edges could not be collapsed.");
        }
                   
        if( splits.Size() > Int(0) )
        {
           wprint(className()+"::UnifyEdgeLengths: "+ToString(splits.Size())+" long  edges could not be split.");
        }

        return (total_collapse_count > zero) || (total_split_count > zero);
        
    } // UnifyEdgeLengths
