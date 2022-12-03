//######################################################################################################
//##    Unify edge lengths
//######################################################################################################

public:
        
    virtual bool UnifyEdgeLengths(
        const Real collapse_threshold,
        const Real split_threshold,
        const Int  max_iter = 100
    ) override
    {
        ptic(className()+"::UnifyEdgeLengths");

        if( split_threshold < collapse_threshold )
        {
            eprint(className()+"::UnifyEdgeLengths: split_threshold < collapse_threshold. Aborting");
            ptoc(className()+"::UnifyEdgeLengths");
            return 0;
        }

        const Real split_threshold_2    = split_threshold * split_threshold;
        const Real collapse_threshold_2 = collapse_threshold * collapse_threshold;

        PairAggregator<Int,Real,Int> splits    (edge_count);
        PairAggregator<Int,Real,Int> collapses (edge_count);
        
        TwoArrayQuickSort<Real,Int,Int> quick_sort;
        
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
            quick_sort( collapses.Get_1().data(), collapses.Get_0().data(), collapses.Size(), false);

            print(className()+"::UnifyEdgeLengths: iteration "+ToString(iter)+":");
            
            collapse_count = CollapseEdges( collapses.Get_0().data(), collapses.Size() );
            total_collapse_count += collapse_count;
            valprint("  collapse_count",collapse_count);
            
            // Order such that longest edges are split first.
            quick_sort( splits.Get_1().data(), splits.Get_0().data(), splits.Size(), true );

            split_count = SplitEdges( splits.Get_0().data(), splits.Size() );
            total_split_count += split_count;
            valprint("  split_count   ",split_count);
            
            ++iter;
        }
        
        if( collapses.Size() > 0 )
        {
            wprint(className()+"::UnifyEdgeLengths: "+ToString(collapses.Size())+" short edges could not be collapsed.");
        }
                   
       if( splits.Size() > 0 )
       {
           wprint(className()+"::UnifyEdgeLengths: "+ToString(splits.Size())+" long  edges could not be split.");
       }

        ptoc(className()+"::UnifyEdgeLengths");

        return (total_collapse_count > zero) || (total_split_count > zero);
        
    } // UnifyEdgeLengths
