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
        
        TwoArraySort<Real,Int,Int,VarSize,std::less   <Real>> sort;
        TwoArraySort<Real,Int,Int,VarSize,std::greater<Real>> reverse_sort;
        
        
        Int total_split_count    = 0;
        Int       split_count    = 1;
        Int total_collapse_count = 0;
        Int       collapse_count = 1;
        Int iter                 = 0;
        
        while( ((split_count > Int(0)) || (collapse_count > Int(0))) && (iter < max_iter) )
        {
            splits.Clear();
            collapses.Clear();

            // Unblock all vertices.
            for( Int v = 0; v < vertex_count; ++v )
            {
                V_state[v] &= (~VertexModifiedMask);
            }
            
            for( Int e = 0; e < edge_count; ++e )
            {
                if( !E_activeQ[e] )
                {
    #ifdef REMESHER_VERBATIM
                    wprint(className()+"::UnifyEdgeLengths: Skipping edge "+ToString(e)+" because it is inactive.");
    #endif
                    continue;
                }
                
                const Vertex_T v_0 = edges(e,0);
                const Vertex_T v_1 = edges(e,1);
                
                // We don't want to mess with edges whose vertices are both pinned.
                if( VertexPinnedQ(v_0) && VertexPinnedQ(v_1) ) { continue; }

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

            collapse_count = CollapseEdges( collapses.data_0(), collapses.Size() );
            total_collapse_count += collapse_count;
            
            // Order such that longest edges are split first.
            reverse_sort( splits.data_1(), splits.data_0(), splits.Size() );
            split_count = SplitEdges( splits.data_0(), splits.Size() );
            total_split_count += split_count;
            
            ++iter;
        }
        
        Int non_collapsed_count = collapses.Size() - collapse_count;
        Int non_split_count     = splits.Size() - split_count;
        
        if( non_collapsed_count > Int(0) )
        {
            wprint(className()+"::UnifyEdgeLengths: "+ToString(non_collapsed_count)+" short edges could not be collapsed.");
        }
                   
        if( non_split_count > Int(0) )
        {
           wprint(className()+"::UnifyEdgeLengths: "+ToString(non_split_count)+" long  edges could not be split.");
        }

        return (total_collapse_count > Int(0)) || (total_split_count > Int(0));
        
    } // UnifyEdgeLengths
