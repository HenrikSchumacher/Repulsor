//######################################################################################################
//##    simplex related
//######################################################################################################

protected:

    //DONE.
    void ComputeSimplexConnectivity( const Simplex_T s )
    {
//            S_active[s] = true;
        
        for( Int i = 0; i < S_vertex_count; ++i )
        {
            Vertex_T v = simplices[s][i];
            
            V_parent_simplices[v].Insert(s);
            
            for( Int j = i+1; j < S_vertex_count; ++j )
            {
                Vertex_T w = simplices[s][j];
                
                Edge_T e = RequireEdge(v,w);
                
                E_parent_simplices[e].Insert(s);
            }
        }
    }
    
    Simplex_T CreateSimplex( const Vertex_T * vertex_list )
    {
        Simplex_T s = simplex_count++;
        
        if( s >= max_simplex_count )
        {
#ifdef REMESHER_VERBATIM
            print(className()+"::CreateSimplex: Reassembling simplex array.");
#endif
            max_simplex_count *= 2;
            
            simplices.template Resize<true>(max_simplex_count,DOM_DIM+1);
            S_active .template Resize<true>(max_simplex_count);
        }
        
        S_active[s] = true;
        
        copy_buffer<S_vertex_count>( vertex_list, simplices.data(s) );

        ComputeSimplexConnectivity(s);

        return s;
    }
    
    Edge_T SimplexFindEdge( const Simplex_T s, const Int k )
    {
        // Find the index of the k-th edge in simplex s.
                
        const Pair_T p { MinMax( simplices[s][Tri_i[k]], simplices[s][Tri_j[k]] ) };
        
        return ( edge_lookup.count(p) > 0 ) ? edge_lookup[p] : -13;
    }
    
    void SimplexComplement(
        const Simplex_T s,
        const Vertex_T v_0,
        const Vertex_T v_1,
        mptr<Vertex_T> vertex_list
    )
    {
        // Fill v_list with the vertices in simplex s that oppose v_0 and v_1.
        
        Int counter = 0;

        for( Int i = 0; i < S_vertex_count; ++i )
        {
            const Vertex_T v = simplices(s,i);
            
            if( (v != v_0) && (v != v_1) )
            {
                vertex_list[counter] = v;
                ++counter;
            }
        }
    }
    
    void DeleteSimplex( const Simplex_T s )
    {
        S_active[s] = false;
        
        for( Int i = 0; i < S_vertex_count; ++i )
        {
            const Vertex_T v = simplices(s,i);
            
            if( v >= zero )
            {
                V_parent_simplices[v].Drop(s);
            }
        }
    }
