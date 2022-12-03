//#####################################################################################################
//##    edge related
//#####################################################################################################
        
protected:

    Edge_T CreateEdge( const Vertex_T v_0, const Vertex_T v_1 )
    {
        Pair_T p ( std::minmax(v_0,v_1) );
        
        Edge_T e = edge_count++;
        
        if( e >= max_edge_count )
        {
#ifdef REMESHER_VERBATIM
            print(className()+"::CreateEdge: Reassembling edge array.");
#endif
            max_edge_count *= two;
            
            edges.Resize( max_edge_count, 2 );
            E_active.Resize( max_edge_count );
        }
    
        edge_lookup.insert( {p,e} );
        
        edges(e,0) = p.first;
        edges(e,1) = p.second;
        
        E_active[e] = true;
        
        E_parent_simplices.push_back( SimplexList_T() );
        
#ifdef REMESHER_DEBUG
        print(className()+"::CreateEdge: Created edge with ID = "+ToString(e)+" and vertices {"+ToString(v_0)+","+ToString(v_1)+"}.");
#endif
        return e;
    }
    
    //DONE.
    Edge_T RequireEdge( const Vertex_T v_0, const Vertex_T v_1 )
    {
        Edge_T e = FindEdge(v_0,v_1);
        
        if( e < zero )
        {
            e = CreateEdge(v_0,v_1);
            
#ifdef REMESHER_DEBUG
            print(className()+"RequireEdge: Created edge "+ToString(e)+" with vertices {"+ToString(v_0)+","+ToString(v_1)+"}.");
#endif
        }
        return e;
    }
    
    void DeleteEdge( const Edge_T e )
    {
#ifdef REMESHER_VERBATIM
        print(className()+"::Edge_Delete("+ToString(e)+")");
#endif
        E_active[e] = false;
        LookupErase(e);
    }
    
    //DONE.
    Edge_T FindEdge( const Vertex_T v_0, const Vertex_T v_1 )
    {
        const Pair_T p ( std::minmax(v_0,v_1) );
        return ( edge_lookup.count(p) > 0) ? edge_lookup[p] : -1;
    }
    
    //DONE.
    void LookupErase( const Edge_T e )
    {
        edge_lookup.erase( std::minmax( edges(e,0), edges(e,1) ) );
    }
    
    void LookupErase( const Vertex_T v_0, const Vertex_T v_1 )
    {
        edge_lookup.erase( std::minmax(v_0, v_1) );
    }
    
    //DONE.
    void LookupInsert( const Edge_T e )
    {
        edge_lookup.insert( { std::minmax( edges(e,0), edges(e,1) ), e } );
    }
    
    void LookupInsert( const Edge_T e, const Vertex_T v_0, const Vertex_T v_1 )
    {
        edge_lookup.insert( { std::minmax(v_0, v_1), e } );
    }
    
    //DONE.
    void FindVerticesOpposingEdge( const Edge_T e, VertexList_T & opp_vertices ) const
    {
        // Go through all simplices containing edge e and collect all the simplices' vertices other than the two vertices of edge e.
        
        opp_vertices.Clear();
        
        const Vertex_T v_0 = edges(e,0);
        const Vertex_T v_1 = edges(e,1);
        
        // Going through the simplices to find opposing vertices.
        for( Simplex_T s : E_parent_simplices[e] )
        {
            for( Int i = 0; i < S_vertex_count; ++i )
            {
                const Vertex_T v = simplices(s,i);
                
                if( (v != v_0) && (v != v_1) )
                {
                    opp_vertices.Insert(v);
                }
            }
        }
    }
    
    //DONE.
    Real SquaredEdgeLength( const Edge_T e ) const
    {
        Real L2 = 0.;
        
        Real const * restrict const V = V_coords.data(edges(e,0));
        Real const * restrict const W = V_coords.data(edges(e,1));
        
        for( Int k = 0; k < AMB_DIM; ++k )
        {
            Real delta = V[k] - W[k];
            L2 += delta * delta;
        }
        
        return L2;
    }
