//################################################################################################
//##    SplitEdge
//################################################################################################

public:

    virtual Int SplitEdges( ptr<Edge_T> e_list, const Int n ) override
    {
        ptic(className()+"::SplitEdges");
        
        Int split_counter = 0;

        for( Int i = 0; i < n; ++i )
        {
            const Int r = SplitEdge(e_list[i]);
            
            if( r >= zero )
            {
                ++split_counter;
            }
        }
        
        ptoc(className()+"::SplitEdges");
        return split_counter;
    }

virtual Int SplitEdges( const std::vector<Int> & e_list ) override
{
    return this->SplitEdges( e_list.data(), int_cast<Int>(e_list.size()) );
}

protected:

    Int SplitEdge( const Vertex_T v_0, const Vertex_T v_1 )
    {
        return SplitEdge( FindEdge(v_0,v_1) );
    }

    Int SplitEdge( const Edge_T e_0 )
    {
        // Returns the index of newly created vertex -- or some error code (negative number).
        
        if( CheckEdge(e_0) < zero )
        {
    #ifdef REMESHER_VERBATIM
            wprint(ClassName()+"::SplitEdge: edge is not splittable. Skipping.");
    #endif
            return -11;
        }
        
        if constexpr( DOM_DIM >= 1 )
        {
            FindVerticesOpposingEdge( e_0, E_opp_vertices );
            
            for( Vertex_T v : E_opp_vertices )
            {
                if( CheckVertex(v) < zero )
                {
    #ifdef REMESHER_VERBATIM
                    wprint(className()+"::SplitEdge: Vertex "+ToString(v)+" opposing edge "+ToString(e_0)+" is blocked. Skipping.");
    #endif
                    return -22;
                }
                
                if constexpr ( DOM_DIM >= 2 )
                {
                    if( V_parent_simplices[v].Size()+1 > V_max_simplex_valence )
                    {
    #ifdef REMESHER_VERBATIM
                        wprint(className()+"::SplitEdge: Simplex valence of vertex "+Tools::ToString(v)+" opposing edge "+Tools::ToString(e_0)+" would grow to "+Tools::ToString(V_max_simplex_valence)+" or higher. Skipping.");
    #endif
                        return -33;
                    }
                }
            }
            
            const Vertex_T v_0 = edges(e_0,0);
            const Vertex_T v_1 = edges(e_0,1);
            const Vertex_T   w = CreateVertex();

            ComputeVertexPosition(v_0,v_1,w);
            
            const Edge_T e_1 = CreateEdge(w,v_1);

            for( Simplex_T s_0 : E_parent_simplices[e_0] )
            {
                SimplexComplement( s_0, v_0, v_1, &simplex_buffer[0] );

                // Tie s_0 and w to each other.
                for( Int i = 0; i < S_vertex_count; ++ i )
                {
                    if( simplices(s_0,i) == v_1 )
                    {
                        simplices(s_0,i) = w;
                    }
                }
                
                V_parent_simplices[w].Insert(s_0);
                
                simplex_buffer[S_vertex_count-2] = w;
                simplex_buffer[S_vertex_count-1] = v_1;

                const Simplex_T s_1 = CreateSimplex( &simplex_buffer[0] );

                V_parent_simplices[w].Insert(s_1);
                
                V_parent_simplices[v_1].Drop(s_0);
                V_parent_simplices[v_1].Insert(s_1);
                
                E_parent_simplices[e_1].Insert(s_1);

                if constexpr ( DOM_DIM > 1 )
                {
                    for( Int i = 0; i < S_vertex_count-2; ++i )
                    {
                        const Vertex_T u_i = simplex_buffer[i];
                        
                        V_parent_simplices[u_i].Insert(s_1);
                        
                        // Create edges containing w but not v_1.
                        // These new edges must be tied to both s_0 and s_1.
                        const Edge_T f = RequireEdge(w,u_i);
                        
                        E_parent_simplices[f].Insert(s_0);
                        E_parent_simplices[f].Insert(s_1);
                        
                        // Create/find edges containing v_1 but not w.
                        // These new edges need to be tied only to s_1.
                        const Edge_T g = RequireEdge(v_1,u_i);
                        
                        E_parent_simplices[g].Drop(s_0);
                        E_parent_simplices[g].Insert(s_1);
                        
                        if constexpr ( DOM_DIM > 2)
                        {
                            for( Int j = i+1; j < S_vertex_count-2; ++j )
                            {
                                Vertex_T u_j = simplex_buffer[j];
                                // Find edges containing neither v_1 nor w.
                                // These new edges must be tied to s_1. (They should already be tried to s_0.)
                                Edge_T h = FindEdge(u_i,u_j);
                                
                                E_parent_simplices[h].Insert(s_1);
                            }
                        }
                        
                    } // for( Int i = 0; i < simplex_vertex_count-2; ++i )
                    
                } // if( DOM_DIM > 1 )
            }
            
            // Redefine old edge. Caution: Must be done _after_ all calls to OppositeVertices.
            LookupErase(e_0);
            
            if( v_0 <= w )
            {
                edges(e_0,1) = w;
            }
            else
            {
                edges(e_0,0) = w;
                edges(e_0,1) = v_0;
            }
            
            LookupInsert(e_0);
            
            MarkVertexAsModified(w);
            MarkVertexAsModified(v_0);
            MarkVertexAsModified(v_1);
            
            compressed = false;
            
            return w;
        }
        else
        {
            return -1001;
        }
        
    } // SplitEdge
