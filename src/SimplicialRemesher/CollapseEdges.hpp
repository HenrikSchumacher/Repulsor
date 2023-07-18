//################################################################################################
//##    CollapseEdge
//################################################################################################
   
public:

    virtual Int CollapseEdges( cptr<Edge_T> e_list, const Int n ) override
    {
        ptic(className()+"::CollapseEdges");
        
        Int collapse_counter = 0;
        
        for( Int i = 0; i < n; ++i )
        {
            const Int r = CollapseEdge(e_list[i]);
            
            if( r >= zero )
            {
                ++collapse_counter;
            }
        }
        ptoc(className()+"::CollapseEdges");
        
        return collapse_counter;
    }


    virtual Int CollapseEdges( cref<std::vector<Int>> e_list ) override
    {
        return this->CollapseEdges( e_list.data(), int_cast<Int>(e_list.size()) );
    }

protected:

    Int CollapseEdge( const Vertex_T v_0, const Vertex_T v_1 )
    {
        return CollapseEdge( FindEdge(v_0, v_1) );
    }

    Int CollapseEdge( const Edge_T e )
    {
        // Returns the vertex to which the edge has been collapsed -- or an error code (negative number).

        if( CheckEdge(e) < zero )
        {
    #ifdef REMESHER_VERBATIM
            wprint(ClassName()+"::CollapseEdge: edge is not collapsible. Skipping.");
    #endif
            return -11;
        }
        
        const Vertex_T v_0 = edges(e,0);
        const Vertex_T v_1 = edges(e,1);
        
        const Int val_0 = V_parent_simplices[v_0].Size();
        
        if constexpr ( DOM_DIM >=2 )
        {
            //TODO: This bound should be adapted for boundary vertices.
            if( val_0 <= 3 )
            {
    #ifdef REMESHER_VERBATIM
                wprint(className()+"::CollapseEdge: Vertex "+ToString(v_0)+" has simplex valence 3 or less. Skipping.");
    #endif
                return -3;
            }
        }
        
        const Int val_1 = V_parent_simplices[v_1].Size();
        
        if constexpr ( DOM_DIM >=2 )
        {
            //TODO: This bound should be adapted for boundary vertices.
            if( val_1 <= 3 )
            {
    #ifdef REMESHER_VERBATIM
                wprint(className()+"::CollapseEdge: Vertex "+ToString(v_1)+" has simplex valence 3 or less. Skipping.");
    #endif
                return -3;
            }
        }
        
        const Int expected_valence = val_0 + val_1 - E_parent_simplices[e].Size() - 2;
        
        if( expected_valence > V_max_simplex_valence )
        {
    #ifdef REMESHER_VERBATIM
            wprint(className()+"::CollapseEdge: Collapse would result in vertex of edge valence > "+Tools::ToString(V_max_simplex_valence)+". Skipping.");
    #endif
            return -3;
        }
                    
        FindVerticesOpposingEdge( e, E_opp_vertices );
        
        if constexpr( DOM_DIM >= 2 )
        {
            auto found = std::find_if(
                E_opp_vertices.begin(),
                E_opp_vertices.end(),
                [&](const Vertex_T w)
                {
                    return ( V_parent_simplices[w].Size() <= 3 );
                }
            );

            if( found < E_opp_vertices.end() )
            {
    #ifdef REMESHER_VERBATIM
                wprint(className()+"::CollapseEdge: Collapsing "+ToString(e)+" would reduce simplex valence of opposing vertex below 3. Skipping.");
    #endif
                return -4;
            }
            
            VertexNeighboringVertices( v_0, V_0_neighbors );
            VertexNeighboringVertices( v_1, V_1_neighbors );
            
            Intersection( V_0_neighbors, V_1_neighbors, E_neighbors );
            
            if( E_neighbors != E_opp_vertices )
            {
    #ifdef REMESHER_VERBATIM
                wprint(className()+"::CollapseEdge: Neighborhood of edge would become nonmanifold after collapse. Skipping.");
    #endif
                return -8;
            }
        }

        // Edge e is collapsible. Start with deleting it.
        DeleteEdge(e);
        MarkVertexAsModified(v_0);
        DeactivateVertex(v_1);

        ComputeVertexPosition(v_0,v_1,v_0);

        // Going through the simplices to delete.
        for( Simplex_T s : E_parent_simplices[e] )
        {
            SimplexComplement( s, v_0, v_1, &opp_buffer[0] );
            
    //        DeleteSimplex(s);
            S_active[s] = false;
            
            V_parent_simplices[v_0].Drop(s);

            // We have to unbind all edges of s, too.
            // We don't have to do that for e as it is already deleted.
            // So there are two kinds of edges to consider:
            //      (i)  Those edges of s that contain either v_0 or v_1 (and not both of them!).
            //      (ii) Those edges of s that contain neither v_0 nor v_1.

            // Case (i):
            // Should not be executed for DOM_DIM <=1. So we ask the compiler to remove this.
            
            if constexpr ( DOM_DIM >= 2 )
            {
                for( Int i = 0; i < S_vertex_count-2; ++i )
                {
                    const Vertex_T w = opp_buffer[i];
                    
                    if( CheckVertex(w) < -one )
                    {
    #ifdef REMESHER_VERBATIM
                        eprint(className()+"::CollapseEdge: Vertex w = "+ToString(w)+" is not present. Skipping it.");
    #endif
                        continue;
                    }
                    
                    MarkVertexAsModified(w);
                    
                    // Unbind s from w.
                    V_parent_simplices[w].Drop(s);
                    
                    Edge_T e_0 = FindEdge(v_0,w);
                    Edge_T e_1 = FindEdge(v_1,w);

                    if( e_0 < zero )
                    {
    #ifdef REMESHER_VERBATIM
                        eprint(className()+"::CollapseEdge: e_0 == {"+ToString(v_0)+","+ToString(w)+"} could not be found in lookup table.");
    #endif
                        return -100;
                    }
                    
                    if( e_1 < zero )
                    {
    #ifdef REMESHER_VERBATIM
                        eprint(className()+"::CollapseEdge: e_1 == {"+ToString(v_1)+","+ToString(w)+"} could not be found in lookup table.");
    #endif
                        return -200;
                    }
                    
                    // Unbind s from e_0.
                    E_parent_simplices[e_0].Drop(s);
                    
                    // We do not have to unbind s from e_1, because we delete e_1 anyways.
                    DeleteEdge(e_1);
                    
                    // Move all the undeleted parent simplices of e_1 to e_0.
                    for( Simplex_T t : E_parent_simplices[e_1] )
                    {
                        if( S_active[t] )
                        {
                            // TODO: In fact, it should be possible to do this without checking for duplicates:
                            E_parent_simplices[e_0].Insert(t);
                        }
                    }
                }
            }

            // TODO: Test this!
            // Case (ii)
            // Unbind all edges e of s that contain neither v_0 nor v_1.
            // Since v_0 and v_1 are already contained in s, S_vertex_count>=4 has to hold for such an edge to exist.
            if constexpr( S_vertex_count >= 3 )
            {
                for( Int i = 0; i < S_vertex_count-2; ++i )
                {
                    for( Int j = i+1; j < S_vertex_count-2; ++j )
                    {
                        Edge_T f = FindEdge(opp_buffer[i],opp_buffer[j]);

                        if( f < zero )
                        {
    #ifdef REMESHER_VERBATIM
                            eprint(className()+"::CollapseEdge: edge { "+ToString(opp_buffer[i]) +", " + ToString(opp_buffer[j])+" } could not be found. (A)");
    #endif
                        }
                        else
                        {
                            E_parent_simplices[f].Drop(s);
                        }
                    }
                }
            }
            
        } // for( Simplex_T s : Edge_ParentSimplices(e) )
        
        // We cycle over the undeleted parent simplices of v_1.
        for( Simplex_T s : V_parent_simplices[v_1] )
        {
            if( S_active[s] )
            {
                // The undeleted parent s of v_1 is from now on also parent of v_0.
                // TODO: It should not be necessary to check for duplicates here as every undeleted parent simplex of v_1 should be visited only once.
                V_parent_simplices[v_0].Insert(s);
                
                // We also have to tell s that v_0 is its child from now on.
                // Moreover, we have to tell the undeleted edges f of s that contain v_1 that v_0 is now their mate.
                // Aaand we also have to correct edge_lookup.
                // We can do all of this by going through the vertices of s:
                
                for( Int i = 0; i < S_vertex_count; ++i )
                {
                    Vertex_T w = simplices(s,i);
                    
                    // Replace v_1 by v_0.
                    if( w == v_1 )
                    {
                        simplices(s,i) = v_0;
                    }
                    else
                    {
                        // We have w != v_1, so {v_1,w} is a proper edge of s.
                        // Finding the edge {v_1,w}.
                        Edge_T f = FindEdge(v_1,w);

                        if( f < zero )
                        {
    //                        eprint(className()+"::CollapseEdge: edge { "+ToString(v_1) +", " + ToString(w)+" } could not be found. (B)");
                            continue;
                        }
                        else
                        {
                            if( E_active[f] )
                            {
                                LookupErase(f);
                                
                                if( edges(f,0) == v_1 )
                                {
                                    edges(f,0) = v_0;
                                }
                                if( edges(f,1) == v_1 )
                                {
                                    edges(f,1) = v_0;
                                }
                                
                                LookupInsert(f);
                            }
                        }
                    }
                }
            }

        } // for( Simplex_T s : VertexParentSimplices(v_1) )
        
        compressed = false;
        
        return v_0;
        
    } // CollapseEdge
