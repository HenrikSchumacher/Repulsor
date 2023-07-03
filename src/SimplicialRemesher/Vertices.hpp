//################################################################################################
//##    vertex related
//################################################################################################

protected:

    Vertex_T CreateVertex()
    {
        Vertex_T w = vertex_count++;
        
        if( vertex_count >= max_vertex_count)
        {
#ifdef REMESHER_VERBATIM
            print(className()+"::CreateVertex: Reassembling vertex array.");
#endif
            max_vertex_count *= two;
            
            V_coords.Resize( max_vertex_count, AMB_DIM );
            V_active.Resize( max_vertex_count );
            V_modified.Resize( max_vertex_count );
            
            if( with_data )
            {
                V_data.Resize( max_vertex_count, V_data.Dimension(1) );
            }
        }
        
        V_parent_simplices.push_back( SimplexList_T() );
        
        V_active  [w] = true;
        V_modified[w] = true;
        
        return w;
    }
    
    void DeactivateVertex( const Vertex_T v )
    {
#ifdef REMESHER_VERBATIM
        print(className()+"::DeactivateVertex("+ToString(v)+")");
#endif
        V_active[v] = false;
    }
    
    void MarkVertexAsModified( const Vertex_T v )
    {
#ifdef REMESHER_VERBATIM
        print(className()+"::MarkVertexAsModified("+ToString(v)+")");
#endif
        V_modified[v] = true;
    }
    
    void ComputeVertexPosition( const Vertex_T v_0, const Vertex_T v_1, const Vertex_T w )
    {
        for( Int k = 0; k < AMB_DIM; ++k )
        {
            V_coords(w,k) = Scalar::Half<Real> * ( V_coords(v_0,k) + V_coords(v_1,k) );
        }
        
        if( with_data )
        {
            const Int data_dim = V_data.Dimension(1);
            
            for( Int k = 0; k < data_dim; ++k )
            {
                V_data(w,k) = Scalar::Half<Real> * ( V_data(v_0,k) + V_data(v_1,k) );
            }
        }
    }
    
    void VertexNeighboringVertices( const Vertex_T v, VertexList_T & neighbors ) const
    {
        neighbors.Clear();
        
        for( Simplex_T s : V_parent_simplices[v] )
        {
            for( Int i = 0; i < S_vertex_count; ++i )
            {
                const Vertex_T w = simplices(s,i);
                
                if( v != w )
                {
                    neighbors.Insert(w);
                }
            }
        }
    }
