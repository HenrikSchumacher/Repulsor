
int CheckVertex( const Vertex_T v ) const
{
    if( v < Vertex_T(0) )
    {
        eprint(className()+"::CheckVertex: Vertex "+ToString(v)+" is invalid.");
        return -11;
    }
    
    if( !VertexActiveQ(v) )
    {
#ifdef REMESHER_VERBATIM
        wprint(className()+"::CheckVertex: Vertex "+ToString(v)+" is already deleted.");
#endif
        return -2;
    }
    
    if( VertexModifiedQ(v) )
    {
#ifdef REMESHER_VERBATIM
        wprint(className()+"::CheckVertex: Vertex "+ToString(v)+" is already modified.");
#endif
        return -1;
    }
    
    return 0;
}

int CheckEdge( const Edge_T e ) const
{
    if( e < Edge_T(0) )
    {
//#ifdef REMESHER_VERBATIM
        eprint(className()+"::CheckEdge: edge "+ToString(e)+" is invalid.");
//#endif
        return -11;
    }
    
    if( !E_activeQ[e] )
    {
#ifdef REMESHER_VERBATIM
        wprint(className()+"::CheckEdge: edge "+ToString(e)+" is already deleted.");
#endif
        return -2;
    }
    
    const Vertex_T v_0 = edges(e,0);
    
    if( v_0 < Vertex_T(0) )
    {
        eprint(className()+"::CheckEdge: Vertex "+ToString(v_0)+" is invalid.");
        return -11;
    }
    
    if( !VertexActiveQ(v_0) )
    {
//#ifdef REMESHER_VERBATIM
        eprint(className()+"::CheckEdge: Vertex "+ToString(v_0)+" is already deleted.");
//#endif
        return -2;
    }
    
    if( VertexModifiedQ(v_0) )
    {
#ifdef REMESHER_VERBATIM
        wprint(className()+"::CheckEdge: Vertex "+ToString(v_0)+" is already modified.");
#endif
        return -1;
    }
    
    const Vertex_T v_1 = edges(e,1);
    
    if( v_1 < Vertex_T(0) )
    {
        eprint(className()+"::CheckEdge: Vertex "+ToString(v_1)+" is invalid.");
        return -11;
    }
    
    if( !VertexActiveQ(v_1) )
    {
//#ifdef REMESHER_VERBATIM
        eprint(className()+"::CheckEdge: Vertex "+ToString(v_1)+" is already deleted.");
//#endif
        return -2;
    }
    
    if( VertexModifiedQ(v_1) )
    {
#ifdef REMESHER_VERBATIM
        wprint(className()+"::CheckEdge: Vertex "+ToString(v_1)+" is already modified.");
#endif
        return -1;
    }
    
    if( v_0 == v_1 )
    {
        eprint(className()+"::CheckEdge: edge "+ToString(e)+" is topologically degenerate.");
        return -4;
    }
    
    if( VertexPinnedQ(v_0) && VertexPinnedQ(v_1) )
    {
#ifdef REMESHER_VERBATIM
        wprint(className()+"::CheckEdge: Both vertices "+ToString(v_0) + " and "+ToString(v_1)+" are pinned.");
#endif
        return -1;
    }
    
    return 0;
}


int CheckSimplex( const Simplex_T s ) const
{
    if( s < Simplex_T(0) )
    {
#ifdef REMESHER_VERBATIM
        eprint(className()+"::CheckSimplex: Simplex "+ToString(s)+" is invalid.");
#endif
        return -11;
    }
    
    if( !S_activeQ[s] )
    {
#ifdef REMESHER_VERBATIM
        eprint(className()+"::CheckSimplex: Simplex "+ToString(s)+" is already deleted.");
#endif
        return -2;
    }
    
    bool okay = true;
    
    for( Int i = 0; i < S_vertex_count; ++i )
    {
        for( Int j = i+1; j < S_vertex_count; ++j )
        {
            okay = okay && (simplices[s][i] != simplices[s][j] );
        }
    }
    
    if( !okay )
    {
//#ifdef REMESHER_VERBATIM
        eprint(className()+"::CheckSimplex: Simplex "+ToString(s)+" containes duplicate vertices.");
//#endif
        return -3;
    }
    
    return 0;
}
