
Int CheckVertex( const Vertex_T v ) const
{
    if( v < zero )
    {
        eprint(className()+"::CheckVertex: Vertex "+ToString(v)+" is invalid.");
        return -11;
    }
    
    if( !V_active[v] )
    {
#ifdef REMESHER_VERBATIM
        wprint(className()+"::CheckVertex: Vertex "+ToString(v)+" is already deleted.");
#endif
        return -2;
    }
    
    if( V_modified[v] )
    {
#ifdef REMESHER_VERBATIM
        wprint(className()+"::CheckVertex: Vertex "+ToString(v)+" is already modified.");
#endif
        return -1;
    }
    
    return 0;
}

Int CheckEdge( const Edge_T e ) const
{
    if( e < zero )
    {
//#ifdef REMESHER_VERBATIM
        eprint(className()+"::CheckEdge: edge "+ToString(e)+" is invalid.");
//#endif
        return -11;
    }
    
    if( !E_active[e] )
    {
#ifdef REMESHER_VERBATIM
        wprint(className()+"::CheckEdge: edge "+ToString(e)+" is already deleted.");
#endif
        return -2;
    }
    
    const Vertex_T v_0 = edges(e,0);
    
    if( v_0 < zero )
    {
        eprint(className()+"::CheckVertex: Vertex "+ToString(v_0)+" is invalid.");
        return -11;
    }
    
    if( !V_active[v_0] )
    {
//#ifdef REMESHER_VERBATIM
        eprint(className()+"::CheckVertex: Vertex "+ToString(v_0)+" is already deleted.");
//#endif
        return -2;
    }
    
    if( V_modified[v_0] )
    {
#ifdef REMESHER_VERBATIM
        wprint(className()+"::CheckVertex: Vertex "+ToString(v_0)+" is already modified.");
#endif
        return -1;
    }
    
    const Vertex_T v_1 = edges(e,1);
    
    if( v_1 < zero )
    {
        eprint(className()+"::CheckVertex: Vertex "+ToString(v_1)+" is invalid.");
        return -11;
    }
    
    if( !V_active[v_1] )
    {
//#ifdef REMESHER_VERBATIM
        eprint(className()+"::CheckVertex: Vertex "+ToString(v_1)+" is already deleted.");
//#endif
        return -2;
    }
    
    if( V_modified[v_1] )
    {
#ifdef REMESHER_VERBATIM
        wprint(className()+"::CheckVertex: Vertex "+ToString(v_1)+" is already modified.");
#endif
        return -1;
    }
    
    if( v_0 == v_1 )
    {
        eprint(className()+"::CheckEdge: edge "+ToString(e)+" is topologically degenerate.");
        return -4;
    }
    
    return 0;
}
