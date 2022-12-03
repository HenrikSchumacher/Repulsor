//################################################################################################
//##    FlipEdge
//################################################################################################

virtual Int DelaunayFlip( const Int max_iter = 128 ) override
{
    ptic(className()+"::DelaunayFlip");
    
    if constexpr ( DOM_DIM != 2 )
    {
        wprint(className()+"::DelaunayFlip: only implemented for domain dimension 2. Skipping.");
        ptoc(className()+"::DelaunayFlip");
        return 0;
    }
    
    // For each edge in e_list perform a check for the Delaunay condition. If successful, perform the flips and return the number of flipped edges.
    Int total_flip_count = 0;
    Int flip_counter = 1;
    Int iter = 0;
    
    while( flip_counter > 0 && iter < max_iter )
    {
        ++iter;
        
        flip_counter = 0;
        
        for( Edge_T e = 0; e < edge_count; ++e )
        {
            const Int r = FlipEdge( e, true );
            
            if( r >= zero )
            {
                ++flip_counter;
            }
            else
            {
#ifdef REMESHER_VERBATIM
                wprint(className()+"::DelaunayFlip: CollapseEdge failed to collapse edge "+ToString(e)+".");
#endif
            }
        }
        
        total_flip_count += flip_counter;
    }
    
    if( iter < max_iter )
    {
        print(className()+"::DelaunayFlip: Success after "+ToString(total_flip_count)+" flips  ("+ToString(iter)+" iterations).");
    }
    else
    {
        wprint(className()+"::DelaunayFlip: Abort after "+ToString(total_flip_count)+" flips  ("+ToString(iter)+" iterations).");
    }
    
    ptoc(className()+"::DelaunayFlip");
    
    return total_flip_count;
}

virtual Int FlipEdges(
    const Edge_T * restrict const e_list, const Int n, const bool check_Delaunay = false
) override
{
    if constexpr ( DOM_DIM != 2 )
    {
        wprint(className()+"::FlipEdges: only implemented for domain dimension 2. Skipping.");
        return 0;
    }
    
    // If check_Delaunay == true, this performs a check for the Delaunay condition for each edge in e_list. If successful, it performs the flips and returns the number of flipped edges.
    Int flip_counter = 0;
    
    for( Int i = 0; i < n; ++i )
    {
        const Int r = FlipEdge( e_list[i], check_Delaunay );
        
        if( r >= zero )
        {
            ++flip_counter;
        }
        else
        {
#ifdef REMESHER_VERBATIM
            wprint(className()+"::FlipEdges: CollapseEdge failed to collapse edge "+ToString(e_list[i])+".");
#endif
        }
    }
    
    return flip_counter;
}



Int FlipEdge( const Vertex_T v_0, const Vertex_T v_1, const bool check_Delaunay = false )
{
    return FlipEdge( FindEdge(v_0, v_1), check_Delaunay );
}

Int FlipEdge( const Edge_T e, const bool check_Delaunay = false )
{
    if constexpr ( DOM_DIM != 2 )
    {
        wprint(ClassName()+"::FlipEdge: only implemented for domain dimension 2. Skipping.");
        return -1;
    }
    
    if( e < 0 )
    {
#ifdef REMESHER_VERBATIM
        eprint(ClassName()+"::FlipEdge: edge is invalid. Skipping.");
#endif
        return -1;
    }
    
    if( !E_active[e] )
    {
#ifdef REMESHER_VERBATIM
        wprint(ClassName()+"::FlipEdge: edge is inactive. Skipping.");
#endif
        return -1;
    }
    
    if( E_parent_simplices[e].Size() != two )
    {
#ifdef REMESHER_VERBATIM
        wprint(ClassName()+"::FlipEdge: edge "+ToString(e)+" cannot be fliped since it has simplex valences != 2. Skipping.");
#endif
        return -1;
    }
    
    const Simplex_T s_0 = E_parent_simplices[e][0];
    const Simplex_T s_1 = E_parent_simplices[e][1];
    
    const Vertex_T v_0 = edges(e,0);
    const Vertex_T v_1 = edges(e,1);
    
    if( v_0 < 0 )
    {
//#ifdef REMESHER_VERBATIM
        wprint(className()+"::FlipEdge: Vertex "+ToString(v_0)+" of edge "+ToString(e)+" is invalid. Skipping.");
//#endif
        return -11;
    }
    
    if( v_1 < 0 )
    {
//#ifdef REMESHER_VERBATIM
        wprint(className()+"::FlipEdge: Vertex "+ToString(v_1)+" of edge "+ToString(e)+" is invalid. Skipping.");
//#endif
        return -11;
    }
    
    if( !V_active[v_0] )
    {
//#ifdef REMESHER_VERBATIM
        wprint(className()+"::FlipEdge: Vertex "+ToString(v_0)+" of edge "+ToString(e)+" is deleted. Skipping.");
//#endif
        return -2;
    }
    
    if( !V_active[v_1] )
    {
//#ifdef REMESHER_VERBATIM
        wprint(className()+"::FlipEdge: Vertex "+ToString(v_1)+" of edge "+ToString(e)+" is deleted. Skipping.");
//#endif
        return -2;
    }
    
    
    Int p_0 = -1; // position of vertex v_1 in simplex s_0 before flip
    
    Vertex_T w_0 = simplices(s_0,0);
    
    if( w_0 == v_1 )
    {
        p_0 = 0;
        w_0 = simplices(s_0,1);
        
        if ( w_0 == v_0 )
        {
            w_0 = simplices(s_0,2);
        }
    }
    else if ( w_0 == v_0 )
    {
        w_0 = simplices(s_0,1);
        
        if( w_0 == v_1 )
        {
            p_0 = 1;
            w_0 = simplices(s_0,2);
        }
        else
        {
            p_0 = 2;
        }
    }
    else
    {
        p_0 = ( simplices(s_0,1) == v_1 ) ? 1 : 2;
    }
    
    
    Int p_1 = -1; // position of vertex v_1 in simplex s_0 before flip
    
    Vertex_T w_1 = simplices(s_1,0);
    
    if( w_1 == v_0 )
    {
        p_1 = 0;
        w_1 = simplices(s_1,1);
        
        if( w_1 == v_1 )
        {
            w_1 = simplices(s_1,2);
        }
    }
    else if( w_1 == v_1 )
    {
        w_1 = simplices(s_1,1);
        
        if( w_1 == v_0 )
        {
            p_1 = 1;
            w_1 = simplices(s_1,2);
        }
        else
        {
            p_1 = 2;
        }
    }
    else
    {
        p_1 = ( simplices(s_1,1) == v_0 ) ? 1 : 2;
    }

    /*           before                           after
    //
    //             v_0                             v_0
    //              +                               +
    //             /|\                             / \
    //            / | \                           /   \
    //           /  |  \                         / s_0 \
    //          /   |   \                       /       \
    //         /    |    \                     /    e    \
    //    w_0 + s_0 | s_1 + w_1           w_0 +-----------+ w_1
    //         \    |    /                     \         /
    //          \   |e  /                       \  s_1  /
    //           \  |  /                         \     /
    //            \ | /                           \   /
    //             \|/                             \ /
    //              +                               +
    //             v_1                             v_1
    */

    
    bool do_flip = false;
    
    if( check_Delaunay )
    {
        // Flip only if the sum of angles opposing e is greater than 180 degree.
        Vector_T a;
        Vector_T b;
        Vector_T c;
        Vector_T d;

        for( Int k = 0; k < AMB_DIM; ++k )
        {
            a[k] = V_coords(v_0,k) - V_coords(w_0,k);
            b[k] = V_coords(v_1,k) - V_coords(w_0,k);
            
            c[k] = V_coords(v_0,k) - V_coords(w_1,k);
            d[k] = V_coords(v_1,k) - V_coords(w_1,k);
        }
        
        Real alpha = Angle(a,b);
        Real beta  = Angle(c,d);
        
        do_flip = ( (alpha + beta) > static_cast<Real>(M_PI) );
    }
    else
    {
        do_flip = true;
    }
    
    
    if( do_flip )
    {
        LookupErase(e);
        
        edges(e,0) = std::min(w_0,w_1);
        edges(e,1) = std::max(w_0,w_1);
        
        // replace v_0 by w_0
        simplices(s_1,p_1) = w_0;
        V_parent_simplices[v_0].Drop(s_1);
        V_parent_simplices[w_0].Insert(s_1);
        
        // replace v_1 by w_1
        simplices(s_0,p_0) = w_1;
        V_parent_simplices[v_1].Drop(s_0);
        V_parent_simplices[w_1].Insert(s_0);
        
        LookupInsert(e);
        

        // Taking care of the hinge's boundary edges.
        {
            Edge_T f = FindEdge(w_0,v_1);
            if( f < 0 )
            {
                eprint(className()+"::FlipEdge: edge {"+ToString(w_0)+", "+ToString(v_1)+"} not found. Aborting.");
                return -222;
            }
            E_parent_simplices[f].Drop(s_0);
            E_parent_simplices[f].Insert(s_1);
        }
        {
            Edge_T f = FindEdge(w_1,v_0);
            if( f < 0 )
            {
                eprint(className()+"::FlipEdge: edge {"+ToString(w_1)+", "+ToString(v_0)+"} not found. Aborting.");
                return -222;
            }
            E_parent_simplices[f].Drop(s_1);
            E_parent_simplices[f].Insert(s_0);
        }
        
        return 0;
    }
    else
    {
        return -1;
    }
    
    
    
} // FlipEdge

