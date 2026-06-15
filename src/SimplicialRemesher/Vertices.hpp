//############################################################
//##    vertex related
//############################################################

protected:


    bool VertexActiveQ( const Vertex_T v ) const
    {
        return (V_state[v] & VertexActiveMask);
    }

    bool VertexPinnedQ( const Vertex_T v ) const
    {
        return (V_state[v] & VertexPinnedMask);
    }

    bool VertexModifiedQ( const Vertex_T v ) const
    {
        return (V_state[v] & VertexModifiedMask);
    }

    Vertex_T CreateVertex()
    {
        Vertex_T w = vertex_count++;
        
        if( vertex_count >= max_vertex_count)
        {
#ifdef REMESHER_VERBATIM
            print(className()+"::CreateVertex: Reassembling vertex array.");
#endif
            max_vertex_count *= Int(2);
            
            V_coords .template Resize<true>( max_vertex_count, AMB_DIM );
            V_charges.template Resize<true>( max_vertex_count );
            V_state  .template Resize<true>( max_vertex_count );
            
            // TODO: Should not be necessary.
            for( Int v = vertex_count; v < max_vertex_count; ++v )
            {
                V_state[v] = VertexState_T{0};
            }
            
            if( with_data )
            {
                V_data.template Resize<true>( max_vertex_count, V_data.Dim(1) );
            }
        }
        
        V_parent_simplices.push_back( SimplexList_T() );
        
        V_charges[w] = Real(1);
        V_state  [w] = (VertexActiveMask | VertexModifiedMask);
        
//        V_activeQ  [w] = true;
//        V_modifiedQ[w] = true;
        
        return w;
    }
    
    void DeactivateVertex( const Vertex_T v )
    {
#ifdef REMESHER_VERBATIM
        print(className()+"::DeactivateVertex("+ToString(v)+")");
#endif
        V_state[v] = Vertex_T(0);
    }
    
    void MarkVertexAsModified( const Vertex_T v )
    {
#ifdef REMESHER_VERBATIM
        print(className()+"::MarkVertexAsModified("+ToString(v)+")");
#endif
        V_state[v] |= VertexModifiedMask;
    }

    void PinVertex( const Vertex_T v )
    {
    #ifdef REMESHER_VERBATIM
        print(className()+"::PinVertex("+ToString(v)+")");
    #endif
        V_state[v] |= VertexPinnedMask;
    }


    void ComputeSplitVertexPosition( const Vertex_T v_0, const Vertex_T v_1, const Vertex_T w )
    {
        ComputeVertexPosition( v_0, Frac<Real>(1,2), v_1, Frac<Real>(1,2), w );
    }
    
    void ComputeCollapseVertexPosition( const Vertex_T v_0, const Vertex_T v_1, const Vertex_T w )
    {
        constexpr Real W [2][2] = {
            { Frac<Real>(1,2), Real(0)         },
            { Real(1),         Frac<Real>(1,2) }
        };
        
        const bool v_0_pinnedQ = VertexPinnedQ(v_0);
        const bool v_1_pinnedQ = VertexPinnedQ(v_1);
        
        // DEBGUGGING
        if( v_0_pinnedQ && v_1_pinnedQ )
        {
            eprint(ClassName()+"::ComputeCollapseVertexPosition: edge with two pinned vertices.");
        }
        
        const Real weight_0 = W[v_0_pinnedQ][v_1_pinnedQ];
        const Real weight_1 = W[v_1_pinnedQ][v_0_pinnedQ];
        
        ComputeVertexPosition( v_0, weight_0, v_1, weight_1, w );
    }

    void ComputeVertexPosition(
        const Vertex_T v_0, const Real weight_0,
        const Vertex_T v_1, const Real weight_1,
        const Vertex_T w
    )
    {
        for( Int k = 0; k < AMB_DIM; ++k )
        {
            V_coords(w,k) = weight_0 * V_coords(v_0,k) + weight_1 * V_coords(v_1,k);
        }

        if( with_data )
        {
            const Int data_dim = V_data.Dim(1);

            for( Int k = 0; k < data_dim; ++k )
            {
                V_data(w,k) = weight_0 *  V_data(v_0,k) + weight_1 * V_data(v_1,k);
            }
        }
    }

//    void ComputeVertexPosition( const Vertex_T v_0, const Vertex_T v_1, const Vertex_T w )
//    {
//        Quadric_T Q   ( V_quadrics.data(v_0) );
//        Quadric_T Q_1 ( V_quadrics.data(v_1) );
//
//        Q += Q_1;
//
//        Q.Write( V_quadrics.data(w) );
//
//        Tiny::SelfAdjointMatrix<AMB_DIM,Real,Int> A;
//        Tiny::Vector<           AMB_DIM,Real,Int> b;
//
//        for( Int i = 0; i < AMB_DIM; ++i )
//        {
//            b[i]  = - Q[AMB_DIM][i];
//
//            for( Int j = i; j < AMB_DIM; ++j )
//            {
//                A[i][j] = Q[i][j];
//            }
//        }
//
//        Tiny::Matrix<AMB_DIM,AMB_DIM,Real,Int> U;
//        Tiny::Matrix<AMB_DIM,AMB_DIM,Real,Int> UT;
//        Tiny::Vector<AMB_DIM,        Real,Int> lambda;
//        Tiny::Vector<AMB_DIM,        Real,Int> x;
//        Tiny::Vector<AMB_DIM,        Real,Int> y;
//
//        A.Eigensystem( UT, lambda, sqrt_eps, 100 );
//
//        Real lambda_max = Abs(lambda[0]);
//        Real lambda_min = Abs(lambda[0]);
//
//        for( Int i = 1; i < AMB_DIM; ++i )
//        {
//            lambda_min = Min( lambda_min, Abs(lambda[i]) );
//            lambda_max = Max( lambda_max, Abs(lambda[i]) );
//        }
//
//        if( lambda_min > sqrt_eps * lambda_max )
//        {
//            UT.Transpose(U);
//
//            Dot<Overwrite>( U, b, y );
//
//            const Real Lambda_tol = lambda_max * sqrt_eps;
//
//            for( Int i = 0; i < AMB_DIM; ++i )
//            {
//                y[i] = (Abs(lambda[i]) > Lambda_tol) ? y[i] / lambda[i] : Scalar::Zero<Real>;
//            }
//
//            Dot<Overwrite>( UT, y, x );
//
//    //        // Check with
//    //        Dot<Overwrite>( A, x, y );
//    //        y -= b;
//    //
//    //        if( y.Norm() > 0.001 * b.Norm() )
//    //        {
//    //            TOOLS_DUMP(w);
//    //            TOOLS_DUMP(lambda);
//    //            TOOLS_DUMP(U);
//    //            TOOLS_DUMP(y.Norm());
//    //            TOOLS_DUMP(Q);
//    //            TOOLS_DUMP(x);
//    //        }
//
//            x.Write( V_coords.data(w) );
//        }
//        else
//        {
//            for( Int k = 0; k < AMB_DIM; ++k )
//            {
//                V_coords(w,k) = Scalar::Half<Real> * ( V_coords(v_0,k) + V_coords(v_1,k) );
//            }
//
//            if( with_data )
//            {
//                const Int data_dim = V_data.Dim(1);
//
//                for( Int k = 0; k < data_dim; ++k )
//                {
//                    V_data(w,k) = Scalar::Half<Real> * ( V_data(v_0,k) + V_data(v_1,k) );
//                }
//            }
//        }
//    }
    
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
