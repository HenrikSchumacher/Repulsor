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
            V_charges.Resize( max_vertex_count );
            V_active.Resize( max_vertex_count );
            V_modified.Resize( max_vertex_count );
            
            if( with_data )
            {
                V_data.Resize( max_vertex_count, V_data.Dimension(1) );
            }
        }
        
        V_parent_simplices.push_back( SimplexList_T() );
        
        V_charges [w] = Scalar::One<Real>;
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
//    //            dump(w);
//    //            dump(lambda);
//    //            dump(U);
//    //            dump(y.Norm());
//    //            dump(Q);
//    //            dump(x);
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
//                const Int data_dim = V_data.Dimension(1);
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
