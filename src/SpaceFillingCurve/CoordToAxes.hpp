//#######################################################################################
//####                                 CoordToAxes                                  #####
//#######################################################################################

    Tensor1<Axes_T,Int> CoordToAxes( cptr<Real> X, const Int point_count )
    {
        ComputeBoundingBox( X, point_count );
        
        ptic(ClassName()+"::CoordToAxes");
        
        Tensor1<Axes_T,Int> A ( point_count );
        
        ParallelDo(
            [&,this]( const Int thread )
            {
                const Int i_begin = JobPointer<Int>( point_count, thread_count, thread     );
                const Int i_end   = JobPointer<Int>( point_count, thread_count, thread + 1 );
                
                for( Int i = i_begin; i < i_end; ++ i )
                {
                    CoordToAxes( &X[n*i], A[i] );
                }
            },
            thread_count
        );
        
        ptoc(ClassName()+"::CoordToAxes");
        
        return A;
    }


    void CoordToAxes( cptr<Real> x, mref<Axes_T> a )
    {
        for( Int k = 0; k < n; ++k )
        {
            a[k] = static_cast<UInt>( std::floor( ( (x[k]-mid[k]) * scale + Scalar::Half<Real> ) * std::numeric_limits<UInt>::max() ) );
        }
    }
