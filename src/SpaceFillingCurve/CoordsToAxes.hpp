
    void CoordsToAxes( cptr<Real> x, mref<Axes_T> a )
    {
        for( Int k = 0; k < n; ++k )
        {
            a[k] = static_cast<UInt>(
                std::floor( ( (x[k]-mid[k]) * scale + Scalar::Half<Real> ) * Scalar::Max<UInt> )
            );
        }
    }

    Tensor1<Axes_T,Int> CoordsToAxes( cptr<Real> X, const Int point_count )
    {
        ComputeBoundingBox( X, point_count );
        
        TOOLS_PTIMER(timer,ClassName()+"::CoordsToAxes");
        
        Tensor1<Axes_T,Int> A ( point_count );
        
        ParallelDo(
            [&,this]( const Int thread )
            {
                const Int i_begin = JobPointer( point_count, thread_count, thread     );
                const Int i_end   = JobPointer( point_count, thread_count, thread + 1 );
                
                for( Int i = i_begin; i < i_end; ++ i )
                {
                    CoordsToAxes( &X[n*i], A[i] );
                }
            },
            thread_count
        );
        
        return A;
    }
