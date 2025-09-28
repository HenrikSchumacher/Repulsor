
    Tensor1<MortonCode_T,Int> CoordsToMorton( cptr<Real> X, const Int point_count )
    {
        ComputeBoundingBox( X, point_count );
        
        TOOLS_PTIMER(timer,ClassName()+"::CoordsToMorton");
        
        Tensor1<MortonCode_T,Int> M ( point_count );
        
        ParallelDo(
            [&,this]( const Int thread )
            {
                Axes_T a;
                
                const Int i_begin = JobPointer( point_count, thread_count, thread     );
                const Int i_end   = JobPointer( point_count, thread_count, thread + 1 );
                
                for( Int i = i_begin; i < i_end; ++ i )
                {
                    CoordsToAxes( &X[n*i], a );

                    AxesToMorton( a, M[i] );
                    
                }
            },
            thread_count
        );
        
        return M;
    }
