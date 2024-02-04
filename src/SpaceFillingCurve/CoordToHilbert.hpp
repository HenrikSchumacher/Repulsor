//#######################################################################################
//####                                CoordToHilbert                                #####
//#######################################################################################

    Tensor1<Morton_T,Int> CoordToHilbert( cptr<Real> X, const Int point_count )
    {
        ComputeBoundingBox( X, point_count );
        
        ptic(ClassName()+"::CoordToHilbert");
        
        Tensor1<Morton_T,Int> M ( point_count );
        
        ParallelDo(
            [&,this]( const Int thread )
            {
                Axes_T a;
                
                const Int i_begin = JobPointer( point_count, thread_count, thread     );
                const Int i_end   = JobPointer( point_count, thread_count, thread + 1 );
                
                for( Int i = i_begin; i < i_end; ++ i )
                {
                    CoordToAxes( &X[n*i], a );
                    
                    AxesToHilbert( a, M[i] );
                }
            },
            thread_count
        );
        
        ptoc(ClassName()+"::CoordToHilbert");
        
        return M;
    }
