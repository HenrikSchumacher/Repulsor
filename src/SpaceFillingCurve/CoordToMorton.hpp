//#######################################################################################
//####                                CoordToMorton                                 #####
//#######################################################################################

    Tensor1<Morton_T,Int> CoordToMorton( cptr<Real> X, const Int point_count )
    {
        ComputeBoundingBox( X, point_count );
        
        ptic(ClassName()+"::CoordToMorton");
        
        Tensor1<Morton_T,Int> M ( point_count );
        
        ParallelDo(
            [&,this]( const Int thread )
            {
                Axes_T a;
                
                const Int i_begin = JobPointer<Int>( point_count, thread_count, thread     );
                const Int i_end   = JobPointer<Int>( point_count, thread_count, thread + 1 );
                
                for( Int i = i_begin; i < i_end; ++ i )
                {
                    CoordToAxes( &X[n*i], a );

//                    dump(static_cast<Int>(a[0]));
//                    dump(static_cast<Int>(a[1]));
//                    dump(static_cast<Int>(a[2]));
//                    
//                    AxesToMorton_Reference( a, M[i] );
//
//                    valprint("correct",static_cast<Int>(M[i][0]));
//                    valprint("correct",static_cast<Int>(M[i][1]));
//                    valprint("correct",static_cast<Int>(M[i][2]));

                    AxesToMorton( a, M[i] );
//
//                    dump(static_cast<Int>(M[i][0]));
//                    dump(static_cast<Int>(M[i][1]));
//                    dump(static_cast<Int>(M[i][2]));
                    
                }
            },
            thread_count
        );
        
        ptoc(ClassName()+"::CoordToMorton");
        
        return M;
    }
