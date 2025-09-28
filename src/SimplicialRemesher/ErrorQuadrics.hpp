public:

    void ComputeErrorQuadrics()
    {
        TOOLS_PTIMER(timer,ClassName()+"::ComputeErrorQuadrics");
        // Error quadrics a la Garland-Heckbert.
        
        Compress();
        
        S_quadrics = Tensor3<Real,Int> ( simplex_count, AMB_DIM+1, AMB_DIM+1 );

        // Compute error quadrics for all faces.
        ParallelDo(
            [&,this]( const Int thread )
            {
                const Int i_begin = JobPointer( simplex_count, thread_count, thread    );
                const Int i_end   = JobPointer( simplex_count, thread_count, thread +1 );

                SimplexDataKernel<DOM_DIM,AMB_DIM,Real,Int> ker( V_coords, simplices, V_charges );

                for( Int i = i_begin; i < i_end; ++i )
                {
                    ker.ReadSimplex(i);

                    ker.WriteErrorQuadric( S_quadrics.data(i) );
                }
            },
            thread_count
        );

        V_quadrics = Tensor3<Real,Int> ( simplex_count, AMB_DIM+1, AMB_DIM+1 );

        ParallelDo(
            [&,this]( const Int thread )
            {
                const Vertex_T v_begin = JobPointer( vertex_count, thread_count, thread    );
                const Vertex_T v_end   = JobPointer( vertex_count, thread_count, thread +1 );

                Quadric_T Q;
                Quadric_T R;

                for( Vertex_T v = v_begin; v < v_end; ++v )
                {
                    Q.SetZero();

                    for( Simplex_T s : V_parent_simplices[v] )
                    {
                        R.Read( S_quadrics.data(s) );

                        Q += R;
                    }

                    Q.Write( V_quadrics.data(v) );
                }
            },
            thread_count
        );
    }
