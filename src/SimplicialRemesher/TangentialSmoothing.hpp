public:
        
    void TangentialSmoothing( const Int max_iter = 1, const Real step_size = Frac<Real>(1,2) ) override
    {
        
        TOOLS_PTIMER(timer,ClassName()+"::TangentialSmoothing");
        
        // For each vertex v compute x_new[v] = (id - P) * ( mean( neighbor x[s] ) - x[v] );
        
        // Here P is the volume-weighted vertex normal projector.
        
        // mean( neighbor x[s] ) is the mean of centers of simplices containing v
        
        Compress();

        // Compute simplex volumes, centers and tangential projectors.
        Tensor1<Real,Int> S_vol                            ( simplex_count                   );
        Tensor2<Real,Int> S_center_times_vol               ( simplex_count, AMB_DIM          );
        Tensor3<Real,Int> S_tangential_projector_times_vol ( simplex_count, AMB_DIM, AMB_DIM );
        
        for( Int iter = 0; iter < max_iter; ++iter )
        {
            ParallelDo(
                [&,this]( const Int thread )
                {
                    const Simplex_T s_begin = JobPointer( simplex_count, thread_count, thread    );
                    const Simplex_T s_end   = JobPointer( simplex_count, thread_count, thread +1 );

                    SimplexDataKernel<DOM_DIM,AMB_DIM,Real,Int> ker( V_coords, simplices, V_charges );
                    
                    Vector_T x;
                    Matrix_T Q;
                    
                    for( Simplex_T s = s_begin; s < s_end; ++s )
                    {
                        ker.ReadSimplex(s);
                        
                        const Real s_vol = ker.Volume();
                        
                        S_vol[s] = s_vol;
                        
                        x = s_vol * ker.Center();
                        x.Write( S_center_times_vol.data(s) );
                        
                        Q.SetIdentity();
                        Q -= ker.NormalProjector();
                        Q *= s_vol;
                        Q.Write( S_tangential_projector_times_vol.data(s) );
                    }
                },
                thread_count
            );
            
            ParallelDo(
                [&,this]( const Int thread )
                {
                    const Vertex_T v_begin = JobPointer( vertex_count, thread_count, thread    );
                    const Vertex_T v_end   = JobPointer( vertex_count, thread_count, thread +1 );

                    Vector_T x;
                    Vector_T u;
                    Matrix_T Q;
                    Matrix_T R;

                    for( Vertex_T v = v_begin; v < v_end; ++v )
                    {
                        if( VertexPinnedQ(v) ) { continue; }
                        
                        Real a = 0;
                        u.SetZero();
                        Q.SetZero();
                        
                        for( Simplex_T s : V_parent_simplices[v] )
                        {
                            const Real vol_s = S_vol[s];
                            
                            a += vol_s;
                            
                            x.Read( S_center_times_vol.data(s) );
                            u += x;
                            
                            R.Read( S_tangential_projector_times_vol.data(s) );
                            Q += R;
                        }
                        
                        
                        Q /= a;
                        // Now Q is is the average of the tangential projectors.
                        // Beware, Q is only a coarse approximation to the tangential projector!
                        // We would have to perform PCA to do it correctly.
                        
                        u /=  a;
                        // Now u is the average of the simplex centers.
                        
                        x.Read( V_coords.data(v) );
                        u -= x;
                        // Now u is the vector that points from x to the average of the simplex centers.
                        
                        u *= step_size;
                        
                        Dot<AddTo>( Q, u, x );

                        x.Write( V_coords.data(v) );
                    }
                },
                thread_count
            );
        }
    }
