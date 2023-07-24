public:
        
    void TangentialSmoothing( const Int max_iter = 1 ) override
    {
        
        ptic(ClassName()+"::TangentialSmoothing");
        
        // For each vertex v compute x_new[v] = (id - P) * ( mean( neighbor x[s] ) - x[v] );
        
        // Here P is the volume-weighted vertex normal projector.
        
        // mean( neighbor x[s] ) is the mean of centers of simplices containing v
        
//            if( !compressedQ )
//            {
//                eprint(ClassName()+"::VertexTangentialProjectors: run this routine only after compressing data structure! Aborting.");
//                ptoc(ClassName()+"::TangentialSmoothing");
//                return;
//            }
        
        Compress();

        // Compute simplex volumes, centers and tangential projectors.
        Tensor1<Real,Int> S_vol                   ( simplex_count                   );
        Tensor2<Real,Int> S_centers               ( simplex_count, AMB_DIM          );
        Tensor3<Real,Int> S_tangential_projectors ( simplex_count, AMB_DIM, AMB_DIM );
        
        for( Int iter = 0; iter < max_iter; ++iter )
        {
            ParallelDo(
                [&,this]( const Int thread )
                {
                    const Simplex_T s_begin = JobPointer( simplex_count, thread_count, thread    );
                    const Simplex_T s_end   = JobPointer( simplex_count, thread_count, thread +1 );

                    SimplexDataKernel<DOM_DIM,AMB_DIM,Real,Int> ker( V_coords, simplices, V_charges );
                    
                    Tiny::Matrix<AMB_DIM,AMB_DIM,Real,Int> Q;
                    
                    for( Simplex_T s = s_begin; s < s_end; ++s )
                    {
                        ker.ReadSimplex(s);
                        
                        Q.SetIdentity();
                        
                        Q -= ker.NormalProjector();
                        
                        Q *= ker.Volume();
                        
                        S_vol[s] = ker.Volume();
                        
                        ker.Center().Write( S_centers.data(s) );
                        
                        Q.Write( S_tangential_projectors.data(s) );
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
                    Vector_T mean;
                    Matrix_T Q;
                    Matrix_T R;

                    for( Vertex_T v = v_begin; v < v_end; ++v )
                    {
                        Real a = 0;
                        mean.SetZero();
                        Q.SetZero();
                        
                        for( Simplex_T s : V_parent_simplices[v] )
                        {
                            a += S_vol[s];
                            
                            x.Read( S_centers.data(s) );
                            
                            mean += x;
                            
                            R.Read( S_tangential_projectors.data(s) );

                            Q += R;
                        }
                        
                        // Beware, Q is only a coarse approximation to the tangential projector!
                        // We would have to perform PCA to do it correctly.
                        Q /= a;
                        
                        mean /=  static_cast<Real>(V_parent_simplices[v].Size());
                        // Now mean is the mean of the simplex centers.
                        
                        
                        x.Read( V_coords.data(v) );

                        mean -= x;

                        mean *= static_cast<Real>(1.5);
                        
                        Dot<AddTo>( Q, mean, x );

                        x.Write( V_coords.data(v) );
                    }
                },
                thread_count
            );
        }
        
        ptoc(ClassName()+"::TangentialSmoothing");
    }
