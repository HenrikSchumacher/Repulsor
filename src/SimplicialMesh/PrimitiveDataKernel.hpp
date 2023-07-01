private:

    struct PrimitiveDataKernel
    {
        using Mesh_T = SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>;
        
        const Mesh_T & M;
        static constexpr Int PROJ_DIM  = (AMB_DIM*(AMB_DIM+1))/2;
        
        
        Tiny::Vector<SIZE,Int,Int>   simplex;
        Tiny::Vector<SIZE,Int,Int> s_simplex;
        
        Tiny::Vector<AMB_DIM,Real,Int> center;
        
        Tiny::Matrix<SIZE,   AMB_DIM,Real,Int> hull;
        Tiny::Matrix<AMB_DIM,DOM_DIM,Real,Int> df;
        Tiny::Matrix<DOM_DIM,AMB_DIM,Real,Int> dfdagger;
        Tiny::Matrix<AMB_DIM,AMB_DIM,Real,Int> P;
        
//        Tiny::Matrix< FAR_DIM,AMB_DIM*SIZE,Real,Int> Dfar;
//        Tiny::Matrix<NEAR_DIM,AMB_DIM*SIZE,Real,Int> Dnear;
        
        Tiny::SelfAdjointMatrix<DOM_DIM,Real,Int> g_inv;
        
        Tiny::Matrix<AMB_DIM,SIZE,Real,Int> D_f;

        Real a;
        Real charge;
        
        SortNet<SIZE,Int> sort;
        
        PrimitiveDataKernel( const Mesh_T & M_ )
        :   M ( M_ )
        {
            if constexpr ( DOM_DIM == 0 )
            {
                D_f.SetZero();
                P.SetZero();
                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    P[k][k] = Scalar::One<Real>;
                }
            }
        }
        
        ~PrimitiveDataKernel() = default;
        
        void ReadPrimitive( const Int i )
        {
              simplex.Read( M.simplices.data(i) );
            s_simplex.Read( simplex.data()      );
          
            // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
            sort( &s_simplex[0] );
            
            for( Int l = 0; l < SIZE; ++l )
            {
                copy_buffer<AMB_DIM>( M.V_coords.data(simplex[l]), hull[l] );
            }
            
            copy_buffer<AMB_DIM>(hull[0],center.data());
            
            charge = M.V_charges[simplex[0]];
            for( Int l = 1; l < SIZE; ++l )
            {
                add_to_buffer<AMB_DIM>( hull[l], center.data() );
                
                charge += M.V_charges[simplex[l]];
            }
            center *= nth;
            charge *= nth;
            
            a = charge * Mesh_T::StandardSimplexVolume();
            
            if constexpr ( DOM_DIM > 0 )
            {
                for( Int l = 0; l < DOM_DIM; ++l )
                {
                    for( Int k = 0; k < AMB_DIM; ++k )
                    {
                        df[k][l] =
                        M.V_coords[s_simplex[l+1]][k]
                        -
                        M.V_coords[s_simplex[0  ]][k];
                    }
                }
                
                df.Transpose( dfdagger );
                
            // Compute inverse metric.
                
                // g = df^T * df.
                // At the moment dfdagger is just the transpose of df.
                Tiny::gemm<Op::Id,Op::Id,DOM_DIM,DOM_DIM,AMB_DIM,
                Scalar::Flag::Plus,Scalar::Flag::Zero
                >(
                    Scalar::One<Real>,  dfdagger.data(), AMB_DIM,
                                        df.data(),       DOM_DIM,
                    Scalar::Zero<Real>, g_inv.data(),    DOM_DIM
                );
                
                // Factorize g in place.
                g_inv.Cholesky();
                
                for( Int l = 0; l < DOM_DIM; ++l )
                {
                    a *= g_inv[l][l];
                }
                
                
            // Compute pseudo-inverse of df
                
                //  dfdagger = g^{-1} * df^T
                g_inv.CholeskySolve( dfdagger );

                
            // Compute normal projektor.
                
                // P = id - df * dfdagger
                Tiny::gemm<Op::Id,Op::Id,AMB_DIM,AMB_DIM,DOM_DIM,
                Scalar::Flag::Minus,Scalar::Flag::Zero
                >(
                    -Scalar::One<Real>,  df.data(),       DOM_DIM,
                                         dfdagger.data(), AMB_DIM,
                     Scalar::Zero<Real>, P.data(),        AMB_DIM
                );
                
                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    P[k][k] += Scalar::One<Real>;
                }
                
            // Compute derivative operator (AMB_DIM x SIZE matrix).
                                
                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    D_f[k][0] = 0;
                    
                    if constexpr ( DOM_DIM > 0 )
                    {
                        for( Int l = 0; l < DOM_DIM; ++l )
                        {
                            D_f[k][0  ] -= dfdagger[l][k];
                            D_f[k][l+1]  = dfdagger[l][k];
                        }
                    }
                }
            }
        }
        
        void WriteNear( mut<Real> near ) const
        {
            near[0] = a;

            hull.Write( &near[1] );
            
            LOOP_UNROLL_FULL
            for( Int k = 0; k < PROJ_DIM; ++k )
            {
                near[1 + AMB_DIM * SIZE + k] = P[ tri_i<AMB_DIM>(k) ][ tri_j<AMB_DIM>(k) ];
            }
        }
        
        void WriteFar( mut<Real> far ) const
        {
            far[0] = a;

            center.Write( &far[1] );
            
            LOOP_UNROLL_FULL
            for( Int k = 0; k < PROJ_DIM; ++k )
            {
                far[1 + AMB_DIM + k] = P[ tri_i<AMB_DIM>(k) ][ tri_j<AMB_DIM>(k) ];
            }
        }
    };
