protected:
    
    void ComputeNearFarData(
        Tensor2<Real,Int> & P_near,
        Tensor2<Real,Int> & P_far
    ) const
    {
        ptic(ClassName()+"::ComputeNearFarData");
        
        ParallelDo(
            [&]( const Int thread )
            {
                Tiny::Vector<AMB_DIM,Real,Int> center;
                
                Tiny::Matrix<SIZE,   AMB_DIM,Real,Int> hull;
                Tiny::Matrix<AMB_DIM,DOM_DIM,Real,Int> df;
                Tiny::Matrix<DOM_DIM,AMB_DIM,Real,Int> dfdagger;
                Tiny::Matrix<AMB_DIM,AMB_DIM,Real,Int> P;
                
                if constexpr ( DOM_DIM == 0 )
                {
                    P.SetZero();
                    for( Int k = 0; k < AMB_DIM; ++k )
                    {
                        P[k][k] = Scalar::One<Real>;
                    }
                }
                
                Tiny::SelfAdjointMatrix<DOM_DIM,Real,Int> g;
                
                Tiny::Vector<SIZE,Int,Int> simplex;
                Tiny::Vector<SIZE,Int,Int> s_simplex;   // Sorted simplex.
                
                const Int i_begin = JobPointer<Int>(SimplexCount(),ThreadCount(),thread  );
                const Int i_end   = JobPointer<Int>(SimplexCount(),ThreadCount(),thread+1);
    
                SortNet<SIZE,Int> sort;
                
                for( Int i = i_begin; i < i_end; ++i )
                {
                    mut<Real> near   = P_near.data(i);
                    mut<Real> far    = P_far.data(i);
                    
                      simplex.Read( simplices.data(i) );
                    s_simplex.Read( simplex.data()    );
                  
                    // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
                    sort( &s_simplex[0] );

                    for( Int l = 0; l < SIZE; ++l )
                    {
                        copy_buffer<AMB_DIM>( V_coords.data(simplex[l]), hull[l] );
                    }

                    hull.Write( &near[1] );
                    
                    Real charge = V_charges[simplex[0]];
                    copy_buffer<AMB_DIM>(hull[0],center.data());
                    for( Int l = 1; l < SIZE; ++l )
                    {
                        add_to_buffer<AMB_DIM>( hull[l], center.data() );
                        
                        charge += V_charges[simplex[l]];
                    }
                    center *= nth;
                    charge *= nth;
                    
                    center.Write( &far[1] );
                    
                    P_charges[i] = charge;
                    
                    Real a = charge * StandardSimplexVolume();
                    
                    if constexpr ( DOM_DIM > 0 )
                    {
                        for( Int l = 0; l < DOM_DIM; ++l )
                        {
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                df[k][l] =
                                    V_coords[s_simplex[l+1]][k]
                                    -
                                    V_coords[s_simplex[0  ]][k];
                            }
                        }
                        
                        df.Transpose( dfdagger );
                        
                        // g = df^T * df.
                        // At the moment dfdagger is just the transpose of df.
                        Tiny::gemm<Op::Id,Op::Id,DOM_DIM,DOM_DIM,AMB_DIM,
                            Scalar::Flag::Plus,Scalar::Flag::Zero
                        >(
                            Scalar::One<Real>,  dfdagger.data(), AMB_DIM,
                                                df.data(),       DOM_DIM,
                            Scalar::Zero<Real>, g.data(),        DOM_DIM
                        );
                        
                        // Factorize g in place.
                        g.Cholesky();
                        
                        for( Int l = 0; l < DOM_DIM; ++l )
                        {
                            a *= g[l][l];
                        }
                    }
                    
                    far[0] = near[0] = a;
                    
                    if constexpr( DOM_DIM > 0 )
                    {
                        //  dfdagger = g^{-1} * df^T
                        g.CholeskySolve( dfdagger );
                        
                        // P = id - df * dfdagger
                        Tiny::gemm<Op::Id,Op::Id,AMB_DIM,AMB_DIM,DOM_DIM,
                        Scalar::Flag::Minus,Scalar::Flag::Zero
                        >(
                            -Scalar::One<Real>, df.data(),       DOM_DIM,
                                                dfdagger.data(), AMB_DIM,
                            Scalar::Zero<Real>, P.data(),        AMB_DIM
                        );
                        
                        for( Int k = 0; k < AMB_DIM; ++k )
                        {
                            P[k][k] += Scalar::One<Real>;
                        }
                    }
                    
                    
                    {
                        // TODO: variable pos introduces a bad dependency in loop.
                        
                        Int pos = 0;
                        
                        for( Int k = 0; k < AMB_DIM; ++k )
                        {
                            for( Int l = k; l < AMB_DIM; ++l )
                            {
                                ++pos;
                                near[AMB_DIM * SIZE + pos] = far[AMB_DIM + pos] = P[k][l];
                            }
                        }
                    }

                }
            },
            ThreadCount()
        );
    
        ptoc(ClassName()+"::ComputeNearFarData");
};
