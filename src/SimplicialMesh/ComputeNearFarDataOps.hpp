protected:
    
    void ComputeNearFarDataOps(
              Tensor2<Real,Int> & P_coords,
              Tensor3<Real,Int> & P_hull_coords,
              Tensor2<Real,Int> & P_near,
              Tensor2<Real,Int> & P_far,
        SparseMatrix_T & DiffOp,
        SparseMatrix_T & AvOp
    ) const
    {
        ptic(ClassName()+"::ComputeNearFarDataOps");
        
        ParallelDo(
            [&,this]( const Int thread )
            {
                
                mut<LInt> Av_outer = AvOp.Outer().data();
                mut< Int> Av_inner = AvOp.Inner().data();
                mut<Real> Av_value = AvOp.Values().data();
    
                mut<LInt> Diff_outer = DiffOp.Outer().data();
                mut< Int> Diff_inner = DiffOp.Inner().data();
                mut<Real> Diff_value = DiffOp.Values().data();
    
                Tiny::Vector<AMB_DIM,Real,Int> center;
                
                Tiny::Matrix<SIZE,   AMB_DIM,Real,Int> hull;
                Tiny::Matrix<AMB_DIM,SIZE,   Real,Int> Df;
                
                Tiny::Matrix<AMB_DIM,DOM_DIM,Real,Int> df;
                Tiny::Matrix<DOM_DIM,AMB_DIM,Real,Int> dfdagger;
                Tiny::Matrix<AMB_DIM,AMB_DIM,Real,Int> P;
                
                if constexpr ( DOM_DIM == 0 )
                {
                    Df.SetZero();
                    P.SetZero();
                    for( Int k = 0; k < AMB_DIM; ++k )
                    {
                        P[k][k] = Scalar::One<Real>;
                    }
                }
                
                Tiny::SelfAdjointMatrix    <DOM_DIM,Real,Int> g;
                
                Tiny::Vector<SIZE,Int,Int> simplex;
                Tiny::Vector<SIZE,Int,Int> s_simplex;   // Sorted simplex.
    
            
                Real std_simplex_volume = 1;
                for( Int l = 2; l < SIZE; ++l )
                {
                    std_simplex_volume *= l;
                }
                std_simplex_volume = Scalar::Inv<Real>(std_simplex_volume);
                
//                    dump(std_simplex_volume);
                
                const Int simplex_count = simplices.Dimension(0);
                const Int i_begin = JobPointer<Int>(simplex_count, ThreadCount(), thread  );
                const Int i_end   = JobPointer<Int>(simplex_count, ThreadCount(), thread+1);
    
                SortNet<SIZE,Int> sort;
                
                // Loop over all simplices.
                for( Int i = i_begin; i < i_end; ++i )
                {
                    mut<Real> near   = P_near.data(i);
                    mut<Real> far    = P_far.data(i);
                    
                      simplex.Read( simplices.data(i) );
                    s_simplex.Read( simplex.data() );

//                        dump(simplex);
                  
                    // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
                    sort( &s_simplex[0] );

//                        dump(s_simplex);
                    
                    Av_outer[i+1] = (i+1) * SIZE;
                    
                    s_simplex.Write( &Av_inner[SIZE * i] );
                    fill_buffer<SIZE>( &Av_value[SIZE * i], nth );
                    
                    for( Int k = 0; k < AMB_DIM; ++k )
                    {
                        const Int row = AMB_DIM * i + k;
                        
                        const Int rp = row * SIZE;
                        
                        Diff_outer[row+1] = rp + SIZE;
                        
                        s_simplex.Write( &Diff_inner[rp] );
                    }
       
                    for( Int l = 0; l < SIZE; ++l )
                    {
                        copy_buffer<AMB_DIM>( V_coords.data(simplex[l]), hull[l] );
                    }
                    
//                        dump(hull);

                    hull.Write( P_hull_coords.data(i) );
                    hull.Write( &near[1]              );
                    
                    
                    Real charge = V_charges[simplex[0]];
                    copy_buffer<AMB_DIM>(hull[0],center.data());
                    for( Int l = 1; l < SIZE; ++l )
                    {
                        add_to_buffer<AMB_DIM>( hull[l], center.data() );
                        
                        charge += V_charges[simplex[l]];
                    }
                    center *= nth;
                    charge *= nth;
                    
                    P_charges[i] = charge;
                    
                    center.Write( &far[1] );
                    center.Write( P_coords.data(i) );
                    
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
                        
                        //                        dump(df);
                        
                        df.Transpose( dfdagger );
                        
                        //                        dump(dfdagger);
                        
                        // g = df^T * df.
                        // At the moment dfdagger is just the transpose of df.
                        Tiny::gemm<Op::Id,Op::Id,DOM_DIM,DOM_DIM,AMB_DIM,
                        Scalar::Flag::Plus,Scalar::Flag::Zero
                        >(
                          Scalar::One<Real>,  dfdagger.data(), AMB_DIM,
                          df.data(),       DOM_DIM,
                          Scalar::Zero<Real>, g.data(),        DOM_DIM
                          );
                        
                        //                        dump(g);
                        
                        // Factorize g in place.
                        g.Cholesky();
                        
                        //                        dump(g);
                        
                        for( Int l = 0; l < DOM_DIM; ++l )
                        {
                            a *= g[l][l];
                        }
                    }
                    
                    far[0] = near[0] = a;
                    
//                        dump(near[0]);
                    
                    if constexpr ( DOM_DIM > 0 )
                    {
                        //  dfdagger = g^{-1} * df^T
                        g.CholeskySolve( dfdagger );

                        //                        dump(dfdagger);

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
                    }
                    
//                        dump(P);
                    
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
//
//                        for( Int k = 0; k < NearDim(); ++k )
//                        {
//                            valprint("near["+ToString(k)+"]",near[k]);
//                        }
//
//                        for( Int k = 0; k < FAR_DIM; ++k )
//                        {
//                            valprint("far["+ToString(k)+"]",far[k]);
//                        }
                    

                    
                    // Create derivative operator  (AMB_DIM x SIZE matrix).
                    
                    if constexpr ( DOM_DIM > 0 )
                    {
                        for( Int k = 0; k < AMB_DIM; ++k )
                        {
                            Df[k][0] = 0;
                            
                            if constexpr ( DOM_DIM > 0 )
                            {
                                for( Int l = 0; l < DOM_DIM; ++l )
                                {
                                    Df[k][0  ] -= dfdagger[l][k];
                                    Df[k][l+1]  = dfdagger[l][k];
                                }
                            }
                        }
                    }
                    
//                        dump(Df)
                    
                    Df.Write( &Diff_value[ HULL_SIZE * i ] );

                }
            },
            ThreadCount()
        );
        
        ptoc(ClassName()+"::ComputeNearFarDataOps");
    };
    
