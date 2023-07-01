

protected:
    
    template<bool addto = true>
    void DFarToHulls(
        Tensor3<Real,Int> & restrict buffer,
    ) const
    {
        ptic(ClassName()+"::DFarToHulls");
        
        if constexpr (!addtoQ)
        {
            buffer.SetZero();
        }
        
        ParallelDo(
            [&]( const Int thread )
            {
                Tiny::Vector<AMB_DIM,Real,Int> center;
                
                Tiny::Matrix<SIZE,   AMB_DIM,Real,Int> hull;
                Tiny::Matrix<AMB_DIM,DOM_DIM,Real,Int> df;
                Tiny::Matrix<DOM_DIM,AMB_DIM,Real,Int> dfdagger;
                Tiny::Matrix<AMB_DIM,AMB_DIM,Real,Int> P;
                
                Tiny::Matrix< FAR_DIM,AMB_DIM*SIZE,Real,Int> Dfar;
                Tiny::Matrix<NEAR_DIM,AMB_DIM*SIZE,Real,Int> Dnear;
                
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
                
                const Int i_begin = JobPointer<Int>(SimplexCount(),ThreadCount(),thread  );
                const Int i_end   = JobPointer<Int>(SimplexCount(),ThreadCount(),thread+1);
                
                for( Int i = i_begin; i < i_end; ++i )
                {

                    for( Int l = 0; l < SIZE; ++l )
                    {
                        copy_buffer<AMB_DIM>( V_coords.data(simplex[l]), hull[l] );
                    }
                    
//                    Real charge = V_charges[simplex[0]];
                    copy_buffer<AMB_DIM>(hull[0],center.data());
                    for( Int l = 1; l < SIZE; ++l )
                    {
                        add_to_buffer<AMB_DIM>( hull[l], center.data() );
                        
                        charge += V_charges[simplex[l]];
                    }
                    center *= nth;
//
//                    const Real charge = P_charges[i];
                    
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
                }
            },
            ThreadCount()
        );
    
        ptoc(ClassName()+"::DFarToHulls");
};
