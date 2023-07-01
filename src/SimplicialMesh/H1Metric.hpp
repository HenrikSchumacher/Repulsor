protected:
    
    virtual SparseMatrix_T H1Metric( const Real c_1, const Real c_0 ) const override
    {
        Tensor2<Int, Int> ilist ( SimplexCount(), SIZE * SIZE );
        Tensor2<Int, Int> jlist ( SimplexCount(), SIZE * SIZE );
        Tensor2<Real,Int> alist ( SimplexCount(), SIZE * SIZE );
        
        ParallelDo(
            [this,&ilist,&jlist,&alist,c_1,c_0]( const Int thread )
            {
                const Int n = SimplexCount();
                const Int k_begin = JobPointer<Int>( n, ThreadCount(), thread     );
                const Int k_end   = JobPointer<Int>( n, ThreadCount(), thread + 1 );
                
                Tiny::Vector<SIZE,Int,Int> simplex;
                Tiny::Matrix<AMB_DIM,DOM_DIM,Real,Int> df;
                Tiny::Matrix<DOM_DIM,AMB_DIM,Real,Int> dftransp;
                Tiny::Matrix<DOM_DIM,DOM_DIM,Real,Int> g_inv;
                Tiny::Matrix<DOM_DIM,DOM_DIM,Real,Int> id;
                
                for( Int i = 0; i < DOM_DIM; ++i )
                {
                    for( Int j = 0; j < DOM_DIM; ++j )
                    {
                        id[i][j] = Int(i==j);
                    }
                }
                
                Tiny::Matrix<SIZE,SIZE,Real,Int> mass;
                
                const Real m_factor = c_0 / StandardSimplexVolume() / Factorial(static_cast<Real>(DOM_DIM + 2));
                
                for( Int i = 0; i < SIZE; ++i )
                {
                    for( Int j = 0; j < SIZE; ++j )
                    {
                        mass[i][j] = m_factor * (Scalar::One<Real> + Delta<Real>(i,j) );
                    }
                }
                
                Tiny::SelfAdjointMatrix<DOM_DIM,Real,Int> g;
                
                Tiny::Matrix<SIZE,SIZE, Int,Int> idx;
                Tiny::Matrix<SIZE,SIZE, Int,Int> jdx;
                
                for( Int k = k_begin; k < k_end; ++k )
                {
                    simplex.Read( simplices.data(k) );
                    
                    Tiny::Matrix<SIZE,SIZE,Real,Int> val ( Scalar::Zero<Real> );
                    
                    for( Int i = 0; i < AMB_DIM; ++i )
                    {
                        for( Int j = 0; j < DOM_DIM; ++j )
                        {
                            df(i,j) = V_coords(simplex[j+1],i) - V_coords(simplex[0],i);
                        }
                    }
                    
                    df.Transpose( dftransp );

                    // g = df^T * df.
                    Dot<Overwrite>( dftransp, df, g_inv );
                    
                    g.Read( g_inv.data() );
                    
                    g.Cholesky();
                    
                    // Compute simplex area
                    Real a = StandardSimplexVolume();
                    
                    for( Int l = 0; l < DOM_DIM; ++l )
                    {
                        a *= g[l][l];
                    }
                    
                    g_inv.Read( id.data() );
                    
                    g.CholeskySolve( g_inv );
                    
                    const Real l_factor = a * c_1;
                    
                    for( Int i = 0; i < DOM_DIM; ++i )
                    {
                        for( Int j = 0; j < DOM_DIM; ++j )
                        {
                            const Real a_ij = l_factor * g_inv[i][j];
                            
                            val[0  ][0  ] += a_ij;
                            val[i+1][0  ] -= a_ij;
                            val[0  ][j+1] -= a_ij;
                            val[i+1][j+1]  = a_ij;
                        }
                    }
                    
                    combine_buffers<
                        Scalar::Flag::Generic, Scalar::Flag::Plus, SIZE * SIZE
                    >(
                        a,                 mass.data(),
                        Scalar::One<Real>, val.data()
                    );
                    
                    val.Write( alist.data(k) );
                    
                    for( Int i = 0; i < SIZE; ++i )
                    {
                        for( Int j = 0; j < SIZE; ++j )
                        {
                            idx[i][j] = simplex[i];
                            jdx[i][j] = simplex[j];
                        }
                    }
                    
                    idx.Write( ilist.data(k) );
                    jdx.Write( jlist.data(k) );
                }
            },
            ThreadCount()
        );
        
        SparseMatrix_T A (
            LInt(SimplexCount()) * LInt(SIZE * SIZE),
            ilist.data(),  jlist.data(),  alist.data(),
            VertexCount(), VertexCount(), ThreadCount(),
            true, false
        );
        
        return A;
    }
    
public:
    
    virtual const SparseMatrix_T & StiffnessMatrix() const override
    {
        std::string tag ("StiffnessMatrix");
        if( !this->InCacheQ(tag))
        {
            ptic(ClassName()+"::"+tag);
            
            this->SetCache( tag,
               std::any(std::move(H1Metric(1,0)))
            );
            
            ptoc(ClassName()+"::"+tag);
        }
        
        return std::any_cast<SparseMatrix_T &>( this->GetCache(tag) );
    }
    
    virtual const SparseMatrix_T & MassMatrix() const override
    {
        std::string tag ("MassMatrix");
        if( !this->InCacheQ(tag))
        {
            ptic(ClassName()+"::"+tag);
            
            this->SetCache( tag,
                std::any(std::move(H1Metric(0,1)))
            );
            
            ptoc(ClassName()+"::"+tag);
        }
        
        return std::any_cast<SparseMatrix_T &>( this->GetCache(tag) );
    }