protected:
    
    void ComputeNearFarDataOps(
        mref<Tensor2<Real,Int>> P_coords,
        mref<Tensor3<Real,Int>> P_hull_coords,
        mref<Tensor2<Real,Int>> P_near,
        mref<Tensor2<Real,Int>> P_far,
        mref<SparseMatrix_T>    DiffOp,
        mref<SparseMatrix_T> AvOp
    ) const
    {
        ptic(ClassName()+"::ComputeNearFarDataOps");
        
        ParallelDo(
            [&,this]( const Int thread )
            {
                
                mptr<LInt> Av_outer   = AvOp.Outer().data();
                mptr< Int> Av_inner   = AvOp.Inner().data();
                mptr<Real> Av_value   = AvOp.Values().data();
    
                mptr<LInt> Diff_outer = DiffOp.Outer().data();
                mptr< Int> Diff_inner = DiffOp.Inner().data();
                mptr<Real> Diff_value = DiffOp.Values().data();
    
                SimplexDataKernel_T ker ( V_coords, simplices, V_charges );
                
                const Int simplex_count = simplices.Dimension(0);
                const Int i_begin = JobPointer<Int>(simplex_count, ThreadCount(), thread  );
                const Int i_end   = JobPointer<Int>(simplex_count, ThreadCount(), thread+1);
                
                // Loop over all simplices.
                for( Int i = i_begin; i < i_end; ++i )
                {
                    ker.ReadSimplex(i);
                    
                    Av_outer[i+1] = (i+1) * SIZE;
                    
                    ker.WriteSortedSimplex( &Av_inner[SIZE * i] );
                    fill_buffer<SIZE>( &Av_value[SIZE * i], nth );
                    
                    for( Int k = 0; k < AMB_DIM; ++k )
                    {
                        const Int row = AMB_DIM * i + k;
                        
                        const Int rp = row * SIZE;
                        
                        Diff_outer[row+1] = rp + SIZE;
                        
                        ker.WriteSortedSimplex( &Diff_inner[rp] );
                    }

                    ker.Hull().Write   ( P_hull_coords.data(i) );
                    ker.Center().Write ( P_coords.data(i) );
                    ker.WriteDiffOp    ( &Diff_value[ HULL_SIZE * i ] );
                    
                    ker.WriteNear      ( P_near.data(i) );
                    ker.WriteFar       ( P_far.data(i)  );
                }
            },
            ThreadCount()
        );
        
        ptoc(ClassName()+"::ComputeNearFarDataOps");
    };
