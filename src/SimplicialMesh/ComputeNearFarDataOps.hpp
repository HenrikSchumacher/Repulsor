protected:
    
    void ComputeNearFarDataOps(
        Tensor2<Real,Int> & restrict P_coords,
        Tensor3<Real,Int> & restrict P_hull_coords,
        Tensor2<Real,Int> & restrict P_near,
        Tensor2<Real,Int> & restrict P_far,
        SparseMatrix_T    & restrict DiffOp,
        SparseMatrix_T    & restrict AvOp
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
    
            
                PrimitiveDataKernel ker ( *this );
                
                const Int simplex_count = simplices.Dimension(0);
                const Int i_begin = JobPointer<Int>(simplex_count, ThreadCount(), thread  );
                const Int i_end   = JobPointer<Int>(simplex_count, ThreadCount(), thread+1);
                
                // Loop over all simplices.
                for( Int i = i_begin; i < i_end; ++i )
                {
                    ker.ReadPrimitive(i);
                    
                    Av_outer[i+1] = (i+1) * SIZE;
                    
                    ker.s_simplex.Write( &Av_inner[SIZE * i] );
                    fill_buffer<SIZE>( &Av_value[SIZE * i], nth );
                    
                    for( Int k = 0; k < AMB_DIM; ++k )
                    {
                        const Int row = AMB_DIM * i + k;
                        
                        const Int rp = row * SIZE;
                        
                        Diff_outer[row+1] = rp + SIZE;
                        
                        ker.s_simplex.Write( &Diff_inner[rp] );
                    }

                    ker.hull.Write   ( P_hull_coords.data(i) );
                    ker.center.Write ( P_coords.data(i) );
                    ker.D_f.Write    ( &Diff_value[ HULL_SIZE * i ] );
                    
                    ker.WriteNear    ( P_near.data(i) );
                    ker.WriteFar     ( P_far.data(i)  );
                }
            },
            ThreadCount()
        );
        
        ptoc(ClassName()+"::ComputeNearFarDataOps");
    };
