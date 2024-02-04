public:

//#######################################################################################
//####                                AxesToMorton                                  #####
//#######################################################################################

    Tensor1<Morton_T,Int> AxesToMorton( cref<Tensor1<Axes_T,Int>> A  )
    {
        ptic(ClassName()+"::AxesToMorton");
        
        Tensor1<Morton_T,Int> M ( A.Dimension(0) );
        
        ParallelDo(
            [&]( const Int thread )
            {
                const Int k_begin = JobPointer( A.Dimension(0), thread_count, thread     );
                const Int k_end   = JobPointer( A.Dimension(0), thread_count, thread + 1 );

                for( Int k = k_begin; k < k_end; ++k )
                {   
                    AxesToMorton( A[k], M[k] );
                }
            },
            thread_count
        );
        
        ptoc(ClassName()+"::AxesToMorton");
        
        return M;
    }

    // The actual worker routine.
    void constexpr AxesToMorton_Reference( cref<Axes_T> a, mref<Morton_T> m ) const
    {
    //        Not needed because we write all bits.
    //        zerofy_buffer(M, n);

        for( Int j = 0; j < bit_count; ++j )
        {
            const Int nj = n * j;

            for( Int i = 0; i < n; ++i )
            {
                const Int i_ = (nj + i) / bit_count;
                const Int j_ = (nj + i) % bit_count;

                set_bit( m[n-1-i_], j_, get_bit( a[i], j ) );
            }
        }
    }


    // The actual worker routine.
    void constexpr AxesToMorton( cref<Axes_T> a, mref<Morton_T> m ) const
    {
        constexpr Int shift = 8 * n;
        
        BitField_T M = 0;

        for( Int i = 0; i < n; ++i )
        {
            const UInt a_i = a[i];
            
            for( Int j = 0; j < byte_count; ++j )
            {
                M |= LUT[get_byte(a_i,j)] << (shift * j + i);
            }
        }
        
        constexpr BitField_T mask ( ~UInt(0) );
        
        for( Int i = 0; i < n; ++i )
        {
            m[n-1-i] = static_cast<UInt>( ( M >> (bit_count * i) & mask ).to_ullong() );
        }
    }
