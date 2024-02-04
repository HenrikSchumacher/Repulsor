public:

//#######################################################################################
//####                                AxesToHilbert                                 #####
//#######################################################################################

    Tensor1<Hilbert_T,Int> AxesToHilbert( cref<Tensor1<Axes_T,Int>> A  )
    {
        ptic(ClassName()+"::AxesToHilbert");
        
        Tensor1<Morton_T,Int> H ( A.Dimension(0) );
        
        ParallelDo(
            [&]( const Int thread )
            {
                const Int k_begin = JobPointer( A.Dimension(0), thread_count, thread     );
                const Int k_end   = JobPointer( A.Dimension(0), thread_count, thread + 1 );

                for( Int k = k_begin; k < k_end; ++k )
                {
                    AxesToHilbert( A[k], H[k] );
                }
            },
            thread_count
        );
        
        ptoc(ClassName()+"::AxesToHilbert");
        
        return H;
    }

// Port of Paul Chernoch's code on https://stackoverflow.com/a/10384110/8248900.
// See also this page for the whole project: https://github.com/paulchernoch/HilbertTransformation


    void constexpr AxesToHilbert( cref<Axes_T> a, mref<Morton_T> h ) const
    {
        Axes_T x = a;

        //See David Moton's Java implementation SmallHilbertCurve.java in https://github.com/davidmoten/hilbert-curve.
        
        constexpr UInt q_max ( one << (bit_count - 1) );

        // Inverse undo
        for( UInt q = q_max; q > one; q >>= 1 )
        {
            const UInt p = q - 1;
            
            for( Int i = 0; i < n; ++i )
            {
                if( (x[i] & q) != 0 )
                {
                    x[0] ^= p; // invert
                }
                else
                {
                    const UInt t = (x[0] ^ x[i]) & p;
                    x[0] ^= t;
                    x[i] ^= t;
                }
            }
        } // exchange
        
        // Gray encode
        for( Int i = 1; i < n; i++ )
        {
            x[i] ^= x[i-1];
        }
        UInt t = 0;
        
        for( UInt q = q_max; q > one; q >>= 1 )
        {
            if( (x[n-1] & q) != 0 )
            {
                t ^= q - 1;
            }
        }
        
        Axes_T y;
        
        for( Int i = 0; i < n; ++i )
        {
            y[n-1-i] = x[i] ^ t;
        }
        
        // Untranspose
        
        AxesToMorton( y, h );
    }
