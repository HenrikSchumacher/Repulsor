protected:

    void PercolateUp_Tasks( const Int C ) const
    {
        ptic(ClassName()+"::PercolateUp_Tasks");
        switch ( buffer_dim )
        {
            case 1:
            {
                percolateUp_Tasks<1>( C, ThreadCount() );
                break;
            }
            case 2:
            {
                percolateUp_Tasks<2>( C, ThreadCount() );
                break;
            }
            case 3:
            {
                percolateUp_Tasks<3>( C, ThreadCount() );
                break;
            }
            case 4:
            {
                percolateUp_Tasks<4>( C, ThreadCount() );
                break;
            }
            case 5:
            {
                percolateUp_Tasks<5>( C, ThreadCount() );
                break;
            }
            case 6:
            {
                percolateUp_Tasks<6>( C, ThreadCount() );
                break;
            }
            case 7:
            {
                percolateUp_Tasks<7>( C, ThreadCount() );
                break;
            }
            case 8:
            {
                percolateUp_Tasks<8>( C, ThreadCount() );
                break;
            }
            case 9:
            {
                percolateUp_Tasks<9>( C, ThreadCount() );
                break;
            }
            case 10:
            {
                percolateUp_Tasks<10>( C, ThreadCount() );
                break;
            }
            case 11:
            {
                percolateUp_Tasks<11>( C, ThreadCount() );
                break;
            }
            case 12:
            {
                percolateUp_Tasks<12>( C, ThreadCount() );
                break;
            }
            case 13:
            {
                percolateUp_Tasks<13>( C, ThreadCount() );
                break;
            }
            case 14:
            {
                percolateUp_Tasks<14>( C, ThreadCount() );
                break;
            }
            case 15:
            {
                percolateUp_Tasks<15>( C, ThreadCount() );
                break;
            }
            case 16:
            {
                percolateUp_Tasks<16>( C, ThreadCount() );
                break;
            }
            case 17:
            {
                percolateUp_Tasks<17>( C, ThreadCount() );
                break;
            }
            case 18:
            {
                percolateUp_Tasks<18>( C, ThreadCount() );
                break;
            }
            case 19:
            {
                percolateUp_Tasks<19>( C, ThreadCount() );
                break;
            }
            case 20:
            {
                percolateUp_Tasks<20>( C, ThreadCount() );
                break;
            }
            default:
            {
                percolateUp_Tasks_gen( C, ThreadCount() );
            }
        }
        ptoc(ClassName()+"::PercolateUp_Tasks");
    }

    template<Int BUFFER_DIM>
    void percolateUp_Tasks( const Int C, const Int free_thread_count ) const
    {
        // C = cluster index
        
        const Int L = C_left [C];
        const Int R = C_right[C];
        
        if( (L >= 0) && (R >= 0) )
        {
            // If not a leaf, compute the values of the children first.
            #pragma omp task final(free_thread_count<1)  shared( L )
                percolateUp_Tasks<BUFFER_DIM>( L, free_thread_count/2 );
            #pragma omp task final(free_thread_count<1)  shared( R )
                percolateUp_Tasks<BUFFER_DIM>( R, free_thread_count-free_thread_count/2 );
            #pragma omp taskwait
            
            // Aftwards, compute the sum of the two children.
            
//            Real * restrict const c = C_in.data();
            
            const Int C_offset = BUFFER_DIM * C;
            const Int L_offset = BUFFER_DIM * L;
            const Int R_offset = BUFFER_DIM * R;
            
            LOOP_UNROLL_FULL
            for( Int k = 0; k < BUFFER_DIM; ++k )
            {
                // Overwrite, not add-into. Thus cleansing is not needed.
                C_in[C_offset + k] = C_in[L_offset + k] + C_in[R_offset + k];
            }
        }
        
    } // PercolateUp_Tasks

    void percolateUp_Tasks_gen( const Int C, const Int free_thread_count ) const
    {
        // C = cluster index
        
        const Int L = C_left [C];
        const Int R = C_right[C];
        
        if( (L >= 0) && (R >= 0) )
        {
            // If not a leaf, compute the values of the children first.
            #pragma omp task final(free_thread_count<1)  shared( L )
                percolateUp_Tasks_gen( L, free_thread_count/2 );
            #pragma omp task final(free_thread_count<1)  shared( R )
                percolateUp_Tasks_gen( R, free_thread_count-free_thread_count/2 );
            #pragma omp taskwait
            
            // Aftwards, compute the sum of the two children.
            
            Real * restrict const c = C_in.data();
            
            const Int C_offset = buffer_dim * C;
            const Int L_offset = buffer_dim * L;
            const Int R_offset = buffer_dim * R;
            
            LOOP_UNROLL(4)
            for( Int k = 0; k < buffer_dim; ++k )
            {
                // Overwrite, not add-into. Thus cleansing is not needed.
                c[C_offset + k] = c[L_offset + k] + c[R_offset + k];
            }
        }
        
    } // PercolateUp_Tasks





    void PercolateDown_Tasks( const Int C ) const
    {
        ptic(ClassName()+"::PercolateDown_Tasks");
        switch ( buffer_dim )
        {
            case 1:
            {
                percolateDown_Tasks<1>( C, ThreadCount() );
                break;
            }
            case 2:
            {
                percolateDown_Tasks<2>( C, ThreadCount() );
                break;
            }
            case 3:
            {
                percolateDown_Tasks<3>( C, ThreadCount() );
                break;
            }
            case 4:
            {
                percolateDown_Tasks<4>( C, ThreadCount() );
                break;
            }
            case 5:
            {
                percolateDown_Tasks<5>( C, ThreadCount() );
                break;
            }
            case 6:
            {
                percolateDown_Tasks<6>( C, ThreadCount() );
                break;
            }
            case 7:
            {
                percolateDown_Tasks<7>( C, ThreadCount() );
                break;
            }
            case 8:
            {
                percolateDown_Tasks<8>( C, ThreadCount() );
                break;
            }
            case 9:
            {
                percolateDown_Tasks<9>( C, ThreadCount() );
                break;
            }
            case 10:
            {
                percolateDown_Tasks<10>( C, ThreadCount() );
                break;
            }
            case 11:
            {
                percolateDown_Tasks<11>( C, ThreadCount() );
                break;
            }
            case 12:
            {
                percolateDown_Tasks<12>( C, ThreadCount() );
                break;
            }
            case 13:
            {
                percolateDown_Tasks<13>( C, ThreadCount() );
                break;
            }
            case 14:
            {
                percolateDown_Tasks<14>( C, ThreadCount() );
                break;
            }
            case 15:
            {
                percolateDown_Tasks<15>( C, ThreadCount() );
                break;
            }
            case 16:
            {
                percolateDown_Tasks<16>( C, ThreadCount() );
                break;
            }
            case 17:
            {
                percolateDown_Tasks<17>( C, ThreadCount() );
                break;
            }
            case 18:
            {
                percolateDown_Tasks<18>( C, ThreadCount() );
                break;
            }
            case 19:
            {
                percolateDown_Tasks<19>( C, ThreadCount() );
                break;
            }
            case 20:
            {
                percolateDown_Tasks<20>( C, ThreadCount() );
                break;
            }
            default:
            {
                percolateDown_Tasks_gen( C, ThreadCount() );
            }
        }
        ptoc(ClassName()+"::PercolateDown_Tasks");
    }


    template<Int BUFFER_DIM>
    void percolateDown_Tasks(const Int C, const Int free_thread_count ) const
    {
        // C = cluster index
        
        const Int L = C_left [C];
        const Int R = C_right[C];
        
        if( ( L >= 0 ) && ( R >= 0 ) )
        {
            Real * restrict const c = C_out.data();
            
            const Int C_offset = BUFFER_DIM * C;
            const Int L_offset = BUFFER_DIM * L;
            const Int R_offset = BUFFER_DIM * R;
            
            LOOP_UNROLL_FULL
            for( Int k = 0; k < BUFFER_DIM; ++k )
            {
                const Real buffer = c[C_offset + k];
                c[L_offset + k] += buffer;
                c[R_offset + k] += buffer;
            }
            
            #pragma omp task final(free_thread_count<1)  shared( L )
                percolateDown_Tasks<BUFFER_DIM>( L, free_thread_count/2 );
            #pragma omp task final(free_thread_count<1)  shared( R )
                percolateDown_Tasks<BUFFER_DIM>( R, free_thread_count-free_thread_count/2 );
            #pragma omp taskwait
        }
    } // percolateDown_Tasks


    void percolateDown_Tasks_gen(const Int C, const Int free_thread_count ) const
    {
        // C = cluster index
        
        const Int L = C_left [C];
        const Int R = C_right[C];
        
        if( ( L >= 0 ) && ( R >= 0 ) )
        {
            Real * restrict const c = C_out.data();
            
            const Int C_offset = buffer_dim * C;
            const Int L_offset = buffer_dim * L;
            const Int R_offset = buffer_dim * R;
            
            LOOP_UNROLL(4)
            for( Int k = 0; k < buffer_dim; ++k )
            {
                const Real buffer = c[C_offset + k];
                c[L_offset + k] += buffer;
                c[R_offset + k] += buffer;
            }
            
            #pragma omp task final(free_thread_count<1)  shared( L )
                percolateDown_Tasks_gen( L, free_thread_count/2 );
            #pragma omp task final(free_thread_count<1)  shared( R )
                percolateDown_Tasks_gen( R, free_thread_count-free_thread_count/2 );
            #pragma omp taskwait
        }
    } // percolateDown_Tasks_gen
