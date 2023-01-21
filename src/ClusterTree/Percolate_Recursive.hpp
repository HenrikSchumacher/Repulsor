protected:
    
    void PercolateUp_Recursive( const Int C ) const
    {
        ptic(ClassName()+"::PercolateUp_Recursive");
        switch ( buffer_dim )
        {
            case 1:
            {
                percolateUp_Recursive<1>(C);
                break;
            }
            case 2:
            {
                percolateUp_Recursive<2>(C);
                break;
            }
            case 3:
            {
                percolateUp_Recursive<3>(C);
                break;
            }
            case 4:
            {
                percolateUp_Recursive<4>(C);
                break;
            }
            case 5:
            {
                percolateUp_Recursive<5>(C);
                break;
            }
            case 6:
            {
                percolateUp_Recursive<6>(C);
                break;
            }
            case 7:
            {
                percolateUp_Recursive<7>(C);
                break;
            }
            case 8:
            {
                percolateUp_Recursive<8>(C);
                break;
            }
            case 9:
            {
                percolateUp_Recursive<9>(C);
                break;
            }
            case 10:
            {
                percolateUp_Recursive<10>(C);
                break;
            }
            case 11:
            {
                percolateUp_Recursive<11>(C);
                break;
            }
            case 12:
            {
                percolateUp_Recursive<12>(C);
                break;
            }
            case 13:
            {
                percolateUp_Recursive<13>(C);
                break;
            }
            case 14:
            {
                percolateUp_Recursive<14>(C);
                break;
            }
            case 15:
            {
                percolateUp_Recursive<15>(C);
                break;
            }
            case 16:
            {
                percolateUp_Recursive<16>(C);
                break;
            }
            case 17:
            {
                percolateUp_Recursive<17>(C);
                break;
            }
            case 18:
            {
                percolateUp_Recursive<18>(C);
                break;
            }
            case 19:
            {
                percolateUp_Recursive<19>(C);
                break;
            }
            case 20:
            {
                percolateUp_Recursive<20>(C);
                break;
            }
            default:
            {
                percolateUp_Recursive_gen(C);
            }
        }
        ptoc(ClassName()+"::PercolateUp_Recursive");
    }
    
    template<Int BUFFER_DIM>
    void percolateUp_Recursive( const Int C ) const
    {
        // C = cluster index
        
        const Int L = C_left [C];
        const Int R = C_right[C];
        
        if( (L >= null) && (R >= null) )
        {
            // If not a leaf, compute the values of the children first.
            percolateUp_Recursive<BUFFER_DIM>(L);
            percolateUp_Recursive<BUFFER_DIM>(R);
            
            // Aftwards, compute the sum of the two children.
//            Real * restrict const c = C_in.data();
            
            const Int C_offset = buffer_dim * C;
            const Int L_offset = buffer_dim * L;
            const Int R_offset = buffer_dim * R;
            
            LOOP_UNROLL_FULL
            for( Int k = 0; k < BUFFER_DIM; ++k )
            {
                // Overwrite, not add-into. Thus cleansing is not required.
                C_in[C_offset+k] = C_in[L_offset+k] + C_in[R_offset+k];
            }
        }
        
    }; // percolateUp_Recursive

    
    void percolateUp_Recursive_gen( const Int C ) const
    {
        // C = cluster index
        
        const Int L = C_left [C];
        const Int R = C_right[C];
        
        if( (L >= null) && (R >= null) )
        {
            // If not a leaf, compute the values of the children first.
            percolateUp_Recursive_gen(L);
            percolateUp_Recursive_gen(R);
            
            // Aftwards, compute the sum of the two children.
//            Real * restrict const c = C_in.data();
            
            const Int C_offset = buffer_dim * C;
            const Int L_offset = buffer_dim * L;
            const Int R_offset = buffer_dim * R;
            
            LOOP_UNROLL(4)
            for( Int k = 0; k < buffer_dim; ++k )
            {
                // Overwrite, not add-into. Thus cleansing is not required.
                C_in[C_offset + k] = C_in[L_offset + k] + C_in[R_offset + k];
            }
        }
        
    }; // percolateUp_Recursive_gen


protected:
    
    void PercolateDown_Recursive( const Int C ) const
    {
        ptic(ClassName()+"::PercolateDown_Recursive");
        switch ( buffer_dim )
        {
            case 1:
            {
                percolateDown_Recursive<1>(C);
                break;
            }
            case 2:
            {
                percolateDown_Recursive<2>(C);
                break;
            }
            case 3:
            {
                percolateDown_Recursive<3>(C);
                break;
            }
            case 4:
            {
                percolateDown_Recursive<4>(C);
                break;
            }
            case 5:
            {
                percolateDown_Recursive<5>(C);
                break;
            }
            case 6:
            {
                percolateDown_Recursive<6>(C);
                break;
            }
            case 7:
            {
                percolateDown_Recursive<7>(C);
                break;
            }
            case 8:
            {
                percolateDown_Recursive<8>(C);
                break;
            }
            case 9:
            {
                percolateDown_Recursive<9>(C);
                break;
            }
            case 10:
            {
                percolateDown_Recursive<10>(C);
                break;
            }
            case 11:
            {
                percolateDown_Recursive<11>(C);
                break;
            }
            case 12:
            {
                percolateDown_Recursive<12>(C);
                break;
            }
            case 13:
            {
                percolateDown_Recursive<13>(C);
                break;
            }
            case 14:
            {
                percolateDown_Recursive<14>(C);
                break;
            }
            case 15:
            {
                percolateDown_Recursive<15>(C);
                break;
            }
            case 16:
            {
                percolateDown_Recursive<16>(C);
                break;
            }
            case 17:
            {
                percolateDown_Recursive<17>(C);
                break;
            }
            case 18:
            {
                percolateDown_Recursive<18>(C);
                break;
            }
            case 19:
            {
                percolateDown_Recursive<19>(C);
                break;
            }
            case 20:
            {
                percolateDown_Recursive<20>(C);
                break;
            }
            default:
            {
                percolateDown_Recursive_gen(C);
            }
        }
        ptoc(ClassName()+"::PercolateDown_Recursive");
    }
    
    template<Int BUFFER_DIM>
    void percolateDown_Recursive( const Int C )  const
    {
        // C = cluster index
        
        const Int L = C_left [C];
        const Int R = C_right[C];
        
        if( ( L >= null ) && ( R >= null ) )
        {
            Real * restrict const c = C_out.data();
            
            const Int C_offset = BUFFER_DIM * C;
            const Int L_offset = BUFFER_DIM * L;
            const Int R_offset = BUFFER_DIM * R;
            
            LOOP_UNROLL_FULL
            for( Int k = 0; k < BUFFER_DIM; ++k )
            {
                const Real buffer = c[C_offset+k];
                c[L_offset + k] += buffer;
                c[R_offset + k] += buffer;
            }
            
            percolateDown_Recursive<BUFFER_DIM>(L);
            percolateDown_Recursive<BUFFER_DIM>(R);
        }
    }; // percolateDown_Recursive
    
    void percolateDown_Recursive_gen( const Int C )  const
    {
        // C = cluster index
        
        const Int L = C_left [C];
        const Int R = C_right[C];
        
        Real * restrict const c = C_out.data();
        
        if( ( L >= null ) && ( R >= null ) )
        {
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
            percolateDown_Recursive_gen(L);
            percolateDown_Recursive_gen(R);
        }
    }; // percolateDown_Recursive_gen
