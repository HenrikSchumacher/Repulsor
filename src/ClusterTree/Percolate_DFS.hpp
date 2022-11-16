protected:
    
    // Sequential variant that uses a stack instead of recursion.
    void PercolateUp_DFS( const Int C ) const
    {
        ptic(ClassName()+"::PercolateUp_DFS");
        switch ( buffer_dim )
        {
            case 1:
            {
                percolateUp_DFS<1>(C);
                break;
            }
            case 2:
            {
                percolateUp_DFS<2>(C);
                break;
            }
            case 3:
            {
                percolateUp_DFS<3>(C);
                break;
            }
            case 4:
            {
                percolateUp_DFS<4>(C);
                break;
            }
            case 5:
            {
                percolateUp_DFS<5>(C);
                break;
            }
            case 6:
            {
                percolateUp_DFS<6>(C);
                break;
            }
            case 7:
            {
                percolateUp_DFS<7>(C);
                break;
            }
            case 8:
            {
                percolateUp_DFS<8>(C);
                break;
            }
            case 9:
            {
                percolateUp_DFS<9>(C);
                break;
            }
            case 10:
            {
                percolateUp_DFS<10>(C);
                break;
            }
            case 11:
            {
                percolateUp_DFS<11>(C);
                break;
            }
            case 12:
            {
                percolateUp_DFS<12>(C);
                break;
            }
            case 13:
            {
                percolateUp_DFS<13>(C);
                break;
            }
            case 14:
            {
                percolateUp_DFS<14>(C);
                break;
            }
            case 15:
            {
                percolateUp_DFS<15>(C);
                break;
            }
            case 16:
            {
                percolateUp_DFS<16>(C);
                break;
            }
            case 17:
            {
                percolateUp_DFS<17>(C);
                break;
            }
            case 18:
            {
                percolateUp_DFS<18>(C);
                break;
            }
            case 19:
            {
                percolateUp_DFS<19>(C);
                break;
            }
            case 20:
            {
                percolateUp_DFS<20>(C);
                break;
            }
            default:
            {
                percolateUp_DFS_gen(C);
            }
        }
        ptoc(ClassName()+"::PercolateUp_DFS");
    }


    template<Int BUFFER_DIM>
    void percolateUp_DFS( const Int C_root ) const
    {
        Int stack   [128] = {};
        Int visited [128] = {false};
        
        Int stack_ptr          = static_cast<Int>(0);
        stack[stack_ptr] = C_root;
        
        const Int  * restrict const left  = C_left.data();
        const Int  * restrict const right = C_right.data();
              Real * restrict const c     = C_in.data();
        
        while( (stack_ptr >= 0) && ( stack_ptr < 126 ) )
        {
            // We are at cluster C.
            const Int C = stack[stack_ptr];
            const Int L = left [C];
            const Int R = right[C];
            
            if( !visited[stack_ptr] )
            {
                visited[stack_ptr] = true;
                
                if( (L >= 0) && (R >= 0) )
                {
                    // If not a leaf, compute the values of the children first.
                    stack[++stack_ptr] = R; // push
                    stack[++stack_ptr] = L; // push
                }
            }
            else
            {
                const Int C_offset = BUFFER_DIM * C;
                const Int L_offset = BUFFER_DIM * L;
                const Int R_offset = BUFFER_DIM * R;
                
                #pragma clang loop unroll(enable)
                for( Int k = 0; k < BUFFER_DIM; ++k )
                {
                    // Overwrite, not add-into. Thus cleansing is not required.
                    c[C_offset + k] = c[L_offset + k] + c[R_offset + k];
                }
                visited[stack_ptr] = false;
                stack_ptr--;  // pop
            }
        }
        
        if( stack_ptr >= 128 )
        {
            eprint(ClassName()+"::percolateUp_DFS: stack overflow detected.");
        }

    } // percolateUp_DFS_gen


    void percolateUp_DFS_gen( const Int C_root ) const
    {
        Int stack   [128] = {};
        Int visited [128] = {false};
        
        Int stack_ptr          = static_cast<Int>(0);
        stack[stack_ptr] = C_root;
        
        const Int  * restrict const left  = C_left.data();
        const Int  * restrict const right = C_right.data();
              Real * restrict const c     = C_in.data();
        
        while( (stack_ptr >= 0) && ( stack_ptr < 126 ) )
        {
            // We are at cluster C.
            const Int C = stack[stack_ptr];
            const Int L = left [C];
            const Int R = right[C];
            
            if( !visited[stack_ptr] )
            {
                visited[stack_ptr] = true;
                
                if( (L >= 0) && (R >= 0) )
                {
                    // If not a leaf, compute the values of the children first.
                    stack[++stack_ptr] = R; // push
                    stack[++stack_ptr] = L; // push
                }
            }
            else
            {
                const Int C_offset = buffer_dim * C;
                const Int L_offset = buffer_dim * L;
                const Int R_offset = buffer_dim * R;
                
                #pragma clang loop unroll_count(4)
                for( Int k = 0; k < buffer_dim; ++k )
                {
                    // Overwrite, not add-into. Thus cleansing is not required.
                    c[C_offset + k] = c[L_offset + k] + c[R_offset + k];
                }
                visited[stack_ptr] = false;
                stack_ptr--;  // pop
            }
        }
        
        if( stack_ptr >= 128 )
        {
            eprint(ClassName()+"::percolateUp_DFS_gen: stack overflow detected.");
        }

    } // percolateUp_DFS_gen






    
    // Sequential variant that uses a stack instead of recursion.
    void PercolateDown_DFS( const Int C ) const
    {
        ptic(ClassName()+"::PercolateDown_DFS");
        switch ( buffer_dim )
        {
            case 1:
            {
                percolateDown_DFS<1>(C);
                break;
            }
            case 2:
            {
                percolateDown_DFS<2>(C);
                break;
            }
            case 3:
            {
                percolateDown_DFS<3>(C);
                break;
            }
            case 4:
            {
                percolateDown_DFS<4>(C);
                break;
            }
            case 5:
            {
                percolateDown_DFS<5>(C);
                break;
            }
            case 6:
            {
                percolateDown_DFS<6>(C);
                break;
            }
            case 7:
            {
                percolateDown_DFS<7>(C);
                break;
            }
            case 8:
            {
                percolateDown_DFS<8>(C);
                break;
            }
            case 9:
            {
                percolateDown_DFS<9>(C);
                break;
            }
            case 10:
            {
                percolateDown_DFS<10>(C);
                break;
            }
            case 11:
            {
                percolateDown_DFS<11>(C);
                break;
            }
            case 12:
            {
                percolateDown_DFS<12>(C);
                break;
            }
            case 13:
            {
                percolateDown_DFS<13>(C);
                break;
            }
            case 14:
            {
                percolateDown_DFS<14>(C);
                break;
            }
            case 15:
            {
                percolateDown_DFS<15>(C);
                break;
            }
            case 16:
            {
                percolateDown_DFS<16>(C);
                break;
            }
            case 17:
            {
                percolateDown_DFS<17>(C);
                break;
            }
            case 18:
            {
                percolateDown_DFS<18>(C);
                break;
            }
            case 19:
            {
                percolateDown_DFS<19>(C);
                break;
            }
            case 20:
            {
                percolateDown_DFS<20>(C);
                break;
            }
            default:
            {
                percolateDown_DFS_gen(C);
            }
        }
        ptoc(ClassName()+"::PercolateDown_DFS");
    }

    template<Int BUFFER_DIM>
    void percolateDown_DFS( const Int C_root ) const
    {
        Int stack [128] = {};

        Int stack_ptr = static_cast<Int>(0);
        stack[stack_ptr] = C_root;
        
        const Int  * restrict const left  = C_left.data();
        const Int  * restrict const right = C_right.data();
              Real * restrict const c     = C_out.data();
        
        while( (stack_ptr >= 0) && ( stack_ptr < 126 ) )
        {
            // We are at cluster C.
            const Int C = stack[stack_ptr--]; //pop
            const Int L = left [C];
            const Int R = right[C];
            
            if( (L >= 0) && (R >= 0) )
            {
                const Int C_offset = BUFFER_DIM * C;
                const Int L_offset = BUFFER_DIM * L;
                const Int R_offset = BUFFER_DIM * R;
                
                // If not a leaf, compute the values of the children first.
                #pragma clang loop unroll(enable)
                for( Int k = 0; k < BUFFER_DIM; ++k )
                {
                    const Real buffer = c[C_offset + k];
                    
                    c[L_offset + k] += buffer;
                    c[R_offset + k] += buffer;
                }
                
                stack[++stack_ptr] = R; // push
                stack[++stack_ptr] = L; // push
            }
        }
        
        if( stack_ptr >= 128 )
        {
            eprint(ClassName()+"::PercolateDown_DFS: stack overflow detected.");
        }
    }; // PercolateDown_DFS


    void percolateDown_DFS_gen( const Int C_root ) const
    {
        Int stack [128] = {};

        Int stack_ptr = static_cast<Int>(0);
        stack[stack_ptr] = C_root;
        
        const Int  * restrict const left  = C_left.data();
        const Int  * restrict const right = C_right.data();
              Real * restrict const c     = C_out.data();
        
        while( (stack_ptr >= 0) && ( stack_ptr < 126 ) )
        {
            // We are at cluster C.
            const Int C = stack[stack_ptr--]; //pop
            const Int L = left [C];
            const Int R = right[C];
            
            if( (L >= 0) && (R >= 0) )
            {
                const Int C_offset = buffer_dim * C;
                const Int L_offset = buffer_dim * L;
                const Int R_offset = buffer_dim * R;
                
                // If not a leaf, compute the values of the children first.
                #pragma clang loop unroll_count(4)
                for( Int k = 0; k < buffer_dim; ++k )
                {
                    const Real buffer = c[C_offset + k];
                    
                    c[L_offset + k] += buffer;
                    c[R_offset + k] += buffer;
                }
                
                stack[++stack_ptr] = R; // push
                stack[++stack_ptr] = L; // push
            }
        }
        
        if( stack_ptr >= 128 )
        {
            eprint(ClassName()+"::PercolateDown_DFS_gen: stack overflow detected.");
        }
    } // PercolateDown_DFS_gen
