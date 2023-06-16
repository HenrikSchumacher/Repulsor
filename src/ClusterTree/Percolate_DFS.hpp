protected:

    // TODO: Since the clusters ought to be in depth-first order, we do not need a stack here!

    // Sequential variant that uses a stack instead of recursion.
    void PercolateUp_DFS( const Int C, const Int max_depth = 64 ) const override
    {
        switch ( buffer_dim )
        {
            case 1:
            {
                percolateUp_DFS<1>(C,max_depth);
                break;
            }
            case AMB_DIM:
            {
                percolateUp_DFS<AMB_DIM>(C,max_depth);
                break;
            }
            case AMB_DIM * AMB_DIM:
            {
                percolateUp_DFS<AMB_DIM * AMB_DIM>(C,max_depth);
                break;
            }
            case (AMB_DIM + 1) * AMB_DIM:
            {
                percolateUp_DFS<(AMB_DIM + 1) * AMB_DIM>(C,max_depth);
                break;
            }
            case FAR_DIM:
            {
                percolateUp_DFS<FAR_DIM>(C,max_depth);
                break;
            }
            default:
            {
                percolateUp_DFS<0>(C,max_depth);
            }
        }
    }


    template<Int BUFFER_DIM>
    void percolateUp_DFS( const Int C_root, const Int max_depth ) const
    {
        Int stack   [128] = {};
        Int visited [128] = {false};
        Int depths  [128] = {};
        
        Int stack_ptr    = null;
        stack[stack_ptr] = C_root;
        
        while( (stack_ptr >= null) && ( stack_ptr < 126 ) )
        {
            // We are at cluster C.
            const Int d = depths[stack_ptr];
            const Int C = stack [stack_ptr];
            const Int L = C_left[C];
            const Int R = C_right[C];
            
            if( !visited[stack_ptr] )
            {
                if( (d < max_depth) && (L >= null) && (R >= null) )
                {
                    visited[stack_ptr] = true;
                    
                    // If not a leaf, compute the values of the children first.
                    
                    // push
                    ++stack_ptr;
                    stack [stack_ptr] = R;
                    depths[stack_ptr] = d+1;
                    
                    // push
                    ++stack_ptr;
                    stack [stack_ptr] = L;
                    depths[stack_ptr] = d+1;
                }
                else
                {
                    // If in leaf node, then backtrack.
                    --stack_ptr;
                }
            }
            else
            {
                const Int C_offset = buffer_dim * C;
                const Int L_offset = buffer_dim * L;
                const Int R_offset = buffer_dim * R;
                
                if constexpr ( BUFFER_DIM > 0 )
                {
                    LOOP_UNROLL_FULL
                    for( Int k = 0; k < BUFFER_DIM; ++k )
                    {
                        // Overwrite, not add-into. Thus cleansing is not required.
                        C_in[C_offset + k] = C_in[L_offset + k] + C_in[R_offset + k];
                    }
                }
                else
                {
                    LOOP_UNROLL(4)
                    for( Int k = 0; k < buffer_dim; ++k )
                    {
                        // Overwrite, not add-into. Thus cleansing is not required.
                        C_in[C_offset + k] = C_in[L_offset + k] + C_in[R_offset + k];
                    }
                }
                visited[stack_ptr] = false;
                --stack_ptr;  // pop
            }
        }
        
        if( stack_ptr >= 128 )
        {
            eprint(ClassName()+"::percolateUp_DFS: stack overflow detected.");
        }

    } // percolateUp_DFS

    
    // Sequential variant that uses a stack instead of recursion.
    void PercolateDown_DFS( const Int C, const Int max_depth = 64 ) const override
    {
//        ptic(ClassName()+"::PercolateDown_DFS");
        switch ( buffer_dim )
        {
            case 1:
            {
                percolateDown_DFS<1>(C,max_depth);
                break;
            }
            case AMB_DIM:
            {
                percolateDown_DFS<AMB_DIM>(C,max_depth);
                break;
            }
            case AMB_DIM * AMB_DIM:
            {
                percolateDown_DFS<AMB_DIM * AMB_DIM>(C,max_depth);
                break;
            }
            case (AMB_DIM + 1) * AMB_DIM:
            {
                percolateDown_DFS<(AMB_DIM + 1) * AMB_DIM>(C,max_depth);
                break;
            }
            case FAR_DIM:
            {
                percolateDown_DFS<FAR_DIM>(C,max_depth);
                break;
            }
            default:
            {
                percolateDown_DFS<0>(C,max_depth);
            }
        }
//        ptoc(ClassName()+"::PercolateDown_DFS");
    }

    template<Int BUFFER_DIM>
    void percolateDown_DFS( const Int C_root, const Int max_depth ) const
    {
        Int stack  [128] = {};
        Int depths [128] = {};

        Int stack_ptr    = null;
        stack[stack_ptr] = C_root;
        
        while( (stack_ptr >= null) && ( stack_ptr < 126 ) )
        {
            // We are at cluster C.
            const Int d = depths[stack_ptr];
            const Int C = stack [stack_ptr];
            const Int L = C_left[C];
            const Int R = C_right[C];
            --stack_ptr; // pop
            
            if( (d < max_depth) && (L >= null) && (R >= null) )
            {
                const Int C_offset = buffer_dim * C;
                const Int L_offset = buffer_dim * L;
                const Int R_offset = buffer_dim * R;
                
                // If not a leaf, compute the values of the children first.
                if constexpr ( BUFFER_DIM > 0 )
                {
                    LOOP_UNROLL_FULL
                    for( Int k = 0; k < BUFFER_DIM; ++k )
                    {
                        const Real buffer = C_out[C_offset + k];
                        
                        C_out[L_offset + k] += buffer;
                        C_out[R_offset + k] += buffer;
                    }
                }
                else
                {
                    LOOP_UNROLL(4)
                    for( Int k = 0; k < buffer_dim; ++k )
                    {
                        const Real buffer = C_out[C_offset + k];
                        
                        C_out[L_offset + k] += buffer;
                        C_out[R_offset + k] += buffer;
                    }
                }
                
                // push
                ++stack_ptr;
                stack [stack_ptr] = R;
                depths[stack_ptr] = d+1;
                
                // push
                ++stack_ptr;
                stack [stack_ptr] = L;
                depths[stack_ptr] = d+1;
            }
        }
        
        if( stack_ptr >= 128 )
        {
            eprint(ClassName()+"::PercolateDown_DFS: stack overflow detected.");
        }
    }; // PercolateDown_DFS
