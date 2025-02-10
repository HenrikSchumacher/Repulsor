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
        // The tree is a binary tree. 
        // Each node needs an identifier of type std::size_t.
        // With contemporary systems, std::size_t is 64 bit in size.
        // There is no way that such a tree can have 128 levels.
        // So it is relatively save here to use a stack of fixed size.
        
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
                combine_buffers3<Scalar::Flag::Plus,Scalar::Flag::Plus,BUFFER_DIM>(
                    Scalar::One<Real>, &C_in[buffer_dim * L],
                    Scalar::One<Real>, &C_in[buffer_dim * R],
                                       &C_in[buffer_dim * C],
                    buffer_dim
                );
                
                visited[stack_ptr] = false;
                --stack_ptr;  // pop
            }
        }
        
        if( stack_ptr >= 126 )
        {
            eprint(ClassName()+"::percolateUp_DFS: Stack overflow.");
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
        // The tree is a binary tree.
        // Each node needs an identifier of type std::size_t.
        // With contemporary systems, std::size_t is 64 bit in size.
        // There is no way that such a tree can have 128 levels.
        // So it is relatively save here to use a stack of fixed size.
        
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
                // If not a leaf, compute the values of the children first.
                combine_buffers<Scalar::Flag::Plus,Scalar::Flag::Plus,BUFFER_DIM>(
                    Scalar::One<Real>, &C_out[buffer_dim * C],
                    Scalar::One<Real>, &C_out[buffer_dim * L],
                    buffer_dim
                );
                
                combine_buffers<Scalar::Flag::Plus,Scalar::Flag::Plus,BUFFER_DIM>(
                    Scalar::One<Real>, &C_out[buffer_dim * C],
                    Scalar::One<Real>, &C_out[buffer_dim * R],
                    buffer_dim
                );
                
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
        
        if( stack_ptr >= 126 )
        {
            eprint(ClassName()+"::PercolateDown_DFS: Stack overflow.");
        }
    }; // PercolateDown_DFS
