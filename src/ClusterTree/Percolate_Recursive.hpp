protected:
    
    void PercolateUp_Recursive( const Int C ) const override
    {
        ptic(ClassName()+"::PercolateUp_Recursive");
        switch ( buffer_dim )
        {
            case 1:
            {
                percolateUp_Recursive<1>(C);
                break;
            }
            case AMB_DIM:
            {
                percolateUp_Recursive<AMB_DIM>(C);
                break;
            }
            case AMB_DIM * AMB_DIM:
            {
                percolateUp_Recursive<AMB_DIM * AMB_DIM>(C);
                break;
            }
            case (AMB_DIM + 1) * AMB_DIM:
            {
                percolateUp_Recursive<(AMB_DIM + 1) * AMB_DIM>(C);
                break;
            }
            case FAR_DIM:
            {
                percolateUp_Recursive<FAR_DIM>(C);
                break;
            }
            default:
            {
                percolateUp_Recursive<0>(C);
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
        
        // Using that children are either both interiors or both leaves.
        if( (L >= null) /*&& (R >= null)*/ )
        {
            // If not a leaf, compute the values of the children first.
            percolateUp_Recursive<BUFFER_DIM>(L);
            percolateUp_Recursive<BUFFER_DIM>(R);
            
            // Aftwards, compute the sum of the two children.
            combine_buffers3<Scalar::Flag::Plus,Scalar::Flag::Plus,BUFFER_DIM>(
                Scalar::One<Real>, &C_in[buffer_dim * L],
                Scalar::One<Real>, &C_in[buffer_dim * R],
                                   &C_in[buffer_dim * C],
                buffer_dim
            );
        }
        
    }; // percolateUp_Recursive

    

protected:
    
    void PercolateDown_Recursive( const Int C ) const override
    {
        ptic(ClassName()+"::PercolateDown_Recursive");
        switch ( buffer_dim )
        {
            case 1:
            {
                percolateDown_Recursive<1>(C);
                break;
            }
            case AMB_DIM:
            {
                percolateDown_Recursive<AMB_DIM>(C);
                break;
            }
            case AMB_DIM * AMB_DIM:
            {
                percolateDown_Recursive<AMB_DIM * AMB_DIM>(C);
                break;
            }
            case (AMB_DIM + 1) * AMB_DIM:
            {
                percolateDown_Recursive<(AMB_DIM + 1) * AMB_DIM>(C);
                break;
            }
            case FAR_DIM:
            {
                percolateDown_Recursive<FAR_DIM>(C);
                break;
            }
            default:
            {
                percolateDown_Recursive<0>(C);
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
        
        // Using that children are either both interiors or both leaves.
        if( (L >= null) /*&& (R >= null)*/ )
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
            
            percolateDown_Recursive<BUFFER_DIM>(L);
            percolateDown_Recursive<BUFFER_DIM>(R);
        }
    }; // percolateDown_Recursive
