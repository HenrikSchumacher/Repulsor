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
        
        if( (L >= null) && (R >= null) )
        {
            // If not a leaf, compute the values of the children first.
            percolateUp_Recursive<BUFFER_DIM>(L);
            percolateUp_Recursive<BUFFER_DIM>(R);
            
            // Aftwards, compute the sum of the two children.
            
            const Int C_offset = buffer_dim * C;
            const Int L_offset = buffer_dim * L;
            const Int R_offset = buffer_dim * R;
            
            if constexpr ( BUFFER_DIM > 0 )
            {
                LOOP_UNROLL_FULL
                for( Int k = 0; k < BUFFER_DIM; ++k )
                {
                    // Overwrite, not add-into. Thus cleansing is not required.
                    C_in[C_offset+k] = C_in[L_offset+k] + C_in[R_offset+k];
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
        
        if( (L >= null) && (R >= null) )
        {
            mut<Real> c = C_out.data();
            
            const Int C_offset = buffer_dim * C;
            const Int L_offset = buffer_dim * L;
            const Int R_offset = buffer_dim * R;
            
            if constexpr ( BUFFER_DIM > 0 )
            {
                LOOP_UNROLL_FULL
                for( Int k = 0; k < BUFFER_DIM; ++k )
                {
                    const Real buffer = c[C_offset+k];
                    c[L_offset + k] += buffer;
                    c[R_offset + k] += buffer;
                }
            }
            else
            {
                LOOP_UNROLL(4)
                for( Int k = 0; k < buffer_dim; ++k )
                {
                    const Real buffer = c[C_offset + k];
                    c[L_offset + k] += buffer;
                    c[R_offset + k] += buffer;
                }
            }
            
            
            percolateDown_Recursive<BUFFER_DIM>(L);
            percolateDown_Recursive<BUFFER_DIM>(R);
        }
    }; // percolateDown_Recursive
