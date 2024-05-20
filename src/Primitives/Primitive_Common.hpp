public:

    using Base_T::serialized_data;


public:

    // Sets the classe's data pointer.
    // We assume that p_ is an array of sufficient size in which the primitive's data is found between
    //      begin = p_ + SIZE * pos
    // and
    //      end   = p_ + SIZE * (pos+1).

    virtual void SetPointer( mptr<SReal> p_ ) override
    {
        serialized_data = p_;
    }

    virtual void SetPointer( mptr<SReal> p_, const Int pos ) override
    {
        SetPointer(&p_[SIZE * pos]);
    }

    virtual void Read( cptr<SReal> p_in, const Int i ) const override
    {
        Read( &p_in[SIZE * i] );
    }

    virtual void Read( cptr<SReal> p_in ) const override
    {
        copy_buffer<SIZE>( p_in, serialized_data );
    }


    virtual void Write( mptr<SReal> q_out, const Int j ) const override
    {
        Write( &q_out[SIZE * j] );
    }

    virtual void Write( mptr<SReal> q_out ) const override
    {
        copy_buffer<SIZE>( serialized_data, q_out );
    }


    virtual void Swap( mptr<SReal> p_out, const Int i, mptr<SReal> q_out, const Int j ) const override
    {
        Swap( &p_out[SIZE * i], &q_out[SIZE * j] );
    }

    virtual void Swap( mptr<SReal> p, mptr<SReal> q ) const override
    {
        std::swap_ranges( q, &q[SIZE], p );
    }

    virtual std::string DataString() const override
    {
        std::stringstream s;
        
        cptr<SReal> a = serialized_data;
        
        s << ClassName();
        s << ": data = { " << a[0];
        
        for( Int k = 1; k < SIZE; ++k )
        {
            s << ", " << a[k];
        }
        s << " }";
        return s.str();
    }
