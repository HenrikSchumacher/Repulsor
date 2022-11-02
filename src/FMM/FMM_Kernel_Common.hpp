public:

    static constexpr Int NonzeroCount()
    {
        return BLOCK_NNZ;
    }

    void Allocate( const LInt nnz )
    {
        if constexpr ( metric_flag )
        {
            if(
               this->metric_values.Dimension(0) != nnz
               ||
               this->metric_values.Dimension(1) != BLOCK_NNZ
               )
            {
                this->metric_values = Values_T( nnz, BLOCK_NNZ );
            }
        }
    }
