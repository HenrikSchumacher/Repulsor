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
               this->OffDiag().Dimension(0) != nnz
               ||
               this->OffDiag().Dimension(1) != BLOCK_NNZ
            )
            {
                this->OffDiag() = Values_T( nnz, BLOCK_NNZ );
            }
        }
    }
