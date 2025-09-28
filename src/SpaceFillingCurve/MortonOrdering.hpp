public:

    template<typename Real, typename ExtInt>
    void MortonOrdering( cptr<Real> X, mptr<ExtInt> perm, const Int point_count )
    {
        TOOLS_PTIMER(timer,ClassName()+"::MortonOrdering");
        Tensor1<MortonCode_T,Int> M = CoordsToMorton( X, point_count );
        Ordering( M.data(), perm, point_count );
    }

    template<typename Real>
    Tensor1<Int,Int> MortonOrdering( cptr<Real> X, const Int point_count )
    {
        Tensor1<Int,Int> perm ( point_count );
        MortonOrdering( X, perm.data(), point_count );
        return perm;
    }
