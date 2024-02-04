public:

    template<typename Real, typename ExtInt>
    void MortonOrdering( cptr<Real> X, mptr<ExtInt> perm, const Int point_count )
    {
        ptic(ClassName()+"::MortonOrdering");

        Tensor1<MortonCode_T,Int> M = CoordsToMorton( X, point_count );
        
        ptic("Ordering");
        Ordering( M.data(), perm, point_count );
        ptoc("Ordering");
        
        ptoc(ClassName()+"::MortonOrdering");
    }


    template<typename Real>
    Tensor1<Int,Int> MortonOrdering( cptr<Real> X, const Int point_count )
    {
        Tensor1<Int,Int> perm ( point_count );
        
        MortonOrdering( X, perm.data(), point_count );
        
        return perm;
    }
