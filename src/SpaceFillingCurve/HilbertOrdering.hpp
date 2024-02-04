
public:


//#######################################################################################
//####                        HilbertOrdering (hybrid sort)                         #####
//#######################################################################################


    template<typename Real, typename ExtInt>
    void HilbertOrdering( cptr<Real> X, mptr<ExtInt> perm, const Int point_count )
    {
        ptic(ClassName()+"::HilbertOrdering");

        Tensor1<Hilbert_T,Int> M = CoordToHilbert( X, point_count );
        
        ptic("Ordering");
        Ordering( M.data(), perm, point_count );
        ptoc("Ordering");
        
        ptoc(ClassName()+"::HilbertOrdering");
    }


    template<typename Real>
    Tensor1<Int,Int> HilbertOrdering( cptr<Real> X, const Int point_count )
    {
        Tensor1<Int,Int> perm ( point_count );
        
        HilbertOrdering( X, perm.data(), point_count );
        
        return perm;
    }
