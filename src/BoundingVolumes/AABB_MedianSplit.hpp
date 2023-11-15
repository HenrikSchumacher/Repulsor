#pragma once

#define BASE  AABB<AMB_DIM,Real,Int,SReal>
#define CLASS AABB_MedianSplit

namespace Repulsor
{
    
    
    // serialized_data is assumed to be an array of size SIZE. Will never be allocated by class! Instead, it is meant to be mapped onto an array of type SReal by calling the member SetPointer.
    
    // DATA LAYOUT
    // serialized_data[0] = squared radius
    // serialized_data[1],...,serialized_data[AMB_DIM] = center
    // serialized_data[AMB_DIM + 1],...,serialized_data[AMB_DIM + AMB_DIM] = half the edge lengths.
    
    template<int AMB_DIM, typename Real, typename Int, typename SReal>
    class CLASS : public BASE
    {
    public:
        
        CLASS() : BASE() {}

        // Copy constructor
        CLASS( const CLASS & other ) : BASE( other ) {}
        
        // Move constructor
        CLASS( CLASS && other ) noexcept : BASE( other ) {}
        
        // Copy assignment
        const CLASS & operator=( const CLASS & rhs)
        {
            return CLASS ( rhs );
        }
        
        // Move assignment
        const CLASS & operator=( CLASS && rhs)
        {
            return CLASS ( std::move(rhs) );
        }
        
        virtual ~CLASS() override = default;
        
    protected:
        
        using BASE::serialized_data;
//        using BASE::self_buffer;
        
    public:
        
        using BASE::Size;
        using BASE::SetPointer;
        using BASE::FromPrimitives;
        
    public:
        
        __ADD_CLONE_CODE__(CLASS)
        
    public:
        
        virtual Int Split(
            mref<PrimitiveSerialized<AMB_DIM,Real,Int,SReal>> P,           // primitive prototype; to be "mapped" over P_serialized, thus not const.
            mptr<SReal> P_serialized, const Int begin, const Int end,      // which _P_rimitives are in question
            mptr<Int>   P_ordering,                                        // to keep track of the permutation of the primitives
            mptr<SReal> C_data,       const Int C_ID,                      // where to get   the bounding volume info for _C_urrent bounding volume
            mptr<SReal> L_serialized, const Int L_ID,                      // where to store the bounding volume info for _L_eft  child (if successful!)
            mptr<SReal> R_serialized, const Int R_ID,                      // where to store the bounding volume info for _R_ight child (if successful!)
            mptr<SReal> score,                                             // some scratch buffer for one scalar per primitive
            mptr<Int>   perm,                                              // some scratch buffer for one Int per primitive (for storing local permutation)
            mptr<Int>   inv_perm,                                          // some scratch buffer for one Int per primitive (for storing inverse of local permutation)
            Int thread_count = 1                                           // how many threads to utilize
        ) override
        {
            return split(
                  P,P_serialized,begin,end,
                  P_ordering,
                  C_data,       C_ID,
                  L_serialized, L_ID,
                  R_serialized, R_ID,
                  score, perm, inv_perm, thread_count
            );
        }
        
        virtual Int Split(
            mref<PolytopeBase<AMB_DIM,Real,Int,SReal>> P,                  // primitive prototype; to be "mapped" over P_serialized, thus not const.
            mptr<SReal> P_serialized, const Int begin, const Int end,      // which _P_rimitives are in question
            mptr<Int>   P_ordering,                                        // to keep track of the permutation of the primitives
            mptr<SReal> C_data,       const Int C_ID,                      // where to get   the bounding volume info for _C_urrent bounding volume
            mptr<SReal> L_serialized, const Int L_ID,                      // where to store the bounding volume info for _L_eft  child (if successful!)
            mptr<SReal> R_serialized, const Int R_ID,                      // where to store the bounding volume info for _R_ight child (if successful!)
            mptr<SReal> score,                                             // some scratch buffer for one scalar per primitive
            mptr<Int>   perm,                                              // some scratch buffer for one Int per primitive (for storing local permutation)
            mptr<Int>   inv_perm,                                          // some scratch buffer for one Int per primitive (for storing inverse of local permutation)
            Int thread_count = 1                                           // how many threads to utilize
        ) override
        {
            return split(
                  P,P_serialized,begin,end,
                  P_ordering,
                  C_data,       C_ID,
                  L_serialized, L_ID,
                  R_serialized, R_ID,
                  score, perm, inv_perm, thread_count
            );
        }
        
    protected:
        
        template<typename Primitive_T>
        Int split(
            mref<Primitive_T> P,                                         // primitive prototype; to be "mapped" over P_serialized, thus not const.
            mptr<SReal> P_serialized, const Int begin, const Int end,    // which _P_rimitives are in question
            mptr<Int>   P_ordering,                                      // to keep track of the permutation of the primitives
            mptr<SReal> C_data,       const Int C_ID,                    // where to get   the bounding volume info for _C_urrent bounding volume
            mptr<SReal> L_serialized, const Int L_ID,                    // where to store the bounding volume info for _L_eft  child (if successful!)
            mptr<SReal> R_serialized, const Int R_ID,                    // where to store the bounding volume info for _R_ight child (if successful!)
            mptr<SReal> score,                                           // some scratch buffer for one scalar per primitive
            mptr<Int>   perm,                                            // some scratch buffer for one Int per primitive (for storing local permutation)
            mptr<Int>   inv_perm,                                        // some scratch buffer for one Int per primitive (for storing inverse of local permutation)
            Int thread_count                                             // how many threads to utilize
        )
        {
            Int P_Size = P.Size();
            
            SetPointer( C_data, C_ID );
            
            // Find the longest axis `split_dir` of primitives's bounding box.
            Int split_dir = 0;
            
            SReal L_max = serialized_data[1 + AMB_DIM + split_dir];
            
            for( Int k = 1; k < AMB_DIM; ++k )
            {
                SReal L_k = serialized_data[1 + AMB_DIM + k];
                if( L_k > L_max )
                {
                    L_max = L_k;
                    split_dir = k;
                }
            }
            
            if( L_max <= Scalar::Zero<SReal> )
            {
                eprint(ClassName()+"Split: longest axis has length <=0.");
                return -1;
            }
            
            
            // Finding the "median". Adapted from https://stackoverflow.com/a/16798127/8248900 (using pointers instead of iterators)

            // Computing score as pojection of the primitives' InteriorPoints on the longest axis.
            // Fill perm with the indices.

            cptr<SReal> p = P_serialized + 1 + split_dir;

            // TODO: Parallelize
            for( Int i = begin; i < end; ++i )
            {
                score[i] = p[ P_Size * i ];
                perm [i] = i;
            }

            Int split_index = begin + ((end-begin)/2);

            Int  * mid = perm + split_index;
            // TODO: Parallelize
            std::nth_element(
                    perm + begin, mid,
                    perm + end,
                    [score](const Int i, const Int j) {return score[i] < score[j];}
            );
            // Now perm contains the desired ordering of score.

            // Invert permutation.
            // TODO: Parallelize
            for( Int i = begin; i < end; ++i )
            {
                inv_perm[perm[i]] = i;
            }

            // https://www.geeksforgeeks.org/permute-the-elements-of-an-array-following-given-order/
            // Reorder primitive according to perm, i.e., write primitive perm[i] to position i.
            // TODO: Parallelize
            for( Int i = begin; i < end; ++i )
            {
                Int next = i;

                while( inv_perm[next] >= 0 )
                {
                    Int temp = inv_perm[next];
                    std::swap( P_ordering   [i]  , P_ordering   [temp] );
                    P.Swap   ( P_serialized, i,    P_serialized, temp  );
                    inv_perm[next] = -1;
                    next = temp;
                }
            }

            // Compute bounding volume of left child.
            SetPointer( L_serialized, L_ID );
            FromPrimitives( P, P_serialized, begin,       split_index,   thread_count );
            // Compute bounding volume of right child.
            SetPointer( R_serialized, R_ID );
            FromPrimitives( P, P_serialized, split_index, end,           thread_count );
            
            return split_index;
            
        } // Split
        
    public:
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+">";
        }
        
    }; // CLASS

} // namespace Repulsor

#undef CLASS
#undef BASE

