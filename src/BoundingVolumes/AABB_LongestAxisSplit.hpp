#pragma once

#define BASE  AABB<AMB_DIM,Real,Int,SReal>
#define CLASS AABB_LongestAxisSplit

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
        
    public:
        
        using BASE::Size;
        using BASE::SetPointer;
        using BASE::FromPrimitives;
        
    public:
        
        __ADD_CLONE_CODE__(CLASS)
    
    public:
        
        virtual Int Split(
            PrimitiveSerialized<AMB_DIM,Real,Int,SReal> & P,               // primitive prototype; to be "mapped" over P_serialized, thus not const.
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
//            ptic(ClassName()+"::Split (PrimitiveSerialized)");
            
            const Int P_Size = P.Size();
            
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
            
//            if( L_max <= Scalar::Zero<SReal> )
//            {
//                eprint(ClassName()+"Split: longest axis has length <=0.");
//                DUMP(this->DataString());
//                //                ptoc(ClassName()+"::Split (PrimitiveSerialized)");
//                return -1;
//            }
            
            // The interesting coordinate of center position.
            const SReal mid = serialized_data[1 + split_dir];
            
            Int split_index = begin;
            
            // Swap primitives according to their interior points in direction split_dir;
            mptr<SReal> p = &P_serialized[1 + split_dir];
            for( Int i = begin; i < end; ++i )
            {
                //                P.SetPointer( P_serialized, i );
                //
                //                SReal x = P.InteriorPoint( split_dir );
                
                // WARNING: For performance reasons, we do NOT use P.InteriorPoint for reading out the coordinate x!
                const SReal x = p[ P_Size * i ];
                
                if( x < mid )
                {
                    std::swap( P_ordering[split_index], P_ordering[i] );
                    
                    P.Swap( P_serialized, split_index, P_serialized, i );
                    
                    ++split_index;
                }
            }
            
            // It can happen then the InteriorPoints of all primitives lie on one side of mid because mid is the center of the bounding box of the primitives --- and not the center of the bounding box of the InteriorPoints! This guarantees that at least one point is split off the cluster.

//            split_index = Min( Max( begin + 1, split_index), end - 1 );
            
            if( (split_index == begin) || (split_index == end) )
            {
                // It can happen then the InteriorPoints of all primitives lie on one side of mid because mid is the center of the bounding box of the primitives --- and not the center of the bounding box of the InteriorPoints! We use median split as fallback.
                
                split_index = begin + ((end-begin)/2);
                
                for( Int i = begin; i < end; ++i )
                {
                    score[i] = p[ P_Size * i ];
                    perm [i] = i;
                }

                Int  * mid_ = perm + split_index;
                std::nth_element( perm + begin, mid_, perm + end,
                        [score](const Int i, const Int j) {return score[i] < score[j];}
                );
                // Now perm contains the desired ordering of score.

                // Invert permutation.
                for( Int i = begin; i < end; ++i )
                {
                    inv_perm[perm[i]] = i;
                }

                // https://www.geeksforgeeks.org/permute-the-elements-of-an-array-following-given-order/
                // Reorder primitive according to perm, i.e., write primitive perm[i] to position i.
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
            }

            // Compute bounding volume of left child.
            SetPointer( L_serialized, L_ID );
            FromPrimitives( P, P_serialized, begin,       split_index,   thread_count );
            // Compute bounding volume of right child.
            SetPointer( R_serialized, R_ID );
            FromPrimitives( P, P_serialized, split_index, end,           thread_count );
            
//            ptoc(ClassName()+"::Split (PrimitiveSerialized)");
            
            return split_index;
        } // Split
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+">";
        }
        
        
    }; // CLASS

} // namespace Repulsor

#undef CLASS
#undef BASE

