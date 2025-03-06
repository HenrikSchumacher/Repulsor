#pragma once

namespace Repulsor
{
    // serialized_data is assumed to be an array of size SIZE. Will never be allocated by class! Instead, it is meant to be mapped onto an array of type SReal by calling the member SetPointer.
    
    // DATA LAYOUT
    // serialized_data[0] = squared radius
    // serialized_data[1],...,serialized_data[AMB_DIM] = center
    // serialized_data[AMB_DIM + 1],...,serialized_data[AMB_DIM + AMB_DIM] = half the edge lengths.
    
    template<int AMB_DIM, typename Real, typename Int, typename SReal>
    class AABB_PreorderedSplit : public AABB<AMB_DIM,Real,Int,SReal>
    {
    public:
        
        using Base_T = AABB<AMB_DIM,Real,Int,SReal>;
        
        AABB_PreorderedSplit() : Base_T() {}

        // Copy constructor
        AABB_PreorderedSplit( const AABB_PreorderedSplit & other ) : Base_T( other ) {}
        
        // Move constructor
        AABB_PreorderedSplit( AABB_PreorderedSplit && other ) noexcept : Base_T( other ) {}
        
        // Copy assignment
        const AABB_PreorderedSplit & operator=( const AABB_PreorderedSplit & rhs)
        {
            return AABB_PreorderedSplit ( rhs );
        }
        
        // Move assignment
        const AABB_PreorderedSplit & operator=( AABB_PreorderedSplit && rhs)
        {
            return AABB_PreorderedSplit ( std::move(rhs) );
        }
        
        virtual ~AABB_PreorderedSplit() override = default;
        
    protected:
        
        using Base_T::serialized_data;
//        using Base_T::self_buffer;
        
    public:
        
        using Base_T::Size;
        using Base_T::SetPointer;
        using Base_T::FromPrimitives;
        
    public:
        
        [[nodiscard]] std::shared_ptr<AABB_PreorderedSplit> Clone () const
        {
            return std::shared_ptr<AABB_PreorderedSplit>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual AABB_PreorderedSplit * CloneImplementation() const override
        {
            return new AABB_PreorderedSplit(*this);
        }
        
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
            mref<Primitive_T> P,                                        // primitive prototype; to be "mapped" over P_serialized, thus not const.
            mptr<SReal> P_serialized, const Int begin, const Int end,   // which _P_rimitives are in question
            mptr<Int>   P_ordering,                                     // to keep track of the permutation of the primitives
            mptr<SReal> C_data,       const Int C_ID,                   // where to get   the bounding volume info for _C_urrent bounding volume
            mptr<SReal> L_serialized, const Int L_ID,                   // where to store the bounding volume info for _L_eft  child (if successful!)
            mptr<SReal> R_serialized, const Int R_ID,                   // where to store the bounding volume info for _R_ight child (if successful!)
            mptr<SReal> score,                                          // some scratch buffer for one scalar per primitive
            mptr<Int>   perm,                                           // some scratch buffer for one Int per primitive (for storing local permutation)
            mptr<Int>   inv_perm,                                       // some scratch buffer for one Int per primitive (for storing inverse of local permutation)
            Int thread_count                                            // how many threads to utilize
        )
        {
            (void)P_ordering;
            (void)C_data;
            (void)C_ID;
            (void)score;
            (void)perm;
            (void)inv_perm;
            
            Int split_index = begin + ((end-begin)/2);
            
            // Compute bounding volume of left child.
            SetPointer( L_serialized, L_ID );
            FromPrimitives( P, P_serialized,   begin,       split_index, thread_count );
            // Compute bounding volume of right child.
            SetPointer( R_serialized, R_ID );
            FromPrimitives( P, P_serialized,   split_index, end,         thread_count );
            
            return split_index;
            
        } // split
        
    public:
        
        virtual std::string ClassName() const override
        {
            return std::string("AABB_PreorderedSplit")+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+">";
        }
        
        
    }; // AABB_PreorderedSplit

} // namespace Repulsor
