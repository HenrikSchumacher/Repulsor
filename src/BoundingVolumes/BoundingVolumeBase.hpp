#pragma once

namespace Repulsor
{
    // Class that adds bounding volume features to a basic serializable primitive, like computing the bounding volume from a set of points.
    template<int AMB_DIM, typename Real, typename Int, typename SReal>
    class BoundingVolumeBase : public PrimitiveSerialized<AMB_DIM,Real,Int,SReal>
    {
    public:
        
        using Base_T = PrimitiveSerialized<AMB_DIM,Real,Int,SReal>;
        
        // Default constructor
        BoundingVolumeBase()
        : Base_T()
        {}
        // Destructor
        virtual ~BoundingVolumeBase() override = default;
        // Copy constructor
        BoundingVolumeBase( const BoundingVolumeBase & other ) = default;
        // Copy assignment operator
        BoundingVolumeBase & operator=( const BoundingVolumeBase & other ) = default;
        // Move constructor
        BoundingVolumeBase( BoundingVolumeBase && other ) = default;
        // Move assignment operator
        BoundingVolumeBase & operator=( BoundingVolumeBase && other ) = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<BoundingVolumeBase> Clone () const
        {
            return std::shared_ptr<BoundingVolumeBase>(CloneImplementation());
        }
        
    private:
        
        [[nodiscard]] virtual BoundingVolumeBase * CloneImplementation() const override = 0;
        
    public:
        
        // array p is suppose to represent a matrix of size N x AMB_DIM
        virtual void FromPointCloud( cptr<SReal> coords, const Int N ) const = 0;
        
        virtual void FromPrimitives(
            mref<PrimitiveSerialized<AMB_DIM,Real,Int,SReal>> P,
            mptr<SReal> P_data,
            const Int begin,
            const Int end,           // which _P_primitives are in question
            Int thread_count = 1     // how many threads to utilize
        ) const = 0;
        
        virtual Int Split(
            mref<PrimitiveSerialized<AMB_DIM,Real,Int,SReal>> P,   // primitive prototype; to be "mapped" over P_data, thus not const.
            mptr<SReal> P_data, const Int begin, const Int end,    // which _P_rimitives are in question
            mptr<Int>   P_ordering,                                // to keep track of the permutation of the primitives
            mptr<SReal> C_data, const Int C_ID,                    // where to get   the bounding volume info for _C_urrent bounding volume
            mptr<SReal> L_data, const Int L_ID,                    // where to store the bounding volume info for _L_eft  child (if successful!)
            mptr<SReal> R_data, const Int R_ID,                    // where to store the bounding volume info for _R_ight child (if successful!)
            mptr<SReal> score,                                     // some scratch buffer for one scalar per primitive
            mptr<Int>   perm,                                      // some scratch buffer for one Int per primitive (for storing local permutation)
            mptr<Int>   inv_perm,                                  // some scratch buffer for one Int per primitive (for storing inverse of local permutation)
            Int thread_count = 1                                   // how many threads to utilize
        ) = 0;
        
        virtual Int Split(
            mref<PolytopeBase<AMB_DIM,Real,Int,SReal>> P,          // primitive prototype; to be "mapped" over P_data, thus not const.
            mptr<SReal> P_data, const Int begin, const Int end,    // which _P_rimitives are in question
            mptr<Int>   P_ordering,                                // to keep track of the permutation of the primitives
            mptr<SReal> C_data, const Int C_ID,                    // where to get   the bounding volume info for _C_urrent bounding volume
            mptr<SReal> L_data, const Int L_ID,                    // where to store the bounding volume info for _L_eft  child (if successful!)
            mptr<SReal> R_data, const Int R_ID,                    // where to store the bounding volume info for _R_ight child (if successful!)
            mptr<SReal> score,                                     // some scratch buffer for one scalar per primitive
            mptr<Int>   perm,                                      // some scratch buffer for one Int per primitive (for storing local permutation)
            mptr<Int>   inv_perm,                                  // some scratch buffer for one Int per primitive (for storing inverse of local permutation)
            Int thread_count = 1                                   // how many threads to utilize
        ) = 0;

        virtual std::string ClassName() const override
        {
            return std::string("BoundingVolumeBase")+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+">";
        }
        
    }; // BoundingVolumeBase

} // namespace Repulsor
