#pragma once

namespace Repulsor
{
    // Adds some common I/O interface for all Polytope types to PrimitiveSerialized.

    template<int AMB_DIM,typename Real,typename Int,typename SReal>
    class PolytopeBase : public PrimitiveSerialized<AMB_DIM,Real,Int,SReal>
    {
        static_assert(FloatQ<SReal>,"");
        
    protected:
        
        // serialized_data is assumed to be an array of size SIZE. Will never be allocated by class! Instead, it is meant to be mapped onto an array of type Real by calling the member SetPointer.
        
        // DATA LAYOUT
        // serialized_data[0] = squared radius
        // serialized_data[1],...,serialized_data[AMB_DIM] = interior_point
        // serialized_data[1+ AMB_DIM],...,serialized_data[SIZE] = data that defines the primitive.
        
    public:
        
        using Base_T = PrimitiveSerialized<AMB_DIM,Real,Int,SReal>;
        
        PolytopeBase() {}
        
        // Copy constructor
        PolytopeBase( const PolytopeBase & other ) : Base_T()
        {
            this->serialized_data = other.serialized_data;
        }

        // Move constructor
        PolytopeBase( PolytopeBase && other ) noexcept : Base_T()
        {
            this->serialized_data = std::move(other.serialized_data);
        }
        
        virtual ~PolytopeBase() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<PolytopeBase> Clone () const
        {
            return std::shared_ptr<PolytopeBase>(CloneImplementation());
        }
        
    private:
        
        [[nodiscard]] virtual PolytopeBase * CloneImplementation() const override = 0;
        
    public:
        
        virtual Int PointCount() const = 0;
        
//        virtual void FromCoordinates( const SReal * const hull_coords_, const Int pos ) const = 0;
//
//        virtual void FromCoordinates( const SReal * const hull_coords_ ) const = 0;
//
//        virtual void FromDoubleCoordinates( const double * const hull_coords_, Int pos ) const = 0;
//
//        virtual void FromDoubleCoordinates( const double * const hull_coords_ ) const = 0;
//
//        virtual void FromIndexList( const SReal * const coords_, const Int * const tuples, const Int pos ) const = 0;
//
//        virtual void FromIndexList( const SReal * const coords_, const Int * const tuple ) const = 0;
//
//        virtual void FromIndexList( const double * const coords_, const Int * const tuples, const Int pos ) const = 0;
//
//        virtual void FromIndexList( const double * const coords_, const Int * const tuple ) const = 0;
        
        
        // Helper function to compute axis-aligned bounding boxes. in the format of box_min, box_max vector.
        // box_min, box_max are supposed to be vectors of size AMB_DIM.
        // BoxMinMax computes the "lower left" lo and "upper right" hi vectors of the primitives bounding box and sets box_min = min(lo, box_min) and box_max = min(h, box_max)
        virtual void BoxMinMax( mptr<SReal> box_min, mptr<SReal> box_max ) const  = 0;

//
//        friend void swap(PolytopeBase &A, PolytopeBase &B) noexcept
//        {
//            // see https://stackoverflow.com/questions/5695548/public-friend-swap-member-function for details
//            using std::swap;
//
//            swap( A.serialized_data, B.serialized_data );
//        }

//        // copy-and-swap idiom
//        virtual PolytopeBase & operator=(PolytopeBase other)
//        {
//            // see https://stackoverflow.com/a/3279550/8248900 for details
//
//            swap(*this, other);
//            return *this;
//
//        }
        
        virtual std::string ClassName() const override
        {
            return std::string("PolytopeBase")+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+">";
        }
        
    }; // PolytopeBase

    
} // namespace Repulsor
