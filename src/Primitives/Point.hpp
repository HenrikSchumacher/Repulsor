#pragma once

namespace Repulsor
{
    template<int AMB_DIM,typename Real,typename Int,typename SReal,
                typename ExtReal = SReal,typename ExtInt = Int>
    class Point : public PrimitiveSerialized<AMB_DIM,Real,Int,SReal>
    {
    public:
        
        static_assert(FloatQ<ExtReal>,"");
        
        static_assert(IntQ<ExtInt>,"");
        
        using Base_T = PrimitiveSerialized<AMB_DIM,Real,Int,SReal>;
        
    protected:
        
        // serialized_data is assumed to be an array of size SIZE. Will never be allocated by class! Instead, it is meant to be mapped onto an array of type Real by calling the member SetPointer.
        
        // DATA LAYOUT
        // serialized_data[0] = squared radius = 0
        // serialized_data[1],...,serialized_data[AMB_DIM] = interior_point
        
    public:
        
        Point() : Base_T() {}
        
        // Copy constructor
        Point( const Point & other ) : Base_T( other ) {}

        // Move constructor
        Point( Point && other ) noexcept : Base_T( other ) {}

        virtual ~Point() override = default;
        
        static constexpr Int SIZE = 1 + AMB_DIM;
        
        virtual constexpr Int Size() const override
        {
            return SIZE;
        }

//    protected:
//
//        mutable SReal self_buffer [SIZE];
        
#include "Primitive_Common.hpp"
        
    public:
        
        [[nodiscard]] std::shared_ptr<Point> Clone () const
        {
            return std::shared_ptr<Point>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual Point * CloneImplementation() const override
        {
            return new Point(*this);
        }
        
    public:
        
        virtual void FromCoordinates( cptr<ExtReal> coords, const Int i = 0 ) const
        {
            serialized_data[0] = Scalar::Zero<SReal>;
            
            copy_buffer<AMB_DIM>( coords, &serialized_data[1] );
        }
        
        
        //Computes support vector supp of dir.
        virtual Real MinSupportVector( cptr<Real> dir, mptr<Real> supp ) const override
        {
            copy_buffer<AMB_DIM>( &serialized_data[1], supp );
            
            return dot_buffer<AMB_DIM>( dir, supp );
        }
        
        //Computes support vector supp of dir.
        virtual Real MaxSupportVector( cptr<Real> dir, mptr<Real> supp ) const override
        {
            return MinSupportVector( dir, supp );
        }

        // Computes only the values of min/max support function. Usefull to compute bounding boxes.
        virtual void MinMaxSupportValue( cptr<Real> dir, mref<Real> min_val, mref<Real> max_val ) const override
        {
            min_val = dot_buffer<AMB_DIM>( &serialized_data[1], dir );
            max_val = min_val;
        }
        
        
        // Helper function to compute axis-aligned bounding boxes. in the format of box_min, box_max vector.
        // box_min, box_max are supposed to be vectors of size AMB_DIM.
        // BoxMinMax computes the "lower left" lo and "upper right" hi vectors of the primitives bounding box and sets box_min = min(lo, box_min) and box_max = min(h, box_max)
        virtual void BoxMinMax( mptr<SReal> box_min_, mptr<SReal> box_max_ ) const
        {
            copy_buffer<AMB_DIM>( &serialized_data[1], box_min_ );
            copy_buffer<AMB_DIM>( &serialized_data[1], box_max_ );
        }
        
        
        virtual std::string ClassName() const override
        {
            return std::string("Point")+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+","+TypeName<ExtReal>+","+TypeName<ExtInt>+">";
        }

    }; // Point
    
} // namespace Repulsor
