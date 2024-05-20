#pragma once

namespace Repulsor
{
    // Pure interface class

    template<int AMB_DIM, typename Real, typename Int,typename SReal, typename ExtReal, typename ExtInt>
    class PolytopeExt : public PolytopeBase<AMB_DIM,Real,Int,SReal>
    {
        static_assert(FloatQ<ExtReal>,"");
        
        static_assert(IntQ<ExtInt>,"");
        
    protected:
        
        // serialized_data is assumed to be an array of size SIZE. Will never be allocated by class! Instead, it is meant to be mapped onto an array of type Real by calling the member SetPointer.
        
        // DATA LAYOUT
        // serialized_data[0] = squared radius
        // serialized_data[1],...,serialized_data[AMB_DIM] = interior_point
        // serialized_data[1+ AMB_DIM],...,serialized_data[SIZE] = data that defines the primitive.
        
    public:
        
        using Base_T = PolytopeBase<AMB_DIM,Real,Int,SReal>;
        
        PolytopeExt() : Base_T () {}
        
        // Copy constructor
        PolytopeExt( const PolytopeExt & other ) : Base_T( other)
        {}

        // Move constructor
        PolytopeExt( PolytopeExt && other ) noexcept : Base_T( other)
        {}
        
        virtual ~PolytopeExt() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<PolytopeExt> Clone () const
        {
            return std::shared_ptr<PolytopeExt>(CloneImplementation());
        }
        
    private:
        
        [[nodiscard]] virtual PolytopeExt * CloneImplementation() const override = 0;
        
    public:
        
        virtual void FromCoordinates( cptr<ExtReal> hull_coords_, const Int i = 0 ) const = 0;
        
        virtual void FromIndexList( cptr<ExtReal> coords_, cptr<ExtInt> tuples, const Int i = 0 ) const = 0;
        
        virtual std::string ClassName() const override
        {
            return std::string("PolytopeExt")+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }
        
    };

} // namespace Repulsor
