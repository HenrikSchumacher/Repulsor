#pragma once

namespace Repulsor
{
    // Real  -  data type that will be handed to GJK; GJK typically needs doubles.
    // Int   -  integer type for return values and loops.
    
    template<int AMB_DIM, typename Real_, typename Int_>
    class alignas(ObjectAlignment) PrimitiveBase // Use this broad alignment to prevent false sharing.
    {
        static_assert(FloatQ<Real_>,"");
        
        static_assert(IntQ<Int_>,"");

    public:
        
        using Real = Real_;
        using Int  = Int_;
        
        PrimitiveBase() = default;
        
        virtual ~PrimitiveBase() = default;
                                 
    public:
        
        [[nodiscard]] std::shared_ptr<PrimitiveBase> Clone () const
        {
            return std::shared_ptr<PrimitiveBase>(CloneImplementation());
        }
        
    private:
        
        [[nodiscard]] virtual PrimitiveBase * CloneImplementation() const = 0;
        
    public:
        
        static constexpr Int AmbDim()
        {
            return AMB_DIM;
        }
        
        // Computes support vector supp of dir.
        virtual Real MaxSupportVector( cptr<Real> dir, mptr<Real> supp ) const = 0;
        
        // Computes support vector supp of dir.
        virtual Real MinSupportVector( cptr<Real> dir, mptr<Real> supp ) const = 0;
        
        // Computes only the values of min/max support function. Usefull to compute bounding boxes.
        virtual void MinMaxSupportValue( cptr<Real> dir, mref<Real> min_val, mref<Real> max_val ) const = 0;
        
        // Returns some point within the primitive and writes it to p.
        virtual void InteriorPoint( mptr<Real> p ) const = 0;
    
        virtual Real InteriorPoint( const Int k ) const = 0;
        
        // Returns some (upper bound of the) squared radius of the primitive as measured from the result of InteriorPoint.
        virtual Real SquaredRadius() const = 0;

        virtual std::string DataString() const = 0;
        
        virtual std::string ClassName() const
        {
            return std::string("PrimitiveBase<")+TypeName<Real>+","+ToString(AMB_DIM)+">";
        }
        
    };
    
} // namespace Repulsor


