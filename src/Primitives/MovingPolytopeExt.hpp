#pragma once

namespace Repulsor
{
    
    template<int AMB_DIM, typename Real,typename Int, typename SReal, typename ExtReal, typename ExtInt>
    class MovingPolytopeExt : public MovingPolytopeBase<AMB_DIM,Real,Int,SReal>
    {
    public:
        
        static_assert(FloatQ<ExtReal>,"");
        
        using Base_T = MovingPolytopeBase<AMB_DIM,Real,Int,SReal>;
        
    protected:

        using Base_T::a;
        using Base_T::b;
        using Base_T::T;
        
    public:
        
        MovingPolytopeExt() {}
        
        // Copy constructor
        MovingPolytopeExt( const MovingPolytopeExt & other )
        :   Base_T(other)
        {}

        // Move constructor
        MovingPolytopeExt( MovingPolytopeExt && other ) noexcept 
        :   Base_T(other)
        {}
        
        virtual ~MovingPolytopeExt() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<MovingPolytopeExt> Clone () const
        {
            return std::shared_ptr<MovingPolytopeExt>(CloneImplementation());
        }
        
    private:
        
        [[nodiscard]] virtual MovingPolytopeExt * CloneImplementation() const override = 0;
        
    public:
        
        virtual void FromCoordinates( const ExtReal * const p_, const Int i = 0 ) = 0;
        
        virtual void FromVelocities ( const ExtReal * const v_, const Int i = 0 ) = 0;
        
        virtual void FromVelocitiesIndexList(
            const ExtReal * const v_, const ExtInt * const tuples, const Int i ) = 0;
        
        virtual std::string ClassName() const override
        {
            return std::string("MovingPolytopeExt")+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }
        
    }; //

} // namespace Repulsor
