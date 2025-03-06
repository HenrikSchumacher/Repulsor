#pragma once

namespace Repulsor
{
    template<int AMB_DIM, typename Real, typename Int, typename SReal>
    class MovingPolytopeBase : public PrimitiveBase<AMB_DIM+1,Real,Int>
    {
    public:
        
        static_assert(FloatQ<SReal>,"");
        
        using Base_T = PrimitiveBase<AMB_DIM+1,Real,Int>;
        
    protected:
    
        mutable SReal a = Scalar::Zero<SReal>;
        mutable SReal b = Scalar::One<SReal>;
        
        mutable SReal T = Scalar::One<SReal>;
        
    public:
        
        MovingPolytopeBase() {}
        
        // Copy constructor
        MovingPolytopeBase( const MovingPolytopeBase & other )
        :   Base_T()
        ,   a(other.a)
        ,   b(other.b)
        ,   T(other.T)
        {}

        // Move constructor
        MovingPolytopeBase( MovingPolytopeBase && other ) noexcept 
        :   Base_T()
        ,   a(other.a)
        ,   b(other.b)
        ,   T(other.T)
        {}
        
        virtual ~MovingPolytopeBase() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<MovingPolytopeBase> Clone () const
        {
            return std::shared_ptr<MovingPolytopeBase>(CloneImplementation());
        }
        
    private:
        
        [[nodiscard]] virtual MovingPolytopeBase * CloneImplementation() const override = 0;
        
    public:
        
        virtual Int PointCount() const = 0;
        
        virtual Int Size() const = 0;
        
        virtual Int CoordinateSize() const = 0;
        
        virtual Int VelocitySize() const  = 0;
        
        virtual void SetFirstTime( const SReal a_ ) const
        {
            a = a_;
        }
        
        virtual void SetSecondTime( const SReal b_ ) const
        {
            b = b_;
        }
        
        virtual void SetTimeScale( const SReal T_ ) const
        {
            T = T_;
        }
        
        virtual SReal TimeScale() const
        {
            return T;
        }
        
        virtual void WriteCoordinatesSerialized(       SReal * const p_ ) const = 0;
        virtual void WriteCoordinatesSerialized(       SReal * const p_ser, const Int i ) const =0;
        
        virtual void ReadCoordinatesSerialized ( const SReal * const p_ ) = 0;
        virtual void ReadCoordinatesSerialized ( const SReal * const p_ser, const Int i ) = 0;

        virtual void WriteVelocitiesSerialized (       SReal * const v_ ) const = 0;
        virtual void WriteVelocitiesSerialized (       SReal * const v_ser, const Int i ) const = 0;
        
        virtual void ReadVelocitiesSerialized  ( const SReal * const v_ ) = 0;
        virtual void ReadVelocitiesSerialized  ( const SReal * const v_ser, const Int i ) = 0;
        
        virtual void WriteDeformedSerialized   (       SReal * const p_,    const SReal t ) const = 0;
        virtual void WriteDeformedSerialized   (       SReal * const p_ser, const SReal t, const Int i ) const = 0;
        
        virtual std::string ClassName() const override
        {
            return std::string("MovingPolytopeBase")+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+">";
        }
        
    }; // MovingPolytopeBase

    
} // namespace Repulsor
