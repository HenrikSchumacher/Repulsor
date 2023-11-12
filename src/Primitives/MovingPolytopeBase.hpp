#pragma once

#define CLASS MovingPolytopeBase
#define BASE  PrimitiveBase<AMB_DIM+1,Real,Int>

namespace Repulsor
{

    template<int AMB_DIM, typename Real, typename Int, typename SReal>
    class CLASS : public BASE
    {
        ASSERT_FLOAT(SReal);
        
    protected:
    
        mutable SReal a = Scalar::Zero<SReal>;
        mutable SReal b = Scalar::One<SReal>;
        
        mutable SReal T = Scalar::One<SReal>;
        
    public:
        
        CLASS() {}
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   BASE()
        ,   a(other.a)
        ,   b(other.b)
        ,   T(other.T)
        {}

        // Move constructor
        CLASS( CLASS && other ) noexcept 
        :   BASE()
        ,   a(other.a)
        ,   b(other.b)
        ,   T(other.T)
        {}
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE_FOR_ABSTRACT_CLASS__(CLASS)
        
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
        
        virtual SReal TimeScale( const SReal T_ ) const
        {
            return T;
        }
        
        virtual void WriteCoordinatesSerialized(       SReal * const p_ser, const Int i = 0 ) const =0;
        virtual void ReadCoordinatesSerialized ( const SReal * const p_ser, const Int i = 0 ) = 0;
        
        virtual void WriteVelocitiesSerialized (       SReal * const v_ser, const Int i = 0 ) const =0;
        virtual void ReadVelocitiesSerialized  ( const SReal * const v_ser, const Int i = 0 ) = 0;
        
        virtual void WriteDeformedSerialized   (       SReal * const p_ser, const SReal t, const Int i = 0 ) const = 0;
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+">";
        }
        
    };

    
} // namespace Repulsor

#undef CLASS
#undef BASE
