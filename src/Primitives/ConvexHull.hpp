#pragma once

//TODO: Test this thoroughly!

namespace Repulsor
{
    template<int HULL_COUNT, int AMB_DIM, typename Real, typename Int>
    class ConvexHull : public PrimitiveBase<AMB_DIM,Real,Int>
    {
    public:
        
        using Base_T = PrimitiveBase<AMB_DIM,Real,Int>;
        
    protected:
        
        const Base_T * primitive [HULL_COUNT];
        
        mutable Real buffer [AMB_DIM * (AMB_DIM+1)];
        mutable Real squared_radius = - Scalar::One<Real>;
        
    public:
        
        ConvexHull() : Base_T()
        {
            for( Int i = 0; i < HULL_COUNT; ++i )
            {
                primitive[i] = nullptr;
            }
        }
        
        virtual ~ConvexHull() override = default;

        void Set( const Int i, const Base_T * const  P )
        {
            if( 0 <= i < HULL_COUNT )
            {
                primitive[i] = P;
            }
        }
        
        bool NullQ( const Int i ) const
        {
            return ( primitive[i] == nullptr );
        }
        
        // Computes some point within the primitive and writes it to p.
        virtual void InteriorPoint( mptr<Real> point ) const override
        {
            Int count = Scalar::Zero<Int>;
            
            zerofy_buffer<AMB_DIM>(point);
            
            for( Int i = 0; i < HULL_COUNT; ++i )
            {
                if( primitive[i] )
                {
                    ++count;
                    primitive[i]->InteriorPoint( &this->buffer[0] );
                    
                    for( Int k = 0; k < AMB_DIM; ++k )
                    {
                        point[k] += this->buffer[k];
                    }
                }
            }
            
            count = Ramp_1(count);
            
            scale_buffer<AMB_DIM>( Inv<Real>(count), point );
        }
        
        
        //Computes support vector supp of dir.
        virtual Real MaxSupportVector( cptr<Real> dir, mptr<Real> supp ) const override
        {
            mptr<Real> b = &this->buffer[0];
            mptr<Real> b_max = &this->buffer[AMB_DIM];
            Real maximum = std::numeric_limits<Real>::lowest();
            Real value;
            
            for( Int i = 0; i < HULL_COUNT; ++i )
            {
                if( primitive[i] )
                {
                    value = primitive[i]->MaxSupportVector(dir, b);
                    
                    if( value > maximum )
                    {
                        maximum = value;
                        std::swap( b, b_max );
                    }
                }
            }
            
            copy_buffer<AMB_DIM>( b_max, supp );

            return maximum;
        }
        
        
        //Computes support vector supp of dir.
        virtual Real MinSupportVector( cptr<Real> dir, mptr<Real> supp ) const override
        {
            mptr<Real> b     = &this->buffer[0];
            mptr<Real> b_min = &this->buffer[AMB_DIM];
            Real minimum = std::numeric_limits<Real>::max();
            Real value;
            
            for( Int i = 0; i < HULL_COUNT; ++i )
            {
                if( primitive[i] )
                {
                    value = primitive[i]->MinSupportVector(dir, b);
                    
                    if( value < minimum )
                    {
                        minimum = value;
                        std::swap( b, b_min );
                    }
                }
            }
            
            copy_buffer<AMB_DIM>( b_min, supp );

            return minimum;
        }
        
        virtual void MinMaxSupportValue( cptr<Real> dir, mref<Real> min_val, mref<Real> max_val ) const override
        {
            min_val = MinSupportVector( dir, &this->buffer[0] );
            max_val = MaxSupportVector( dir, &this->buffer[AMB_DIM] );
        }
        
        virtual Real SquaredRadius() const override
        {
            // Computes the sum of the squared axis lengths which equals the sum of squared singular values which equals the square of the Frobenius norm.
            
            if( squared_radius >= Scalar::Zero<Real> )
            {
                return squared_radius;
            }
            
            squared_radius = Scalar::Zero<Real>;
            
            for( Int i = 0; i < HULL_COUNT; ++i )
            {
                if( primitive[i] )
                {
                    squared_radius = Max( squared_radius, primitive[i]->SquaredRadius() );
                }
            }
            
            return squared_radius;
        }
        
        
        // Copy constructor
        ConvexHull( const ConvexHull & B )
        {
            for( Int i = 0; i < HULL_COUNT; ++ i)
            {
                primitive[i] = B.primitive[i];
            }
            squared_radius = B.squared_radius;
        }
        
        // Move constructor
        ConvexHull( ConvexHull && B ) noexcept
        {
            for( Int i = 0; i < HULL_COUNT; ++ i )
            {
                primitive[i] = B.primitive[i];
                B.primitive[i] = nullptr;
            }
            squared_radius = std::move(B.squared_radius);
        }
        
        
    public:
        
        virtual std::string ClassName() const override
        {
            return std::string("ConvexHull")+"<"+ToString(HULL_COUNT)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+">";
        }
        
    }; // ConvexHull
    
} // namespace Repulsor
