#pragma once

#define CLASS Parallelepiped
#define BASE  PrimitiveSerialized<AMB_DIM,Real,Int,SReal>

// TODO: Test MinMaxSupportValue

namespace Repulsor
{
    
    template<int AMB_DIM, typename Real, typename Int, typename SReal>
    class CLASS : public BASE
    {
    protected:
        
        // The Parallelepiped is the image of Cuboid[ {-1,...,-1}, {1,...,1} ] under the affine mapping x \mapsto transform * x+center.
        
        
        // serialized_data is assumed to be an array of size SIZE. Will never be allocated by class! Instead, it is meant to be mapped onto an array of type Real by calling the member SetPointer.
        
        // DATA LAYOUT
        // serialized_data[0] = squared radius
        // serialized_data[1],...,serialized_data[AMB_DIM] = center
        // serialized_data[AMB_DIM+1],...,serialized_data[AMB_DIM+POINT_COUNT x AMB_DIM] = transform.
        
    
    public:
        
        CLASS() : BASE() {}
        
        // Copy constructor
        CLASS( const CLASS & other ) : BASE( other ) {}

        // Move constructor
        CLASS( CLASS && other ) noexcept : BASE( other ) {}

        virtual ~CLASS() override = default;
        
        
        static constexpr Int SIZE = 1+AMB_DIM+AMB_DIM*AMB_DIM;
        
        virtual constexpr Int Size() const override
        {
            return SIZE;
        }
        
//    private:
//
//        mutable SReal self_buffer [SIZE];
        
    public:
        
#include "Primitive_Common.hpp"
        
        __ADD_CLONE_CODE__(CLASS)
        
        void FromTransform( cptr<SReal> center, cptr<SReal> transform ) const
        {
            mref<SReal> r2 = this->serialized_data[0];
            
            mptr<SReal> x  = this->serialized_data+1;
            mptr<SReal> A  = this->serialized_data+1+AMB_DIM;
            
            copy_buffer<AMB_DIM          >( center,    x );
            copy_buffer<AMB_DIM * AMB_DIM>( transform, A );

            
            // Computes the sum of the squared half-axis lengths which equals the sum of squared singular values which equals the square of the Frobenius norm.
            
            r2 = Scalar::Zero<SReal>;

            for( Int i = 0; i < AMB_DIM; ++i )
            {
                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    r2 += A[AMB_DIM*i+j]*A[AMB_DIM*i+j];
                }
            }
            // CAUTION: Up to this point, this only really the squared radius if the frame has orthogonal columns. Otherwise, r2 might be too small.
            
            
            // Adding a safety term that vanishes when frame has orthogonal column. This is to guarantee that r2 is indeed greater or equal to the actual squared radius.
            // This is superflous if the Parallelepiped is supposed to represent an OBB.
            SReal dots = Scalar::Zero<SReal>;

            for( Int i = 0; i < AMB_DIM; ++i )
            {
                for( Int j = i+1; j < AMB_DIM; ++j )
                {
                    SReal dot = A[i]*A[j];

                    for( Int k = 1; k < AMB_DIM; ++k )
                    {
                        dot += A[AMB_DIM*k+i]*A[AMB_DIM*k+j];
                    }
                    
                    dots += Abs( dot );
                }

            }
            r2 += Scalar::Two<SReal>*dots;
        }
        
        
        //Computes support vector supp of dir.
        virtual Real MaxSupportVector( cptr<Real> dir, mptr<Real> supp ) const override
        {
            cptr<SReal> x = this->serialized_data+1;
            cptr<SReal> A = this->serialized_data+1+AMB_DIM;

            Real R1;
            Real R2;
            Real R3 = Scalar::Zero<Real>;
                        
            copy_buffer<AMB_DIM>( x, supp );
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                R1 = static_cast<Real>(A[i]) * dir[0];
                
                for( Int j = 1; j < AMB_DIM; ++j )
                {
                    R1 += dir[j] * static_cast<Real>(A[AMB_DIM*j+i]);
                }
                
                R2 = (R1 >= Scalar::Zero<Real>) ? Scalar::One<Real> : -Scalar::One<Real>;
                
                R3 += static_cast<Real>(x[i]) * dir[i] + R1 * R2;
                
                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    supp[j] += static_cast<Real>(A[AMB_DIM*j+i])*R2;
                }
            }

            return R3;
        }
        
        
        //Computes support vector supp of dir.
        virtual Real MinSupportVector( cptr<Real> dir, mptr<Real> supp ) const override
        {
            const SReal * restrict const x = this->serialized_data+1;
            const SReal * restrict const A = this->serialized_data+1+AMB_DIM;

            Real R1;
            Real R2;
            Real R3 = Scalar::Zero<Real>;
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                supp[i] = static_cast<Real>(x[i]);
            }
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                R1 = static_cast<Real>(A[i]) * dir[0];
                
                for( Int j = 1; j < AMB_DIM; ++j )
                {
                    R1 += dir[j] * static_cast<Real>(A[AMB_DIM*j+i]);
                }
                
                R2 = (R1 >= Scalar::Zero<Real>) ? -Scalar::One<Real> : Scalar::One<Real>;
                
                R3 += static_cast<Real>(x[i]) * dir[i] + R1 * R2;
                
                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    supp[j] += static_cast<Real>(A[AMB_DIM*j+i])*R2;
                }
            }

            return R3;
        }
        
        // Computes only the values of min/max support function. Usefull to compute bounding boxes.
        virtual void MinMaxSupportValue( cptr<Real> dir, mref<Real> min_val, mref<Real> max_val ) const override
        {
            cptr<Real> x = this->serialized_data+1;
            cptr<Real> A = this->serialized_data+1+AMB_DIM;

            Real R1;
            Real R2;
            min_val = Scalar::Zero<Real>;
            max_val = Scalar::Zero<Real>;
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                R1 = static_cast<Real>(A[i]) * dir[0];
                
                for( Int j = 1; j < AMB_DIM; ++j )
                {
                    R1 += dir[j] * static_cast<Real>(A[AMB_DIM*j+i]);
                }
                
                R2 = (R1 >= Scalar::Zero<Real>) ? Scalar::One<Real> : -Scalar::One<Real>;
                Real x_i = static_cast<Real>(x[i]);
                min_val += x_i * dir[i] - R1 * R2;
                max_val += x_i * dir[i] + R1 * R2;
            }
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+">";
        }
        
    }; // CLASS
    
} // namespace Repulsor

#undef CLASS
#undef BASE
