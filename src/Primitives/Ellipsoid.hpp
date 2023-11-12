#pragma once

// TODO: Implement MinMaxSupportValue

#define BASE PrimitiveSerialized<AMB_DIM,Real,Int,SReal>
#define CLASS Ellipsoid

namespace Repulsor
{
    
    template<int AMB_DIM, typename Real, typename Int, typename SReal>
    class CLASS : public BASE
    {
    protected:
        
        // The Ellipsoid is the image of unit sphere under the affine mapping x \mapsto transform * x + center.
        
        // serialized_data is assumed to be an array of size SIZE. Will never be allocated by class! Instead, it is meant to be mapped onto an array of type Real by calling the member SetPointer.
        
        // DATA LAYOUT
        // serialized_data[0] = squared radius
        // serialized_data[1],...,serialized_data[AMB_DIM] = center
        // serialized_data[AMB_DIM + 1],...,serialized_data[AMB_DIM + POINT_COUNT x AMB_DIM] = transform.
        
    public:
        
        CLASS() : BASE() {}
        
        // Copy constructor
        CLASS( const CLASS & other ) : BASE( other ) {}

        // Move constructor
        CLASS( CLASS && other ) noexcept : BASE( other ) {}

        virtual ~CLASS() override = default;
        
        static constexpr Int SIZE = 1 + AMB_DIM + AMB_DIM * AMB_DIM;
        
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
        
        void FromTransform(  cptr<SReal> center,  cptr<SReal> transform ) const
        {
            mref<SReal> r2 = this->serialized_data[0];
            
            mptr<SReal> x = serialized_data + 1;
            mptr<SReal> A = serialized_data + 1 + AMB_DIM;
            
            copy_buffer<AMB_DIM          >( center,    x );
            copy_buffer<AMB_DIM * AMB_DIM>( transform, A );

            
            // Computes the maximum of the squared lengths of the columns of transform.
            // CAUTION: This is only really the squared radius if transform has orthogonal columns!!!
            // TODO: Find safer way to do this!

            for( Int i = 0; i < AMB_DIM; ++i )
            {
                this->SReal_buffer[i] = A[AMB_DIM * i] * A[AMB_DIM * i];
                for( Int k = 1; k < AMB_DIM; ++k )
                {
                    this->SReal_buffer[i] += A[ AMB_DIM * i + k] * A[ AMB_DIM * i + k];
                }
            }
            
            r2 = this->SReal_buffer[0];

            for( Int i = 1; i < AMB_DIM; ++i )
            {
                r2 = Max( r2, this->SReal_buffer[i]);
            }
            
        }
        
        //Computes support vector supp of dir.
        virtual Real MaxSupportVector( cptr<Real> dir, mptr<Real> supp ) const override
        {
            const SReal * restrict const x = this->serialized_data+ 1;
            const SReal * restrict const A = this->serialized_data+ 1 + AMB_DIM;
                   Real * restrict const b = this->Real_buffer;

            Real R1 = Scalar::Zero<Real>;
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                b[i] = static_cast<Real>(A[AMB_DIM * i]) * dir[0];

                for( Int j = 1; j < AMB_DIM; ++j )
                {
                    b[i] += dir[j] * static_cast<Real>(A[AMB_DIM * j + i]);
                }
                
                R1 += b[i] * b[i];
            }
            
            R1 = InvSqrt(R1);

            for( Int i = 0; i < AMB_DIM; ++i )
            {
                b[i] *= R1;
            }
            
            // Now this->Real_buffer is the max support vector on the unit sphere belonging to the director transform * dir.
            
            R1 = Scalar::Zero<Real>;
            
            // Transform the point back to the ellipsoid.
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                supp[i] = x[i] + static_cast<Real>(A[i]) * b[0];

                for( Int j = 1; j < AMB_DIM; ++j )
                {
                    supp[i] += static_cast<Real>(A[AMB_DIM * i + j]) * b[j];
                }
                
                R1 += supp[i] * dir[i];
            }

            return R1;
        }
        
        
        //Computes support vector supp of dir.
        virtual Real MinSupportVector( cptr<Real> dir, mptr<Real> supp ) const override
        {
            const SReal * restrict const x = this->serialized_data+ 1;
            const SReal * restrict const A = this->serialized_data+ 1 + AMB_DIM;
                   Real * restrict const b = this->Real_buffer;

            Real R1 = Scalar::Zero<Real>;
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                b[i] = A[AMB_DIM * i] * dir[0];

                for( Int j = 1; j < AMB_DIM; ++j )
                {
                    b[i] += dir[j] * static_cast<Real>(A[AMB_DIM * j + i]);
                }
                
                R1 += b[i] * b[i];
            }
            
            R1 = -InvSqrt(R1);

            for( Int i = 0; i < AMB_DIM; ++i )
            {
                b[i] *= R1;
            }
            
            // Now this->Real_buffer is the max support vector on the unit sphere belonging to the director transform * dir.
            
            R1 = Scalar::Zero<Real>;
            
            // Transform the point back to the ellipsoid.
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                supp[i] = static_cast<Real>(x[i]) + static_cast<Real>(A[i]) * b[0];

                for( Int j = 1; j < AMB_DIM; ++j )
                {
                    supp[i] += A[AMB_DIM * i + j] * b[j];
                }
                
                R1 += supp[i] * dir[i];
            }

            return R1;
        }
        
        virtual void MinMaxSupportValue( cptr<Real> dir, mref<Real> min_val, mref<Real> max_val ) const override
        {
            min_val = MinSupportVector( dir, &this->Real_buffer[0] );
            max_val = MaxSupportVector( dir, &this->Real_buffer[AMB_DIM] );
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+">";
        }

    }; // CLASS
    
} // namespace Repulsor

#undef CLASS
#undef BASE
