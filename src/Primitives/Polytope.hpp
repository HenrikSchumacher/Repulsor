#pragma once

namespace Repulsor
{
    template<int POINT_COUNT, int AMB_DIM_, typename Real ,typename Int, typename SReal,
                typename ExtReal = SReal,typename ExtInt = Int>
    class Polytope : public PolytopeExt<AMB_DIM_,Real,Int,SReal,ExtReal,Int>
    {
        static_assert(FloatQ<ExtReal>,"");
        
        static_assert(IntQ<ExtInt>,"");
        
    protected:
        
        // this->serialized_data is assumed to be an array of size SIZE. Will never be allocated by class! Instead, it is meant to be mapped onto an array of type Real by calling the member SetPointer.
        
        // DATA LAYOUT
        // serialized_data[0] = squared radius
        // serialized_data[1],...,serialized_data[AMB_DIM] = interior_point
        // serialized_data[AMB_DIM + 1],...,serialized_data[AMB_DIM + POINT_COUNT x AMB_DIM] = points whose convex hull defines the polytope.
        
    public:
        
        using Base_T = PolytopeExt<AMB_DIM_,Real,Int,SReal,ExtReal,Int>;
        
        static constexpr Int AMB_DIM = AMB_DIM_;
        
        Polytope() : Base_T() {}
        
        // Copy constructor
        Polytope( const Polytope & other ) : Base_T( other ) {}

        // Move constructor
        Polytope( Polytope && other ) noexcept : Base_T( other ) {}

        virtual ~Polytope() override = default;
        
        static constexpr Int SIZE = 1 + AMB_DIM + POINT_COUNT * AMB_DIM;
        
        virtual Int Size() const override
        {
            return SIZE;
        }
        
        virtual Int PointCount() const override
        {
            return POINT_COUNT;
        }
        
#include "Primitive_Common.hpp"
        
    public:
        
        [[nodiscard]] std::shared_ptr<Polytope> Clone () const
        {
            return std::shared_ptr<Polytope>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual Polytope * CloneImplementation() const override
        {
            return new Polytope(*this);
        }
        
    public:
        
        virtual void FromCoordinates( cptr<ExtReal> hull_coords_, const Int i = 0 ) const override
        {
            constexpr SReal w = Inv<SReal>(POINT_COUNT);
            
            mref<SReal> r2 = this->serialized_data[0];
            
            cptr<ExtReal> A           = &hull_coords_[i * (POINT_COUNT * AMB_DIM)];
            mptr<SReal>   center      = &this->serialized_data[1];
            mptr<SReal>   hull_coords = &this->serialized_data[1 + AMB_DIM];
            
            // Copy the hull coordinates and compute average.
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                SReal x = static_cast<SReal>(A[ AMB_DIM * 0 + k ]);
                hull_coords[ AMB_DIM * 0 + k ] = x;

                center[k] = x;

                for( Int j = 1; j < POINT_COUNT; ++j )
                {
                    x = static_cast<SReal>(A[ AMB_DIM * j + k]);
                    hull_coords[ AMB_DIM * j + k ] = x;
                    center[k] += x;
                }
                center[k] *= w;
            }
            
            // Compute radius.
            r2 = Scalar::Zero<SReal>;

            for( Int j = 0; j < POINT_COUNT; ++j )
            {
                SReal diff = A[ AMB_DIM * j] - center[0];
                SReal square = diff * diff;

                for( Int k = 1; k < AMB_DIM; ++k )
                {
                    diff = A[ AMB_DIM * j + k] - center[k];
                    square += diff * diff;
                }
                r2 = Max( r2, square );
            }
        }
        
        virtual void FromIndexList( cptr<ExtReal> coords_, cptr<ExtInt> tuples, const Int i = 0 ) const override
        {
            cptr<ExtInt> s = tuples + POINT_COUNT * i;
            
            constexpr SReal w = Inv<SReal>(POINT_COUNT);
            
            mref<SReal> r2          = this->serialized_data[0];
            mptr<SReal> center      = &this->serialized_data[1];
            mptr<SReal> hull_coords = &this->serialized_data[1 + AMB_DIM];
            
            for( Int j = 0; j < POINT_COUNT; ++j )
            {
                copy_buffer<AMB_DIM>( &coords_[AMB_DIM * s[j]], &hull_coords[AMB_DIM * j] );
            }
            
            // Compute average.
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                center[k] = hull_coords[ AMB_DIM * 0 + k ];

                for( Int j = 1; j < POINT_COUNT; ++j )
                {
                    center[k] += hull_coords[ AMB_DIM * j + k ];
                }
                
                center[k] *= w;
            }
            
            // Compute radius.
            r2 = Scalar::Zero<Real>;

            for( Int j = 0; j < POINT_COUNT; ++j )
            {
                SReal diff = hull_coords[ AMB_DIM * j] - center[0];
                SReal square = diff * diff;
                
                for( Int k = 1; k < AMB_DIM; ++k )
                {
                    diff = hull_coords[ AMB_DIM * j + k] - center[k];
                    square += diff * diff;
                }
                r2 = Max( r2, square );
            }
        }
        
        
        
        //Computes support vector supp of dir.
        virtual Real MinSupportVector( cptr<Real> dir, mptr<Real> supp ) const override
        {
            // apply dot product with direction to all the points that span the convex hull
            cptr<SReal> A = &this->serialized_data[1 + AMB_DIM];

            Real value = dot_buffers<AMB_DIM>( A, dir );

            Int pos = 0;
            Real minimum = value;

            for( Int j = 1; j < POINT_COUNT; ++j )
            {
                value = dot_buffers<AMB_DIM>( &A[ AMB_DIM * j ], dir );

                if( value < minimum )
                {
                    pos = j;
                    minimum = value;
                }
            }

            copy_buffer<AMB_DIM>( &A[ AMB_DIM * pos ], supp );

            return minimum;
        }
        
        //Computes support vector supp of dir.
        virtual Real MaxSupportVector( cptr<Real> dir, mptr<Real> supp ) const override
        {
            // apply dot product with direction to all the points that span the convex hull
            cptr<SReal> A = &this->serialized_data[1 + AMB_DIM];

            Real value = dot_buffers<AMB_DIM>( A, dir );

            Int pos = 0;
            Real maximum = value;

            for( Int j = 1; j < POINT_COUNT; ++j )
            {
                value = dot_buffers<AMB_DIM>( &A[ AMB_DIM * j ], dir );

                if( value > maximum )
                {
                    pos = j;
                    maximum = value;
                }
            }

            copy_buffer<AMB_DIM>( &A[ AMB_DIM * pos ], supp );
            
            return maximum;
        }

        // Computes only the values of min/max support function. Usefull to compute bounding boxes.
        virtual void MinMaxSupportValue( cptr<Real> dir, mref<Real> min_val, mref<Real> max_val ) const override
        {
            // apply dot product with direction to all the points that span the convex hull
            cptr<SReal> A = &this->serialized_data[1 + AMB_DIM];

            Real value = dot_buffers<AMB_DIM>( A, dir );

            min_val = value;
            max_val = value;
            
            for( Int j = 1; j < POINT_COUNT; ++j )
            {
                value = dot_buffers<AMB_DIM>( &A[ AMB_DIM * j ], dir );

                min_val = Min( min_val, value );
                max_val = Max( max_val, value );
            }
        }
        
        
        // Helper function to compute axis-aligned bounding boxes. in the format of box_min, box_max vector.
        // box_min, box_max are supposed to be vectors of size AMB_DIM.
        // BoxMinMax computes the "lower left" lo and "upper right" hi vectors of the primitives bounding box and sets box_min = min(lo, box_min) and box_max = min(h, box_max)
        virtual void BoxMinMax( mptr<SReal> box_min, mptr<SReal> box_max ) const override
        {
            cptr<SReal> p = &this->serialized_data[1 + AMB_DIM];
            
            Tiny::Vector<AMB_DIM,SReal,Int> L;
            Tiny::Vector<AMB_DIM,SReal,Int> U;
            
            L.Read(p);
            U.Read(p);
            
            for( Int j = 1; j < POINT_COUNT; ++j )
            {
                for( Int k = 0; k < AMB_DIM; ++k )
                {
                    const SReal x = p[ AMB_DIM * j + k ];
                    L[k] = Min( L[k], x );
                    U[k] = Max( U[k], x );
                }
            }
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                box_min[k] = Min( L[k], box_min[k] );
                box_max[k] = Max( U[k], box_max[k] );
            }
        }
        
        virtual std::string ClassName() const override
        {
            return std::string("Polytope")+"<"+ToString(POINT_COUNT)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<SReal>+","+TypeName<ExtReal>+","+TypeName<ExtInt>+">";
        }

    }; // Polytope
    
    template <int AMB_DIM_, typename Real, typename Int, typename SReal,
        typename ExtReal = SReal, typename ExtInt = Int>
    [[nodiscard]] 
    std::shared_ptr<PolytopeExt<AMB_DIM_,Real,Int,SReal,ExtReal,Int>>
    MakePolytope( const Int P_size )
    {
        using Base_T = PolytopeExt<AMB_DIM_,Real,Int,SReal,ExtReal,Int>;
        Base_T * r;
        switch(  P_size  )
        {
            case 1:
            {
                r = new Polytope<1,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 2:
            {
                r = new Polytope<2,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 3:
            {
                r = new Polytope<3,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 4:
            {
                r = new Polytope<4,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 5:
            {
                r = new Polytope<5,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 6:
            {
                r = new Polytope<6,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 7:
            {
                r = new Polytope<7,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 8:
            {
                r = new Polytope<8,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 9:
            {
                r = new Polytope<9,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 10:
            {
                r = new Polytope<10,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 11:
            {
                r = new Polytope<11,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 12:
            {
                r = new Polytope<12,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 13:
            {
                r = new Polytope<13,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 14:
            {
                r = new Polytope<14,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 15:
            {
                r = new Polytope<15,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 16:
            {
                r = new Polytope<16,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 17:
            {
                r = new Polytope<17,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 18:
            {
                r = new Polytope<18,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 19:
            {
                r = new Polytope<19,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            case 20:
            {
                r = new Polytope<20,AMB_DIM_,Real,Int,SReal,ExtReal,ExtInt>();
                break;
            }
            default:
            {
                eprint("MakePolytope: Number of vertices of polytope = " + ToString( P_size ) + " is not in the range from 1 to 20. Returning nullptr.");
                return nullptr;
            }
        }
        
        return std::shared_ptr<Base_T>(r);
        
    } // MakePolytope
    
} // namespace Repulsor
