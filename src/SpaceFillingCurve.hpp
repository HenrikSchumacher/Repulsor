#pragma once


#include <bit>

namespace Repulsor
{

    // TODO: Make sure that only (signed?) integral types are used.
    // TODO: Make sure that the types and n and bit_count_ match.
    // TODO: Represent also the Hilbert code as Uint[n]?
    // TODO: Generate the reordings.
    // TODO: Input checks.
    // TODO: Guards against overflow.
    // TODO: Quality measures?
    // TODO: Input of floating point types? Translation and scaling.
    
    template<int n_, typename Real_, int bit_count_ = 16, typename Int_ = long long>
    class SpaceFillingCurve
    {
        /// n is the ambient dimensions of space filling curve.
        /// Pick bit_count_ so that 2^bit_count_ exceeds the largest value in any coordinate.
        
        static_assert(FloatQ<Real_>,"");
        
        static_assert(IntQ<Int_>,"");
        
        static_assert( bit_count_ <= 64, "Maximally 64 bits are supported."  );
        
        static_assert( bit_count_ > 0, "Positive number of bits is required."  );
        
    public:
        
        using Real = Real_;
        using Int  = Int_;
        
        using UInt = typename std::conditional<
            bit_count_ <= 8,
            uint8_t,
            typename std::conditional<
                bit_count_ <= 16,
                uint16_t,
                typename std::conditional<
                    bit_count_ <= 32,
                    uint32_t,
                    uint64_t
                >::type
            >::type
        >::type;
        
        static constexpr Int n = n_;
        
        static constexpr Int bit_count = std::numeric_limits<UInt>::digits;
        
        static constexpr Int total_bit_count  = bit_count * n;
        
        static constexpr Int byte_count = sizeof(UInt);
        
        static constexpr Int total_byte_count = sizeof(UInt) * n;
        

        using MortonCode_T = std::array<UInt,n>;
        
        using HilbertCode_T = std::array<UInt,n>;
        
        using Axes_T   = std::array<UInt,n>;
        
        using Vector_T = Tiny::Vector<n,Real,Int>;
        
        using BitField_T = std::bitset<total_bit_count>;
        
        using SmallBitField_T = std::bitset<n * 8>;
        
        
    protected:
        
        static constexpr UInt zero = static_cast<UInt>(0);
        static constexpr UInt one  = static_cast<UInt>(1);
        static constexpr UInt two  = static_cast<UInt>(2);
        
        
        const Int thread_count;
        
//        const UInt q_max;
        
        Vector_T lower;
        Vector_T upper;
        Vector_T mid;
        
        Tiny::Vector<256,BitField_T,Int> LUT;

        
        Real scale;
        
        Axes_T X_buffer = {};   // Axes
//        UInt Y_buffer [n] = {};   // TransposedIndex
//        UInt Z_buffer [n] = {};   // Index
//        UInt H_buffer [n] = {};   // HilbertIndex
        MortonCode_T M_buffer = {};   // MortonIndex
        
//        UInt * restrict const Y_scratch = &Y_buffer[0];
//        UInt * restrict const Z_scratch = &Z_buffer[0];
//        UInt * restrict const H_scratch = &M_buffer[0];
        
    private:
        
        void Init()
        {
//            TOOLS_DUMP(sizeof(BitField_T));
//            TOOLS_DUMP(sizeof(Axes_T));
//            TOOLS_DUMP(sizeof(MortonCode_T));

            TOOLS_PTIC("Creating lookup table");
            for( Int k = 0; k < 256; ++k )
            {
                BitField_T m = 0;
                
                for( Int j = 0; j < 8; ++j )
                {
//                    set_bit( m, n * j, get_bit( k, j ) );
                    
                    m[n * j] = get_bit( k, j );
                }

                LUT[k] = m;
            }
            TOOLS_PTOC("Creating lookup table");
            
//            TOOLS_DUMP(LUT);
        }
        
    public:
        
        
        SpaceFillingCurve( const Int thread_count_  )
        :   thread_count( thread_count_ )
//        ,   q_max( one << (bit_count - 1) )
        {
            Init();
        }
        
        // No default constructor
        SpaceFillingCurve()
        :   thread_count( Int(1) )
//        ,   q_max( one << (bit_count - 1) )
        {
            Init();
        }
        // Destructor
        ~SpaceFillingCurve() = default;
        // Copy constructor
        SpaceFillingCurve( const SpaceFillingCurve & other ) = default;
        // Copy assignment operator
        SpaceFillingCurve & operator=( const SpaceFillingCurve & other ) = default;
        // Move constructor
        SpaceFillingCurve( SpaceFillingCurve && other ) = default;
        // Move assignment operator
        SpaceFillingCurve & operator=( SpaceFillingCurve && other ) = default;
        
        // Conventions:
        // X - Coordinate vector with entries in [ 0, 2^bit_count -1 [.
        // Y - transposed Hilbert index as UInt [n]
        // H - Hilbert index as UInt [n].
        // index - Hilbert index as BigInt
        
        static constexpr Int BitCount()
        {
            return bit_count;
        }
        
        static constexpr Int AmbDim()
        {
            return n;
        }
        
        static constexpr Int TotalBitCount()
        {
            return total_bit_count;
        }
        
        Int ThreadCount() const
        {
            return thread_count;
        }
        
    public:
        
    
        
        void ComputeBoundingBox( cptr<Real> X, const Int point_count )
        {
            TOOLS_PTIC(ClassName()+"::ComputeBoundingBox");
            
            lower.Fill(Scalar::Max<Real>);
            upper.Fill(Scalar::Min<Real>);
            
            std::vector<Vector_T> thread_lower ( thread_count );
            std::vector<Vector_T> thread_upper ( thread_count );

            ParallelDo(
                [&]( const Int thread )
                {
                    Vector_T lower_loc = { Scalar::Max<Real> };
                    Vector_T upper_loc = { Scalar::Min<Real> };
                    
                    const Int i_begin = JobPointer( point_count, thread_count, thread     );
                    const Int i_end   = JobPointer( point_count, thread_count, thread + 1 );
                    
                    for( Int i = i_begin; i < i_end; ++ i )
                    {
                        Vector_T x ( &X[n*i] );
                    
                        for( Int k = 0; k < n; ++k )
                        {
                            lower_loc[k] = Min( lower_loc[k], x[k] );
                            upper_loc[k] = Max( upper_loc[k], x[k] );
                        }
                    }
                    
                    thread_lower[thread].Read( lower_loc.data() );
                    thread_upper[thread].Read( upper_loc.data() );
                },
                thread_count
            );
            
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                for( Int k = 0; k < n; ++k )
                {
                    lower[k] = Min( lower[k], thread_lower[thread][k] );
                    upper[k] = Max( upper[k], thread_upper[thread][k] );
                }
            }

            Real L = upper[0] - lower[0];
            
            mid[0] = Scalar::Half<Real> * ( lower[0] + upper[0] );
            
            for( Int k = 1; k < n; ++k )
            {
                mid[k] = Scalar::Half<Real> * ( lower[k] + upper[k] );
                
                L = Max( L, upper[k] - lower[k] );
            }
            
            scale = (Scalar::One<Real> - static_cast<Real>(128) * Scalar::eps<Real> ) / L;
            
            TOOLS_PTOC(ClassName()+"::ComputeBoundingBox");
        }
        
#include "SpaceFillingCurve/CoordsToAxes.hpp"
#include "SpaceFillingCurve/AxesToMorton.hpp"
#include "SpaceFillingCurve/MortonToAxes.hpp"
#include "SpaceFillingCurve/CoordsToMorton.hpp"
#include "SpaceFillingCurve/MortonOrdering.hpp"
//#include "SpaceFillingCurve/MortonOrdering_Radix.hpp"
//#include "SpaceFillingCurve/MortonOrdering_QuickSort.hpp"


#include "SpaceFillingCurve/AxesToHilbert.hpp"
#include "SpaceFillingCurve/CoordsToHilbert.hpp"
#include "SpaceFillingCurve/HilbertOrdering.hpp"
        
    public:
        
        static std::string ClassName()
        {
            return std::string("SpaceFillingCurve")+"<"+ToString(n)+","+ToString(bit_count)+","+TypeName<Int>+">";
        }
    };  // class SpaceFullingCurve
    
} // namespace Repulsor

