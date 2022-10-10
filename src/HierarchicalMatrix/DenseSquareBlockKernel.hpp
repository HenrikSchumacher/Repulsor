#pragma once

#define CLASS DenseSquareBlockKernel
#define CLASS SquareBlockKernel<SIZE,T_,Int_,T_in_,T_out_>

namespace Repulsor
{
    template<int SIZE, typename T_, typename Int_, typename T_in_, typename T_out_>
    class CLASS
    {
    public:

        using T     = T_;
        using Int   = Int_;
        using T_out = T_out_;
        using T_in  = T_in_;

        constexpr Int NONZERO_COUNT = SIZE * SIZE;
        
        CLASS() = delete;
        
        CLASS(
            const T     * restrict const a_,
        )
        :   BASE(a_)
        {}
        
        CLASS(
            const T     * restrict const A_,
            const T_out                  alpha_
            const T_in  * restrict const X_,
            const T_out                  beta_
            const T_out * restrict const Y_
        )
        :   BASE( A_, alpha_, X_, beta_, Y_ )
        {}
        
        // Copy constructor
        CLASS( const CLASS & other ) : BASE(other) {}
        
        virtual ~CLASS() override = default;
     
    protected:

        using BASE::A;
        using BASE::X;
        using BASE::Y;
        using BASE::z;
        
    public:
        
        virtual Int NonzeroCount() const override
        {
            return NONZERO_COUNT;
        };
        
        virtual void TransposeBlock( const Int from, const Int to ) override
        {
            const Real * restrict const a_from = &A[ NONZERO_COUNT * from];
                  Real * restrict const a_to   = &A[ NONZERO_COUNT * to  ];
            
            for( Int j = 0; j < SIZE; ++j )
            {
                for( Int i = 0; i < SIZE; ++i )
                {
                    a_to[SIZE * j + i ] = a_from[SIZE * i + j ];
                }
            }
        }
        
        virtual void ApplyBlock( const Int block_id, const Int j_global ) override
        {
            const T * restrict const a  = &A[NONZERO_COUNT * block_id];
            // Since we need the casted vector ROWS times, it might be a good idea to do the conversion only one.
            T x [ SIZE ];
            
            copy_cast_buffer( &X[SIZE * j_global], &x[0], COLS );
      
            // TODO: SIMDization or offloading to a BLAS implementation.
            for( Int i = 0; i < SIZE; ++i )
            {
                for( Int j = 0; j < SIZE; ++j )
                {
                    z[i] += a[SIZE * i + j ] * x[j];
                }
            }
        }
        
    public:
        
        virtual std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(SIZE)+","+TypeName<T>::Get()+","+TypeName<Int>::Get()+","+TypeName<T_in>::Get()+","+TypeName<T_out>::Get()+">";
        }

    };

} // namespace Repulsor

#undef BASE
#undef CLASS


