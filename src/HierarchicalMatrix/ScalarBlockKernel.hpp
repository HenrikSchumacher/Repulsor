#pragma once

#define CLASS ScalarBlockKernel
#define CLASS SquareBlockKernel<SIZE,T_,Int_,T_in_,T_out_>

namespace Repulsor
{
    template<int SIZE, typename T_, typename Int_, typename T_in_, typename T_out_ >
    class CLASS
    {
    public:

        using T    = T_;
        using Int  = Int_;
        using T_in = T_in_;
        using T_in = T_in_;

        constexpr Int NONZERO_COUNT = 1;
        
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
            A[to] = A[from];
        }
        
        virtual void ApplyBlock( const Int k, const Int j ) override
        {
            const T a_k = A[k];
            
            for( Int l = 0; l < SIZE; ++l )
            {
                z[l] += factor * static_cast<T>(X[ SIZE * j + l ]);
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


