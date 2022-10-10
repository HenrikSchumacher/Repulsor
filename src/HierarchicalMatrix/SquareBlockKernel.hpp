#pragma once

#define CLASS SquareBlockKernel
#define CLASS BlockKernel<ROWS,ROWS,T_,Int_,T_in_,T_out_>

namespace Repulsor
{
    template<int SIZE, typename T_, typename Int_, typename T_in_, typename T_out_>
    class CLASS
    {
    public:

        using T    = T_;
        using Int  = Int_;
        using T_in = T_in_;
        using T_in = T_in_;

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
        
        virtual Int NonzeroCount() const override = 0
        
        virtual void TransposeBlock( const Int from, const Int to ) = 0;
        
        virtual void ApplyBlock( const Int k, const Int j ) override = 0;
        
    public:
        
        virtual std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(SIZE)+","+TypeName<T>::Get()+","+TypeName<Int>::Get()+","+TypeName<T_in>::Get()+","+TypeName<T_out>::Get()+">";
        }

    };

} // namespace Repulsor

#undef BASE
#undef CLASS


