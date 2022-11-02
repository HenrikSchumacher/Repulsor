#pragma once

#define CLASS TP0_Kernel_MultiplyMetric
#define BASE  BlockKernel_fixed<                            \
    AMB_DIM_+1,AMB_DIM_+1,MAX_RHS_COUNT_,true,              \
    Scalar_, Scalar_in_, Scalar_out_, Int_, LInt_,          \
    alpha_flag, beta_flag,                                  \
    x_RM, false, false, false,                              \
    y_RM, false,                                            \
    false                                                   \
>

namespace Repulsor
{
    template<
        int AMB_DIM_, int MAX_RHS_COUNT_,
        typename Scalar_, typename Scalar_in_, typename Scalar_out_, typename Int_, typename LInt_,
        bool x_RM, bool y_RM,
        int alpha_flag, int beta_flag
    >
    class CLASS : public BASE
    {

    public:

        using Scalar     = Scalar_;
        using Scalar_out = Scalar_out_;
        using Scalar_in  = Scalar_in_;
        using Int        = Int_;
        using LInt       = LInt_;
        
        static constexpr Int AMB_DIM = AMB_DIM_;
        
        using BASE::ROWS;
        using BASE::COLS;
        using BASE::MAX_RHS_COUNT;
        
        static constexpr LInt BLOCK_NNZ = 2;
        static constexpr LInt DIAG_NNZ  = 2;
        
    protected:
        
        using BASE::A;
        using BASE::A_const;
        using BASE::X;
        using BASE::Y;
        using BASE::x;
        using BASE::y;
        
        using BASE::ReadX;
        
    public:
        
        CLASS() = delete;
        
        explicit CLASS( Scalar * restrict const A_ )
        :   BASE( A_ )
        {}
        
        CLASS(
            const Scalar     * restrict const A_,
            const Scalar_out                  alpha_,
            const Scalar_in  * restrict const X_,
            const Scalar_out                  beta_,
                  Scalar_out * restrict const Y_,
            const Int                         rhs_count_
        )
        :   BASE( A_, alpha_, X_, beta_, Y_, rhs_count_ )
        {}
        
        // Copy constructor
        CLASS( const CLASS & other ) : BASE(other) {}
        
        ~CLASS() = default;
        
    public:
        
        LInt NonzeroCount() const
        {
            return BLOCK_NNZ;
        }
                
        void TransposeBlock( const LInt from, const LInt to ) const
        {
            const Scalar * restrict const a_from = &A[ BLOCK_NNZ * from];
                  Scalar * restrict const a_to   = &A[ BLOCK_NNZ * to  ];
            
            a_to[0] = a_from[0];
            a_to[1] = a_from[1];
        }
        
        
        
        force_inline void ApplyBlock( const LInt k_global, const Int j_global )
        {
            ReadX( j_global );
            
            const Scalar * restrict const a = &A_const[BLOCK_NNZ * k_global];
//            The metric block looks like this for AMB_DIM == 3:
//
//              /                               \
//              |  a[0]     0       0       0   |
//              |                               |
//              |   0      a[1]     0       0   |
//              |                               |
//              |   0       0      a[1]     0   |
//              |                               |
//              |   0       0       0      a[1] |
//              \                               /
            
            #pragma unroll
            for( Int k = 0; k < MAX_RHS_COUNT; ++k )
            {
                y[k][0] += a[0] * x[k][0];
                
                #pragma unroll
                for( Int j = 1; j < COLS; ++j )
                {
                    y[k][j] += a[1] * x[k][j];
                }
            }
        }
        
    public:
        
        std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"
                +ToString(AMB_DIM)
            +","+ToString(MAX_RHS_COUNT)
            +","+TypeName<Scalar>::Get()+","+TypeName<Scalar_in>::Get()+","+TypeName<Scalar_out>::Get()
            +","+TypeName<Int>::Get()+","+TypeName<LInt>::Get()
            +","+ToString(x_RM)
            +","+ToString(y_RM)
            +","+ToString(alpha_flag)
            +","+ToString(beta_flag)
            +">";
        }

    };
} // namespace Repulsor

#undef BASE
#undef CLASS


