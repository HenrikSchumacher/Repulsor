#pragma once

#define BASE BlockKernel_fixed<                                             \
        AMB_DIM_+1,AMB_DIM_+1,MAX_RHS_COUNT_, true,                         \
        Scalar_, Scalar_in_, Scalar_out_, Int_, LInt_,                      \
        alpha_flag, beta_flag,                                              \
        true, true, false, true,                                            \
        true, true,                                                         \
        false                                                               \
    >

namespace Repulsor
{
    template<
        int AMB_DIM_, int MAX_RHS_COUNT_,
        typename Scalar_, typename Scalar_in_, typename Scalar_out_, typename Int_, typename LInt_,
        int alpha_flag, int beta_flag
    >
    class TP0_Kernel_MultiplyMetric : public BASE
    
    {
    private:
        
        using Base_T = BASE;

    public:

        using Scalar     = Scalar_;
        using Scalar_out = Scalar_out_;
        using Scalar_in  = Scalar_in_;
        using Int        = Int_;
        using LInt       = LInt_;
        
        static constexpr Int AMB_DIM = AMB_DIM_;
        
        using Base_T::ROWS;
        using Base_T::COLS;
        using Base_T::MAX_RHS_COUNT;
        
        static constexpr LInt BLOCK_NNZ = 2;
        static constexpr LInt DIAG_NNZ  = 2;
        
    protected:
        
        using Base_T::A;
        using Base_T::A_const;
        using Base_T::X;
        using Base_T::Y;
        using Base_T::x;
        using Base_T::y;
        
        using Base_T::ReadX;
        using Base_T::get_x;
        using Base_T::get_y;
        using Base_T::FMA;
        
    public:
        
        TP0_Kernel_MultiplyMetric() = delete;
        
        explicit TP0_Kernel_MultiplyMetric( Scalar * restrict const A_ )
        :   Base_T( A_ )
        {}
        
        TP0_Kernel_MultiplyMetric(
            const Scalar     * restrict const A_,
            const Scalar_out                  alpha_,
            const Scalar_in  * restrict const X_,
            const Scalar_out                  beta_,
                  Scalar_out * restrict const Y_,
            const Int                         rhs_count_
        )
        :   Base_T( A_, alpha_, X_, beta_, Y_, rhs_count_ )
        {}
        
        // Copy constructor
        TP0_Kernel_MultiplyMetric( const TP0_Kernel_MultiplyMetric & other ) : Base_T(other) {}
        
        ~TP0_Kernel_MultiplyMetric() = default;
        
    public:
        
        static constexpr LInt NonzeroCount()
        {
            return BLOCK_NNZ;
        }
                

        force_inline void TransposeBlock( const LInt from, const LInt to ) const
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
            
/*            The metric block looks like this for AMB_DIM == 3:
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
*/
            
            const Scalar a_0 ( a[0] );
            
            LOOP_UNROLL_FULL
            for( Int k = 0; k < MAX_RHS_COUNT; ++k )
            {
                FMA( a_0, get_x(0,k), get_y(0,k) );
            }
            
            const Scalar a_1 ( a[1] );
            
            LOOP_UNROLL_FULL
            for( Int j = 1; j < COLS; ++j )
            {
                LOOP_UNROLL_FULL
                for( Int k = 0; k < MAX_RHS_COUNT; ++k )
                {
                    FMA( a_1, get_x(j,k), get_y(j,k) );
                }
            }
        }
        
    public:
        
        std::string ClassName() const
        {
            return "TP0_Kernel_MultiplyMetric<"
                +ToString(AMB_DIM)
            +","+ToString(MAX_RHS_COUNT)
            +","+TypeName<Scalar>::Get()+","+TypeName<Scalar_in>::Get()+","+TypeName<Scalar_out>::Get()
            +","+TypeName<Int>::Get()+","+TypeName<LInt>::Get()
            +","+ToString(alpha_flag)
            +","+ToString(beta_flag)
            +">";
        }

    }; // class TP0_Kernel_MultiplyMetric
    
} // namespace Repulsor


#undef BASE
