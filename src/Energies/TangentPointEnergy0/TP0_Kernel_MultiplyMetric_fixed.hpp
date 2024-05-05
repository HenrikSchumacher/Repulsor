#pragma once

#define BASE BlockKernel_fixed<                                       \
        AMB_DIM_+1,AMB_DIM_+1,MAX_RHS_COUNT_, true,                   \
        Real_, Real_in_, Real_out_, Int_, LInt_,                      \
        alpha_flag, beta_flag,                                        \
        true, true, false, true,                                      \
        true, true,                                                   \
        false                                                         \
    >

namespace Repulsor
{
    template<
        int AMB_DIM_, int MAX_RHS_COUNT_,
        typename Real_, typename Real_in_, typename Real_out_, typename Int_, typename LInt_,
        Scalar::Flag alpha_flag, Scalar::Flag beta_flag
    >
    class TP0_Kernel_MultiplyMetric_fixed : public BASE
    
    {
    private:
        
        using Base_T = BASE;

    public:

        using Real     = Real_;
        using Real_out = Real_out_;
        using Real_in  = Real_in_;
        using Int      = Int_;
        using LInt     = LInt_;
        
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
        using Base_T::x_from;
        using Base_T::y;
        
        using Base_T::ReadX;
        using Base_T::get_x;
        using Base_T::get_y;
        using Base_T::FMA;
        
    public:
        
        TP0_Kernel_MultiplyMetric_fixed() = delete;
        
        explicit TP0_Kernel_MultiplyMetric_fixed( mptr<Real> A_ )
        :   Base_T( A_ )
        {}
        
        TP0_Kernel_MultiplyMetric_fixed(
            cptr<Real>     A_,
            cref<Real_out> alpha_,  cptr<Real_in> X_,
            cref<Real_out> beta_,   mptr<Real>    Y_,
            const Int      rhs_count_
        )
        :   Base_T( A_, alpha_, X_, beta_, Y_, rhs_count_ )
        {}
        
        // Copy constructor
        TP0_Kernel_MultiplyMetric_fixed( cref<TP0_Kernel_MultiplyMetric_fixed> other ) : Base_T(other) {}
        
        ~TP0_Kernel_MultiplyMetric_fixed() = default;
        
    public:
        
        static constexpr LInt NonzeroCount()
        {
            return BLOCK_NNZ;
        }
                

        force_inline void TransposeBlock( const LInt from, const LInt to ) const
        {
            copy_buffer<2>( &A[ BLOCK_NNZ * from], &A[ BLOCK_NNZ * to  ] );
        }
        
        
        force_inline void ApplyBlock( const LInt k_global, const Int j_global )
        {
            ReadX( j_global );
            
            cptr<Real> a = &A_const[BLOCK_NNZ * k_global];
            
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

            
//            if constexpr ( vecQ )
//            {
//                std::array<vec_T<MAX_RHS_COUNT,Real>,COLS> x_vec;
//                std::array<vec_T<MAX_RHS_COUNT,Real>,ROWS> y_vec;
//
//                // This explicity copying is somewhat insane but will probably optimized a bit by the compiler.
//                for( Int j = 0; j < COLS; ++j )
//                {
//                    for( Int k = 0; k < MAX_RHS_COUNT; ++k )
//                    {
//                        x_vec[j][k] = get_x(j,k);
//                    }
//                }
//
//                // This explicity copying is somewhat insane but will probably optimized a bit by the compiler.
//                for( Int j = 0; j < ROWS; ++j )
//                {
//                    for( Int k = 0; k < MAX_RHS_COUNT; ++k )
//                    {
//                        y_vec[j][k] = get_y(j,k);
//                    }
//                }
//
//                const Real a_0 ( a[0] );
//
//                y_vec[0] += a_0 * x_vec[0];
//
//                const Real a_1 ( a[1] );
//
//                for( Int j = 1; j < COLS; ++j )
//                {
//                    y_vec[j] += a_1 * x_vec[j];
//                }
//
//                // This explicity copying is somewhat insane but will probably optimized a bit by the compiler.
//                for( Int j = 0; j < ROWS; ++j )
//                {
//                    for( Int k = 0; k < MAX_RHS_COUNT; ++k )
//                    {
//                        get_y(j,k) = y_vec[j][k];
//                    }
//                }
//            }
//            else
            {
                
                const Real a_0 ( a[0] );

                for( Int k = 0; k < MAX_RHS_COUNT; ++k )
                {
                    FMA( a_0, get_x(0,k), get_y(0,k) );
                }

                const Real a_1 ( a[1] );

                for( Int j = 1; j < COLS; ++j )
                {
                    for( Int k = 0; k < MAX_RHS_COUNT; ++k )
                    {
                        FMA( a_1, get_x(j,k), get_y(j,k) );
                    }
                }
            }
        }
        
    public:
        
        std::string ClassName() const
        {
            return "TP0_Kernel_MultiplyMetric_fixed<"
                +ToString(AMB_DIM)
            +","+ToString(MAX_RHS_COUNT)
            +","+TypeName<Real>+","+TypeName<Real_in>+","+TypeName<Real_out>
            +","+TypeName<Int>+","+TypeName<LInt>
            +","+ToString(alpha_flag)
            +","+ToString(beta_flag)
            +">";
        }

    }; // class TP0_Kernel_MultiplyMetric_fixed
    
} // namespace Repulsor


#undef BASE

