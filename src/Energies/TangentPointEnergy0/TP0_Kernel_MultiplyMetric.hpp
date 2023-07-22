#pragma once

#define BASE BlockKernel_fixed_2<                                     \
        AMB_DIM_+1,AMB_DIM_+1,NRHS,                                   \
        Scal_, Scal_in_, Scal_out_, Int_, LInt_,                      \
        alpha_flag, beta_flag,                                        \
        true, true,                                                   \
        true                                                          \
    >

namespace Repulsor
{
    template<
        int AMB_DIM_, int NRHS,
        typename Scal_, typename Scal_in_, typename Scal_out_, typename Int_, typename LInt_,
        Scalar::Flag alpha_flag, Scalar::Flag beta_flag
    >
    class TP0_Kernel_MultiplyMetric : public BASE
    
    {
    private:
        
        using Base_T = BASE;

    public:

        using Scal     = Scal_;
        using Scal_out = Scal_out_;
        using Scal_in  = Scal_in_;
        using Int      = Int_;
        using LInt     = LInt_;
        
        static constexpr Int AMB_DIM = AMB_DIM_;
        
        using Base_T::ROWS;
        using Base_T::COLS;
        using Base_T::MAX_RHS_COUNT;
        
        using BASE::vecQ;
        
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
        
    public:
        
        TP0_Kernel_MultiplyMetric() = delete;
        
        explicit TP0_Kernel_MultiplyMetric( mptr<Scal> A_ )
        :   Base_T( A_ )
        {}
        
        TP0_Kernel_MultiplyMetric(
            cptr<Scal>     A_,
            cref<Scal_out> alpha_,  cptr<Scal_in> X_,
            cref<Scal_out> beta_,   mptr<Scal>    Y_,
            const Int      rhs_count_
        )
        :   Base_T( A_, alpha_, X_, beta_, Y_, rhs_count_ )
        {}
        
        // Copy constructor
        TP0_Kernel_MultiplyMetric( cref<TP0_Kernel_MultiplyMetric> other ) : Base_T(other) {}
        
        ~TP0_Kernel_MultiplyMetric() = default;
        
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
            
            cptr<Scal> a = &A_const[BLOCK_NNZ * k_global];
            
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

            
            if constexpr ( vec_enabledQ )
            {
                
                const Scal a_0 ( a[0] );

                y[0] += a_0 * x[0];
                
                const Scal a_1 ( a[1] );
                
                for( Int j = 1; j < COLS; ++j )
                {
                    y[j] += a_1 * x[j];
                }
            }
            else
            {
                
                const Scal a_0 ( a[0] );

                for( Int k = 0; k < MAX_RHS_COUNT; ++k )
                {
                    y[0][k] += a_0 * x[0][k];
                }

                const Scal a_1 ( a[1] );

                for( Int j = 1; j < COLS; ++j )
                {
                    for( Int k = 0; k < MAX_RHS_COUNT; ++k )
                    {
                        y[j][k] += a_1 * x[j][k];
                    }
                }
            }
        }
        
    public:
        
        std::string ClassName() const
        {
            return "TP0_Kernel_MultiplyMetric<"
                +ToString(AMB_DIM)
            +","+ToString(MAX_RHS_COUNT)
            +","+TypeName<Scal>+","+TypeName<Scal_in>+","+TypeName<Scal_out>
            +","+TypeName<Int>+","+TypeName<LInt>
            +","+ToString(alpha_flag)
            +","+ToString(beta_flag)
            +">";
        }

    }; // class TP0_Kernel_MultiplyMetric
    
} // namespace Repulsor


#undef BASE
