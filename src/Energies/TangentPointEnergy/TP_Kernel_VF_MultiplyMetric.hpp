#pragma once

#define CLASS TP_Kernel_VF_MultiplyMetric
#define BASE  BlockKernel_fixed<                            \
    AMB_DIM_+1,AMB_DIM_+1,MAX_RHS_COUNT_,true,              \
    Scalar_,Scalar_in_,Scalar_out_,                         \
    Int_,LInt_,                                             \
    alpha_flag, beta_flag,                                  \
    x_RM, false, true, true,                                \
    y_RM, false,                                            \
    true                                                    \
>

namespace Repulsor
{
    template<
        int AMB_DIM_, int MAX_RHS_COUNT_,
        typename Scalar_, typename Scalar_in_, typename Scalar_out_,
        typename Int_, typename LInt_,
        bool x_RM, bool y_RM,
        int alpha_flag, int beta_flag
    >
    class CLASS : public BASE
    {
        
    public:

        using Scalar     = Scalar_;
        using Int        = Int_;
        using LInt       = LInt_;
        using Scalar_out = Scalar_out_;
        using Scalar_in  = Scalar_in_;

        static constexpr Int AMB_DIM = AMB_DIM_;
        
        using BASE::ROWS;
        using BASE::COLS;
        using BASE::MAX_RHS_COUNT;
        
        static constexpr Int BLOCK_NNZ = 1 + 2 * AMB_DIM;
        static constexpr Int DIAG_NNZ  = (1 + AMB_DIM) * (1 + AMB_DIM);
        
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
        
        virtual ~CLASS() override = default;
        
    public:
        
        virtual Int NonzeroCount() const override
        {
            return BLOCK_NNZ;
        }
                
        virtual void TransposeBlock( const Int from, const Int to ) const override
        {
            const Scalar * restrict const a_from = &A[ BLOCK_NNZ * from];
                  Scalar * restrict const a_to   = &A[ BLOCK_NNZ * to  ];
            
            a_to[0] = a_from[0];
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                a_to[1    + k] = a_from[ROWS + k];
                a_to[ROWS + k] = a_from[1    + k];
            }
        }
        
        
        
        virtual force_inline void apply_block( const Int k_global, const Int j_global ) override
        {
            // Since we need the casted vector ROWS times, it might be a good idea to do the conversion only once.
            ReadX( j_global );
            // It's a bit mysterious to me why copying to a local array makes this run a couple of percents faster.
            // Probably the copy has to be done anyways and this way the compiler has better guarantees.
            
            const Scalar * restrict const a = &A_const[BLOCK_NNZ * k_global];
//            The metric block looks like this for AMB_DIM == 3:
//
//              /                                                                 \
//              |   - K_xy - K_yx     Kyx * v[0]     Kyx * v[1]     Kyx * v[2]    |
//              |                                                                 |
//              |   - K_xy * v[0]         0              0              0         |
//              |                                                                 |
//              |   - K_xy * v[1]         0              0              0         |
//              |                                                                 |
//              |   - K_xy * v[2]         0              0              0         |
//              \                                                                 /
//
//            This are 1 + 2 * AMB_DIM nonzero values.
//            It is tempting to compress this to 2 + AMB_DIM values.
//            BUT we have to add the local matrices from several subtriangles!
//            Thus this structure cannot be exploited.
            
            for( Int k = 0; k < MAX_RHS_COUNT; ++k )
            {
                y[k][0] += a[0] * x[k][0];

                for( Int i = 1; i < ROWS; ++i )
                {
                    y[k][0] += a[i] * x[k][i];
                    y[k][i] += a[AMB_DIM+i] * x[k][0];
                }
            }
        }
        
        virtual force_inline void begin_row( const Int i_global ) override
        {}

        virtual force_inline void end_row( const Int i_global ) override
        {
            // TODO: Multiply diagonal block
        }
        
    public:
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"
                +ToString(AMB_DIM)
            +","+ToString(MAX_RHS_COUNT)
            +","+TypeName<Scalar>::Get()
            +","+TypeName<Int>::Get()
            +","+TypeName<LInt>::Get()
            +","+TypeName<Scalar_in>::Get()
            +","+TypeName<Scalar_out>::Get()
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

