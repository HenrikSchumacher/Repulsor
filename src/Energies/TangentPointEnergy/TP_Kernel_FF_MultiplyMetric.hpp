#pragma once

#define CLASS TP_Kernel_FF_MultiplyMetric
#define BASE  BlockKernel_gen<AMB_DIM_,AMB_DIM_,MAX_RHS_COUNT_,Scalar_,Int_,Scalar_in_,Scalar_out_,   \
    x_RM,  y_RM,                                                                                \
    alpha_flag, beta_flag                                                                       \
>

namespace Repulsor
{
    template<
        int AMB_DIM_, int MAX_RHS_COUNT_,
        typename Scalar_, typename Int_, typename Scalar_in_, typename Scalar_out_,
        bool x_RM, bool y_RM,
        int alpha_flag, int beta_flag
    >
    class CLASS : public BASE
    {
    public:

        using Scalar     = Scalar_;
        using Int        = Int_;
        using Scalar_out = Scalar_out_;
        using Scalar_in  = Scalar_in_;

        static constexpr Int AMB_DIM = AMB_DIM_;
        using BASE::ROWS;
        using BASE::COLS;
        using BASE::MAX_RHS_COUNT;
        
        static constexpr Int NONZERO_COUNT = 2 + AMB_DIM;
        
    protected:
        
        using BASE::A;
        using BASE::A_const;
        using BASE::X;
        using BASE::Y;
        using BASE::x;
        using BASE::z;
        using BASE::rhs_count;
        
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
        {
//            #pragma omp single
//            {
//                print(ClassName());
//                valprint("alpha     ",alpha);
//                valprint("alpha_flag",alpha_flag);
//                valprint("beta      ",beta);
//                valprint("beta_flag ",beta_flag);
//            }
        }
        
        // Copy constructor
        CLASS( const CLASS & other ) : BASE(other) {}
        
        virtual ~CLASS() override = default;
        
    public:
        
        virtual Int NonzeroCount() const override
        {
            return NONZERO_COUNT;
        }
                
        virtual force_inline void TransposeBlock( const Int from, const Int to ) const override
        {
            const Scalar * restrict const a_from = &A[ NONZERO_COUNT * from];
                  Scalar * restrict const a_to   = &A[ NONZERO_COUNT * to  ];
            
            a_to[0] = a_from[1];
            a_to[1] = a_from[0];
            for( int k = 2; k < NONZERO_COUNT; ++ k )
            {
                a_to[k] = -a_from[k];
            }
        }
        
        virtual force_inline void ApplyBlock( const Int block_id, const Int j_global ) override __attribute__ ((hot))
        {
            // Since we need the casted vector ROWS times, it might be a good idea to do the conversion only once.
            this->ReadVector( j_global );
            // It's a bit mysterious to me why copying to a local array makes this run a couple of percents faster.
            // Probably the copy has to be done anyways and this way the compiler has better guarantees.
            
            const Scalar * restrict const a = &A_const[NONZERO_COUNT * block_id];
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
//              |                                                                 |
//              \                                                                 /
//
//            This are 1 + 2 * AMB_DIM nonzero values.
//            With sparse matrix multiplication the bottleneck will be the bandwidth.
//            Hence instead, we store this data in only 2 + AMB_DIM values by just storing
//
//                   K_xy, K_yx, v[0], ..., v[AMB_DIM-1]!
            
            for( Int k = 0; k < rhs_count; ++k )
            {
                z[k][0] -= (a[0] + a[1]) * x[k][0];
                
                for( Int i = 1; i < AMB_DIM; ++i )
                {
                    z[k][0] += a[1] * a[i+1] * x[k][i];
                    z[k][i] -= a[0] * a[i+1] * x[k][0];
                }
            }
        }
        
    public:
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"
                +ToString(AMB_DIM)
            +","+ToString(MAX_RHS_COUNT)
            +","+TypeName<Scalar>::Get()
            +","+TypeName<Int>::Get()
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
