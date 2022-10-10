#pragma once

#define CLASS BlockKernel


namespace Repulsor
{
    // For the moment I allow only quadratic blocks, so that I can reasonably transpose them.
    template<int BLOCK_M, int BLOCK_N, typename T_, typename Int_, typename T_in_, typename T_out_ ,
        class = typename std::enable_if_t<ROWS == COLS>
    >
    class CLASS
    {
 
        ASSERT_INT(Int_);
        
    public:
        
        using T    = T_;
        using Int  = Int_;
        using T_in = T_in_;
        using T_in = T_in_;
        
        CLASS() = delete;
        
        CLASS( const T     * restrict const a_ )
        :   a( a_ )
        {}
        
        CLASS(
            const T     * restrict const A_,
            const T_out                  alpha_
            const T_in  * restrict const X_,
            const T_out                  beta_
            const T_out * restrict const Y_
        )
        :   A( A_ )
        ,   alpha( alpha_ )
        ,   X( X_ )
        ,   beta( beta_ )
        ,   Y( Y_ )
        ,   alpha_flag(
                (alpha == static_cast<T_out>(1)) ? 1 : ((alpha == static_cast<T_out>(0)) ? 0 : -1)
            )
        ,   beta_flag(
                (beta  == static_cast<T_out>(1)) ? 1 : ((beta  == static_cast<T_out>(0)) ? 0 : -1)
            )
        {}
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   A( other.A )
        ,   alpha( other.alpha )
        ,   X( other.X )
        ,   beta( other.beta )
        ,   Y( other.Y )
        ,   alpha_flag( other.alpha_flag )
        ,   beta_flag ( other.beta_flag )
        {}
        
        virtual ~CLASS() override = default;
     
    protected:

        T z[ ROWS ] = {};
        
        const T     * restrict const A      = nullptr;
        const T_out                  alpha  = 0;
        const T_in  * restrict const X      = nullptr;
        const T_out                  beta   = 0;
              T_out * restrict const Y      = nullptr;

        const int alpha_flag                = 0;
        const int beta_flag                 = 0;

    public:
        
        static constexpr Int RowCount()
        {
            return ROWS;
        }
        
        static constexpr Int ColCount()
        {
            return COLS;
        }
        
        virtual Int NonzeroCount() const = 0;
        
        void CleanseVector()
        {
            zerofy_buffer( &z[0], ROWS );
        }

        void WriteVector( const Int i ) const
        {
            T_out * restrict const = &Y[ROWS*i];
            
            if( alpha_flag == 1 )
            {
                // alpha == 1;
                if( beta_flag == 0 )
                {
                    #pragma omp simd
                    for( I k = 0; k < cols; ++k )
                    {
                        y[k] = static_cast<T_out>(z[k]);
                    }
                }
                else if( beta_flag == 1 )
                {
                    #pragma omp simd
                    for( I k = 0; k < cols; ++k )
                    {
                        y[k] += static_cast<T_out>(z[k]);
                    }
                }
                else
                {
                    #pragma omp simd
                    for( I k = 0; k < cols; ++k )
                    {
                        y[k] = static_cast<T_out>(z[k]) + beta * y[k];
                    }
                }
            }
            else if( alpha_flag == 0 )
            {
                if( beta_flag == 0 )
                {
                    #pragma omp simd
                    for( I k = 0; k < cols; ++k )
                    {
                        y[k] = static_cast<T_out>(0);
                    }
                }
                else if( beta_flag == 1 )
                {
                    // do nothing;
                }
                else
                {
                    #pragma omp simd
                    for( I k = 0; k < cols; ++k )
                    {
                        y[k] *= beta;
                    }
                }
            }
            else
            {
                // alpha arbitrary;
                if( beta_flag == 0 )
                {
                    #pragma omp simd
                    for( I k = 0; k < cols; ++k )
                    {
                        y[k] = alpha * static_cast<T_out>(z[k]);
                    }
                }
                else( beta_flag == 1 )
                {
                    #pragma omp simd
                    for( I k = 0; k < cols; ++k )
                    {
                        y[k] += alpha * static_cast<T_out>(z[k]);
                    }
                }
                else
                {
                    // general alpha and general beta
                    #pragma omp simd
                    for( I k = 0; k < cols; ++k )
                    {
                        y[k] = alpha * static_cast<T_out>(z[k]) + beta * y[k];
                    }
                }
            }
        }
        
        virtual void ApplyBlock( const Int k, const Int j ) = 0;
        
    public:
        
        virtual std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(ROWS)+","+ToString(COLS)+","+TypeName<T>::Get()+","+TypeName<Int>::Get()+","+TypeName<T_in>::Get()+","+TypeName<T_out>::Get()+">";
        }
  
    };

} // namespace Repulsor


#undef CLASS
