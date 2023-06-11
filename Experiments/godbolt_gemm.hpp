// see https://stackoverflow.com/questions/35620853/how-to-write-a-matrix-matrix-product-that-can-compete-with-eigen/35637007#35637007

template <typename T>
struct BlockSize
{};


template <>
struct BlockSize<double>
{
    static constexpr int MC     = BS_D_MC;
    static constexpr int KC     = BS_D_KC;
    static constexpr int NC     = BS_D_NC;
    static constexpr int MR     = BS_D_MR;
    static constexpr int NR     = BS_D_NR;

    static constexpr int rwidth = 128;
    static constexpr int align  = rwidth / 8;
    static constexpr int vlen   = rwidth / (8*sizeof(double));

    static_assert(MC>0 && KC>0 && NC>0 && MR>0 && NR>0, "Invalid block size.");
    static_assert(MC % MR == 0, "MC must be a multiple of MR.");
    static_assert(NC % NR == 0, "NC must be a multiple of NR.");
    static_assert(rwidth % sizeof(double) == 0, "SIMD register width not sane.");
};

#include <list>
#include <vector>
#include <memory>
#include <cstring>
#include <algorithm>

#define restrict __restrict__

#define force_inline inline

using Scal = double;
using Int    = int;

constexpr Int n = 4;

template<Int M, Int N, Int K>
inline void dot3(
    const Scal * restrict A,
    const Scal * restrict B,
          Scal * restrict C
)
{
    typedef Scal v_T __attribute__((vector_size (BlockSize<T>::rwidth/8)));

    static constexpr Int vlen = BlockSize<Scal>::vlen;
    static constexpr Int MR   = M;
    static constexpr Int NR   = N/vlen;

    v_T c[MR][NR] = {};
    A = (const Scal * ) __builtin_assume_aligned (A, BlockSize<double>::align);
    B = (const Scal * ) __builtin_assume_aligned (B, BlockSize<double>::align);

    for( size_t k = 0; k < K; ++k )
    {
        const v_T * b = (const v_T * ) B;
        for( size_t i = 0; i < MR; ++i )
        {
            for( size_t j = 0; j < NR; ++j )
            {
   
                c[i][j] += A[i] * b[j];
            }
        }
        A+=MR;
        B+=NR*vlen;
    }

    std::copy_n( (Scal *)&c[0][0], M*N, C );
}

template<Int M, Int N, Int K>
inline void dot2(
    const Scal * restrict const A_,
    const Scal * restrict const B_,
          Scal * restrict const C
)
{
    Scal c[M][N] = {};

    const Scal * restrict A = A_;
    const Scal * restrict B = B_;
    
    for( Int k = 0; k < K; ++k )
    {
        for( Int i = 0; i < M; ++i)
        {
            for( Int j = 0; j < N; ++j)
            {
                c[i][j] += A[i] * B[j];
            }
        }
        A+=M;
        B+=N;
    }

    std::copy_n( &c[0][0], M*N, C );
}


template<int M, int N, int K>
inline void dot(
    const Scal * restrict const a,
    const Scal * restrict const b,
          Scal * restrict const c
)
{
    Scal A[M][K];
    Scal B[N][K];
    Scal C[M][N] = {{}};

    std::copy_n( a, M*K, &A[0][0] );
    std::copy_n( b, K*N, &B[0][0] );

    for( size_t i = 0; i < M; ++i )
    {
        Scal C_i [M] = {};
        for( size_t j = 0; j < N; ++j )
        {
            for( size_t k = 0; k < K; ++k )
            {
                C[i][j] += A[i][k] * B[j][k];
            }
        }
    }

    std::copy_n( &C[0][0], M*N, c );
}

template void dot3<n,n,n>( const Scal * restrict const a, const Scal * restrict const b, Scal * restrict const c );
template void dot2<n,n,n>( const Scal * restrict const a, const Scal * restrict const b, Scal * restrict const c );
template void dot <n,n,n>( const Scal * restrict const a, const Scal * restrict const b, Scal * restrict const c );
