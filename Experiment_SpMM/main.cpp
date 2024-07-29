
#include <iostream>

#define TOOLS_ENABLE_PROFILER
#define TOOLS_AGGRESSIVE_INLINING

#ifdef __APPLE__
/// Use these while on a mac. Don't forget to issue the compiler flag `-framework Accelerate`.
///
    #include "../submodules/Tensors/Accelerate.hpp"
#else
/// This should work for OpenBLAS.
    #include "../submodules/Tensors/OpenBLAS.hpp"
#endif

#include "../submodules/Tensors/Tensors.hpp"
//#include "Repulsor.hpp"

//#include <locale>

//using namespace Repulsor;
using namespace Tensors;
using namespace Tools;

using Real    = Real64;
using Complex = Complex64;
using Int     = Int32;

using Scal    = Real;

int main(int argc, const char * argv[])
{

    const Int thread_count = 8;
    
    Sparse::MatrixCSR<Real,Int,Size_T> A;
    
    
    tic("LoadFromMatrixMarket");
//    A.LoadFromMatrixMarket( HomeDirectory() / "Triceratops_Laplacian.mtx", thread_count );
    A.LoadFromMatrixMarket( HomeDirectory() / "Triceratops_BiLaplacian.mtx", thread_count );
//    A.LoadFromMatrixMarket( HomeDirectory() / "MatrixMarket" / "bcsstk17.mtx", thread_count );
    toc("LoadFromMatrixMarket");
    
//    tic("WriteToMatrixMarket");
//    A.WriteToMatrixMarket( HomeDirectory() / "b.mtx" );
//    toc("WriteToMatrixMarket");
    
    
    dump( A.RowCount() );
    dump( A.ColCount() );
    
    
    print("");

    
    constexpr Int NRHS = 1;
    
    const Int ldX = 32;
    const Int ldY = 32;
    
    const Int p = 16;
    
    
    
    Tensor2<Scal,Int> X_wide ( A.RowCount(), ldX );
    Tensor2<Scal,Int> Y_wide ( A.RowCount(), ldY );
    
    X_wide.Random();
    Y_wide.Random();
    
    Tensor2<Scal,Int> X ( A.RowCount(), NRHS );
    Tensor2<Scal,Int> Y ( A.RowCount(), NRHS );
    Tensor2<Scal,Int> Z ( A.RowCount(), NRHS );
    
    X.Read( &X_wide[0][p], ldX );
    Y.Read( &Y_wide[0][p], ldY );
    
    Scal alpha = 1.9;
    Scal beta  = 1.2;
    
    Profiler::Clear();
    
    tic("Y = alpha * A.X + beta * Y");
    A.Dot<NRHS>( alpha, X.data(), NRHS, beta, Y.data(), NRHS, NRHS );
    toc("Y = alpha * A.X + beta * Y");
    
    
    tic("Y_wide = alpha * A.X_wide + beta * Y_wide");
    A.Dot<NRHS>( alpha, &X_wide[0][p], ldX, beta, &Y_wide[0][p], ldY, NRHS );
    toc("Y_wide = alpha * A.X_wide + beta * Y_wide");

    tic("Y = 1 * A.X + 0 * Y");
    A.Dot<NRHS>( 
        Scalar::One<Real>,  X.data(), NRHS,
        Scalar::Zero<Real>, Y.data(), NRHS,
        NRHS
    );
    toc("Y = alpha * A.X + beta * Y");
    
    tic("Y_wide = 1 * A.X_wide + 0 * Y_wide");
    A.Dot<NRHS>( 
        Scalar::One<Real>,  &X_wide[0][p], ldX,
        Scalar::Zero<Real>, &Y_wide[0][p], ldY,
        NRHS
    );
    toc("Y_wide = alpha * A.X_wide + beta * Y_wide");
    
    
    
    Z.Read(Y.data());
    
    combine_matrices<Scalar::Flag::Minus,Scalar::Flag::Plus,VarSize,NRHS,Par>
    ( 
     -1., &Y_wide[0][p], ldY,
      1., Z.data(),       NRHS,
        A.RowCount(), NRHS, thread_count
    );
    
    print("");
    dump( Z.MaxNorm() );
    print("");

    ptic("With load+write");
    X.Read( &X_wide[0][p], ldX, thread_count );
    Y.Read( &Y_wide[0][p], ldY, thread_count );
    A.Dot(
        Scalar::One <Real>, X.data(), NRHS,
        Scalar::Zero<Real>, Y.data(), NRHS,
        NRHS
    );
    X.Write( &X_wide[0][p], ldX, thread_count );
    Y.Write( &Y_wide[0][p], ldY, thread_count );
    ptoc("With load+write");

}
