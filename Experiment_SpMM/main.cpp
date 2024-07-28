
#include <iostream>

#define TOOLS_ENABLE_PROFILER

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

int main(int argc, const char * argv[]) 
{
    Profiler::Clear();
    
    Sparse::MatrixCSR<Real64,Int32,UInt64> A;
    
//    tic("LoadFromMatrixMarket");
//    A.LoadFromMatrixMarket( HomeDirectory() / "Triceratops_BiLaplacian.mtx", 4 );
//    toc("LoadFromMatrixMarket");
    
    tic("LoadFromMatrixMarket");
//    A.LoadFromMatrixMarket( HomeDirectory() / "Triceratops_BiLaplacian.mtx", 8 );
    A.LoadFromMatrixMarket( HomeDirectory() / "MatrixMarket" / "bcsstk17.mtx", 8);
    toc("LoadFromMatrixMarket");
    
    tic("WriteToMatrixMarket");
    A.WriteToMatrixMarket( HomeDirectory() / "b.mtx" );
    toc("WriteToMatrixMarket");

    
    dump( A.RowCount() );
    dump( A.ColCount() );

}
