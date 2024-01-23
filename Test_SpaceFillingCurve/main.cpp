#include <iostream>


//#define TOOLS_DEBUG

#define TOOLS_ENABLE_PROFILER // enable profiler

#ifdef __APPLE__
/// Use these while on a mac. Don't forget to issue the compiler flag `-framework Accelerate`.
///
    #define ACCELERATE_NEW_LAPACK
    #include <Accelerate/Accelerate.h>
#else
/// Use these instead of Accelerate under Windows or Linux, e.g. together with OpenBLAS, Intel oneMKL, AMD AOCL-BLAS. Of course, your path variables or compiler flags should hint the compiler to these files. And you also have to link the according libraries.

/// This should work for OpenBLAS.

    #include <cblas.h>
    #include <lapack.h>

/// Use this with Intel oneMKL. Check their documentation if you are unsure.

    //#include <mkl_cblas.h>
    //#include <mkl_lapack.h>

#endif

#include "../Repulsor.hpp"
#include "../src/SpaceFillingCurve.hpp"

using namespace Repulsor;
using namespace Tensors;
using namespace Tools;


// We have to toggle which domain dimensions and ambient dimensions shall be supported by runtime polymorphism before we load Repulsor.hpp
// You can activate everything you want, but compile times might increase substatially.
using Int     = Int32;
using LInt    = Int64;
using ExtInt  = Int64;

using Real    = Real64;
using SReal   = Real64;
using ExtReal = Real64;


int main(int argc, const char * argv[])
{
    print("");
    print("###############################################################");
    print("###             Test program SpaceFillingCurve             ####");
    print("###############################################################");
    print("");
    
    Profiler::Clear( getenv("HOME") );

    
//    int thread_count = 8;
    int thread_count = 1;
    

    SpaceFillingCurve<3,Real> C ( thread_count );
    
    
    return 0;
}
