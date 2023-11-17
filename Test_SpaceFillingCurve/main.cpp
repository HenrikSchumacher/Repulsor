#include <iostream>

#include <sys/types.h>
#include <pwd.h>



//#define REMESHER_VERBATIM

#define NDEBUG
#define TOOLS_ENABLE_PROFILER // enable profiler
#define TOOLS_DEBUG

#define LAPACK_DISABLE_NAN_CHECK

#ifdef __APPLE__
    #define ACCELERATE_NEW_LAPACK
    #include <Accelerate/Accelerate.h>
#else
    #include <cblas.h>
    #include <lapack.h>

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
    
    
    const char * homedir = getenv("HOME");

    if( homedir == nullptr)
    {
        homedir = getpwuid(getuid())->pw_dir;
    }
    std::string path ( homedir );
    
    Profiler::Clear( path );

    
//    int thread_count = 8;
    int thread_count = 1;
    

    SpaceFillingCurve<3,Real> C ( thread_count );
    
    
    return 0;
}
