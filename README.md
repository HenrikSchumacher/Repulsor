# Repulsor

by Henrik Schumacher 

This header-only _C++_ library allows you to work with the so-called _generalized tangent-point energy_ of curves and surface in 2-, 3-, and 4-dimensional Euclidean space.

The library provide:

- a simple data structure for simplicial meshes --> `SimplicialMesh`;
- an implementation of the multipole-accelerate discrete tangent-point energy and its derivative (it can also be used to compute Coulomb-like potentials when setting the exponent of the numerator to 0)  --> `TangentPointEnergy0`;
- an implementation of the multipole-accelerate discrete tangent-point energy between a manifold and an obstacle manifold (of potentially different dimens)  --> `TangentPointObstacleEnergy`;
- facilities to compute the matrix-vector product with an (almost) Riemannian metric --> `TangentPointMetric0`;
- an interative solver for this Riemannian metric that can be used to compute a good Sobolev gradient for steepest descent --> `TangentPointMetric0::Solve`;
- a line search implementation for finding collision-free step sizes in the direction into a given displacement field along a curve or surface --> `SimplicialMesh::MaximumSafeStepSize`.

Moreover, the library also provides 
- an implementation of the naive (all-pairs) discrete tangent-point energy and its derivative --> `TangentPointEnergy_AllPairs`;
- a simple remesher based in edge-splits and edge-contractions --> `SimplicialRemesher` (there are certainly more sophisticated and more performant remeshers out there; but this one gets the job done, too).

For code examples please see the directories `Example` and `Example2`.  


# Download/installation

Either clone with

    git clone --recurse-submodules

or clone as usual and then run the following to connect all submodules to their repos.

    git submodule update --init --recursive
    

Pull changes from the remote repositories of any submodule by executing

    git pull && git submodule foreach --recursive "git checkout main && git pull"

    
# Usage

## Within the _C++_ code

_Repulsor_ is header-only, so you do not have to precompile anything and thus you also find no makefile here. Just include

    #include "Repulsor.hpp"
    
and tell your compiler where to find it: either put the parent directory of this file onto the search path or submit it with the `-I` compiler command (_clang_ or _gcc_ ).

However, you also need implementations of _CBLAS_ and _LAPACK_. There are various choices for these and they may have a heavy impact on the performance. That's why I did not hard-code any implementation into `Repulsor.hpp`. See the instructions below for details. I recommend to use Intel oneMKL for maximal performance on Intel machines and Apple Accelerate on Apple Silicon. I think _AMD AOCL-BLAS_ should be the best option on AMD hardware, but I have no experience to back that up. OpenBLAS seems to provide also a very good performance (at least on Apple Silicon). And it has the advantage of being the most portable implementation.


## CBLAS/LAPACK

Most CBLAS and LAPACK implementations come with files `cblas.h` and `lapack.h`. In that case, just put

    #include <cblas.h>
    #include <lapack.h>

into your _C++_ code _before_ you include `Repulsor.hpp`. Of course, your path variables or compiler flags should hint the compiler to these files. And you also have to link the according libraries. But the actual library name and so the linker command may differ from implementation to implementation, so please refer to their documentations for details.

If you use Intel oneMKL, then you would rather have to use the following:

    #include <mkl_cblas.h>
    #include <mkl_lapack.h>
    
If in doubt, then better consult the [Link Line Advisor](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html).

Under macos you can use the _Accelerate_ framework. Make sure to use the most recent API by inserting this

    #define ACCELERATE_NEW_LAPACK
    #include <Accelerate/Accelerate.h>
    
into your code before you include `Repulsor.hpp`. Then you also have to issue the compile option `-framework Accelerate`.


## Compiling and linking

The code should be largely portable, but for a lack of test sytems, I tested it only in the following configurations:

- under macos on Apple Silicon with _Apple clang_ as compiler and 
- under Linux with x86 processors with _gcc_ as compiler. 

This is why I give instructions only for _clang_ and _gcc_.
In particular, I have no experience with the Microsoft Visual C++ Compiler. So beware that the options and flags there might come under different names. 


### Compiler options

_Repulsor_ uses several _C++ 20_ features, so make sure to use a compatible _C++_ implementation, e.g., by issueing the compiler option `-std=c++20`.

Older versions of _Repulsor_ employd _OpenMP_ for parallelization. That turned out to be a maintenance hell because of various incompatible implementations (`libomp`, `libiomp5`, `libgomp`,...)  that might be used by the calling program. Nowadays _Repulsor_ simply spawns a couple of `std::thread`s when needed. So make sure to use the `-pthread` option.
Btw., _Repulsor_ does not use any sophisticated thread pools like _OpenMP_ does, so I do not expect any problems with other parallelization frameworks. (Curiously, this did not come with _any_ performance penalty.)

Optimization flags like `-O3` or even `-Ofast` are certainly a good idea. I also found that using `-flto` can make a measurable difference.

With _clang_ as compiler you also have to issue `-fenable-matrix` to enable the clang matrix extension.

If you use _Accelerate_, then you also have to add the compiler option `-framework Accelerate`.

### Linker flags

Only the linker flags needed for CBLAS and LAPACK; nothing if you use _Accelerate_. 
