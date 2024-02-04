#include <filesystem>
#include <iostream>

/// This is a simple program meant to make you familiar with how to
/// - generate a simplicial mesh from data _from file_;
/// - compute the tangent-point energy (and derivative) of such a mesh;
/// - perform matrix-vector operations with the tangent-point metric;
/// - compute the gradient of the tangent-point energy w.r.t. that metric; and
/// - find a collision-free step size for a displacement field along the manifold.


/// Enable the shipped profiler.
#define TOOLS_ENABLE_PROFILER

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

#include "../Repulsor.hpp"

using namespace Tools;

using Real = double;
using Int  = int;
using LInt = std::size_t;


#include "../Meshes/TwoSpheres.hpp"

const double * vertex_coords = &vertex_coordinates[0][0];
const int    * triangles     = &simplices[0][0];

const int triangle_count = simplex_count;



int main(void)
{
    constexpr int n = 2;
    constexpr int m = 3;
    
    double p = 6;

    using Mesh_T = Repulsor::SimplicialMesh<n,m,double,int>;
    using TPE_T  = Repulsor::TangentPointEnergy<Mesh_T>;
    
    Mesh_T mesh ( vertex_coords, vertex_count, triangles, triangle_count );
    
    TPE_T tpe ( p );
    
    double energy = tpe.Value( mesh );
    
    std::cout << "tangent-point energy = " << energy << "." << std::endl;
    
    return 0;
}
