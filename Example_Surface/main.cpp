#include <filesystem>
#include <iostream>

/// This is a simple program meant to make you familiar with how to
/// - generate a simplicial mesh from data _from file_;
/// - compute the tangent-point energy (and derivative) of such a mesh;
/// - perform matrix-vector operations with the tangent-point metric;
/// - compute the gradient of the tangent-point energy w.r.t. that metric; and
/// - find a collision-free step size for a displacement field along the manifold.


#ifdef __APPLE__
/// Use these while on a mac. Don't forget to issue the compiler flag `-framework Accelerate`.
///
    #include "../submodules/Tensors/Accelerate.hpp"
#else
/// This should work for OpenBLAS.
    #include "../submodules/Tensors/OpenBLAS.hpp"
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
