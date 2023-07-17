#include <iostream>

#include <sys/types.h>
#include <pwd.h>

// We have to toggle which domain dimensions and ambient dimensions shall be supported by runtime polymorphism before we load Repulsor.hpp
// You can activate everything you want, but compile times might increase substatially.

//#define NDEBUG
//#define TOOLS_DEBUG

#define TOOLS_ENABLE_PROFILER // enable profiler


#define LAPACK_DISABLE_NAN_CHECK
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>

//#define LAPACK_DISABLE_NAN_CHECK
//#include <cblas.h>
//#include <lapack.h>

#include "../Repulsor.hpp"


using namespace Repulsor;
using namespace Tensors;
using namespace Tools;

using Int     = Int64;
using LInt    = Int64;
using ExtInt  = Int64;

using Real    = Real64;
using SReal   = Real64;
using ExtReal = Real64;


int main(int argc, const char * argv[])
{
    const char * homedir = getenv("HOME");

    if( homedir == nullptr)
    {
        homedir = getpwuid(getuid())->pw_dir;
    }
    std::string home ( homedir );
    
    Profiler::Clear( home );

    std::string path = home + "/github/Repulsor/CoulombPotential/";
    std::string name = "Sphere_Points.txt";
    
    
    using Mesh_T      = SimplicialMesh<0,3,Real,Int,LInt,SReal,ExtReal>;
    using MeshBase_T  = SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>;

    const Int thread_count = 1;

    SimplicialMesh_Factory<MeshBase_T,0,0,3,3> mesh_factory;

    tic("Initializing mesh");

    std::unique_ptr<MeshBase_T> M_ptr = mesh_factory.Make_FromFile( path + name, thread_count );

    Tensor1<Real,Int> vertex_charges ( M_ptr->VertexCount(), 1. );

    Mesh_T M (
        M_ptr->VertexCoordinates().data(), M_ptr->VertexCount(),  false,
        M_ptr->Simplices().data(),         M_ptr->SimplexCount(), false,
        vertex_charges.data(),
        M_ptr->ThreadCount()
    );

    M_ptr = nullptr;

    dump(M.ThreadCount());

    toc("Initializing mesh");

    print("");

    // Some quite decent settings for 2-dimensional surfaces.
    M.cluster_tree_settings.split_threshold                        =  2;
    M.cluster_tree_settings.thread_count                           =  0; // take as many threads as there are used by SimplicialMesh M
    M.block_cluster_tree_settings.far_field_separation_parameter   =  0.25;
    M.adaptivity_settings.theta                                    = 10.0;

    tic("GetClusterTree");
    M.GetClusterTree();           // Not necessary. Will automatically called by all routines that require it.
    toc("GetClusterTree");

    tic("GetBlockClusterTree");
    M.GetBlockClusterTree();      // Not necessary. Will automatically called by all routines that require it.
    toc("GetBlockClusterTree");

    const Real alpha = 1;

//    CoulombEnergy<Mesh_T> E( alpha );

    TangentPointEnergy<Mesh_T> E( 0, alpha );

    ExtReal en;
    // Mesh_T::CotangentVector_T is Tensor2<REAL,INT> in this case. It is a simple container class for heap-allocated matrices.
    Mesh_T::CotangentVector_T diff;
//
    tic("tpe.Energy(M)");
    en = E.Value(M);
    toc("tpe.Energy(M)");

    dump(en);

    tic("tpe.Differential(M)");
    diff = E.Differential(M);
    toc("tpe.Differential(M)");
    
    return 0;
}
