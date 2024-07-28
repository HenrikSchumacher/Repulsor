#include <filesystem>
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


int main(void)
{
    /// Set up profiler to write to `~/Tools_Profile.tsv` and `~/Tools_Log.txt`
    Profiler::Clear();
    
    using Mesh_T      = SimplicialMesh<0,3,Real,Int,LInt,SReal,ExtReal>;
    using MeshBase_T  = SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>;

    const Int thread_count = 1;

    SimplicialMesh_Factory<MeshBase_T,0,0,3,3> mesh_factory;

    tic("Initializing mesh");

    std::filesystem::path this_file { __FILE__ };
    std::filesystem::path repo_dir  = this_file.parent_path().parent_path();
    std::filesystem::path mesh_file = repo_dir / "Meshes" / "Sphere_Points.txt";
    
    std::unique_ptr<MeshBase_T> M_ptr = mesh_factory.Make_FromFile(
        mesh_file, thread_count
    );

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
    M.GetClusterTree();           // Not necessary. Will automatically be called by all routines that require it.
    toc("GetClusterTree");

    tic("GetBlockClusterTree");
    M.GetBlockClusterTree();      // Not necessary. Will automatically be called by all routines that require it.
    toc("GetBlockClusterTree");

    const Real alpha = 1;

    TangentPointEnergy<Mesh_T> E( 0, alpha );

    ExtReal en;
    Mesh_T::CotangentVector_T diff ( M.VertexCount(), M.AmbDim() );

    tic("tpe.Differential(M)");
    en = E.Differential(M, diff.data());
    toc("tpe.Differential(M)");

    dump(en);
    
    return 0;
}
