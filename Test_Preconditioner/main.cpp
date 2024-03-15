#include <filesystem>
#include <iostream>

//#define TOOLS_DEBUG

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


// We have to toggle which domain dimensions and ambient dimensions shall be supported by runtime polymorphism before we load Repulsor.hpp
// You can activate everything you want, but compile times might increase substatially.
using Int     = Int64;
using LInt    = Int64;
using ExtInt  = Int64;

using Real    = Real64;
using SReal   = Real64;
using ExtReal = Real64;


int main(void)
{
    print("");
    print("###############################################################");
    print("###              Test program preconditioner               ####");
    print("###############################################################");
    print("");
    
    Profiler::Clear( getenv("HOME") );

    
    int thread_count = 8;
    
    
    using MeshBase_T = SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>;
    using Mesh_T     = SimplicialMesh<2,3,Real,Int,LInt,SReal,ExtReal>;

    SimplicialMesh_Factory<MeshBase_T,2,2,3,3> mesh_factory;


    
    tic("Initializing mesh");
    
    std::filesystem::path this_file { __FILE__ };
    std::filesystem::path repo_dir  = this_file.parent_path().parent_path();
    std::filesystem::path mesh_file = repo_dir / "Meshes" / "TorusMesh_00038400T.txt";
    
    std::unique_ptr<MeshBase_T> M_ptr = mesh_factory.Make_FromFile(
        mesh_file.c_str(), thread_count
    );
    
    auto & M = *M_ptr;
    
    // Some quite decent settings for 2-dimensional surfaces.
    M.cluster_tree_settings.split_threshold                        =  2;
    M.cluster_tree_settings.thread_count                           =  0; // take as many threads as there are used by SimplicialMesh M
    M.block_cluster_tree_settings.far_field_separation_parameter   =  0.25;
    M.adaptivity_settings.theta                                    = 10.0;
    

    dump(M.ThreadCount());

    toc("Initializing mesh");
    
    const Real q = 6;
    const Real p = 12;
    
    TangentPointEnergy0<Mesh_T> tpe (q,p);
    TangentPointMetric0<Mesh_T> tpm (q,p);

    
    tic("Creating ClusterTree");
    M.GetClusterTree();
    toc("Creating ClusterTree");

    tic("Creating BlockClusterTree");
    M.GetBlockClusterTree();
    toc("Creating BlockClusterTree");

    print("");

    
    double en;

    tic("tpe.Energy");
    en = tpe.Value(M);
    toc("tpe.Energy");

    dump(en);

    
    print("");

    constexpr Int NRHS = 3;
    
    Tensor2<ExtReal,Int> B  ( M.VertexCount(), NRHS );
    Tensor2<ExtReal,Int> X  ( M.VertexCount(), NRHS );
    Tensor2<ExtReal,Int> Y  ( M.VertexCount(), NRHS );

    X.SetZero();
    
    tic("tpe.Differential");
    tpe.Differential( M, B.data() );
    toc("tpe.Differential");
    
    print("");
    
    const Int  max_iter  = 100;
    const Real tolerance = 0.00001;
    
    tic("Solving for gradient");
    tpm.Solve(M, B.data(), X.data(), NRHS, max_iter, tolerance );
    toc("Solving for gradient");
    
    dump(tpm.CG_IterationCount());
    dump(tpm.CG_RelativeResiduals());
    
    tpm.MultiplyMetric( M, 1., X.data(), 0, Y.data(), NRHS );

    
    dump(B.MaxNorm());
    
    Y-=B;
    
    dump(Y.MaxNorm());
    
    return 0;
}
