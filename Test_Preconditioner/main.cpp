#include <filesystem>
#include <iostream>

//#define TOOLS_DEBUG

#define TOOLS_ENABLE_PROFILER
#define REPULSOR_USE_AMD
//#define REPULSOR_USE_METIS

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

    
    int thread_count = 1;
    
    
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
    
    M.SetH1SolverID( 1 ); // Try to use AMD reordering if available.

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

    constexpr Int NRHS = 6;
    
    const Int ldB = NRHS;
    const Int ldX = NRHS;
    const Int ldY = NRHS;
    
    Tensor2<ExtReal,Int> B ( M.VertexCount(), ldB );
    Tensor2<ExtReal,Int> X ( M.VertexCount(), ldX );
    Tensor2<ExtReal,Int> Y ( M.VertexCount(), ldY );
    
    B.SetZero();
    X.SetZero();
    
    dump( B.Dimension(0) );
    dump( B.Dimension(1) );
    
    tic("tpe.Differential");
    tpe.Differential( M, Real(1), Real(0), &B[0][0], ldB );
    toc("tpe.Differential");
    
    
    dump( ArrayToString(B.data(), {4,6}) );
    print("");
    
    dump(&B[0][3]);
    dump(&B[M.VertexCount()-1][6]);
    dump(&(B.data()[M.VertexCount() * 6]));
    
    dump(B[M.VertexCount()-1][5]);
    
    
    
    dump(ldB);
    
    tic("tpe.Differential");
    tpe.Differential( M, Real(2), Real(0), &B[0][3], ldB );
    toc("tpe.Differential");
    
    dump( ArrayToString(B.data(), {4,6}) );
    print("");
    
    
    
    
//    tpm.MultiplyPreconditioner( M,
//        Real(1), &B[0][0], ldB,
//        Real(0), &X[0][0], ldX,
//        3
//    );
//    
//    dump( ArrayToString(X.data(), {4,6}) );
//    
//    dump( X.MaxNorm() );
//    print("");
    
    tpm.MultiplyPreconditioner( M,
        Real(1), &B[0][3], ldB,
        Real(0), &X[0][3], ldX,
        3
    );
    
    dump( ArrayToString(X.data(), {4,6}) );
    
    dump( X.MaxNorm() );
    print("");
    
    
//    X.SetZero();
//    
//    tpm.MultiplyMetric( M,
//        Real(1), &B[0][0], ldB,
//        Real(0), &X[0][0], ldX,
//        6
//    );
//    
//    dump( ArrayToString(X.data(), {4,6}) );
//    
//    dump( X.MaxNorm() );
//    print("");
    
    // Computing X = alpha * A^{-1} . B + beta * X using CG method.
    const Int  max_iter  = 100;
    const Real tolerance = 0.0001;
    const Real alpha = 1;
    const Real beta  = 0;

    tic("Solving for gradient");
    tpm.Solve( M, 
        alpha, &B[0][0], ldB,
        beta,  &X[0][0], ldX,
        3, max_iter, tolerance
    );
    toc("Solving for gradient");
    
    dump( B.CountNaNs() );
    dump( X.CountNaNs() );
    
    dump(tpm.CG_IterationCount());
    dump(tpm.CG_RelativeResiduals());
    
    print("");
    
    tic("Solving for gradient");
    tpm.Solve( M,
        alpha, &B[0][3], ldB,
        beta,  &X[0][3], ldX,
        3, max_iter, tolerance
    );
    toc("Solving for gradient");
    
    dump( B.CountNaNs() );
    dump( X.CountNaNs() );
    
    dump(tpm.CG_IterationCount());
    dump(tpm.CG_RelativeResiduals());
    
    print("");
    
    const Real alpha_inv = Real(1) / alpha;
    const Real zero      = 0;
    
    tpm.MultiplyMetric( M,
        alpha_inv, &X[0][0], ldX,
        zero,      &Y[0][0], ldY,
        3
    );
    
    tpm.MultiplyMetric( M,
        alpha_inv, &X[0][3], ldX,
        zero,      &Y[0][3], ldY,
        3
    );

    
    dump(B.MaxNorm());
    
    Y-=B;
    
    dump(Y.MaxNorm());
    
    return 0;
}
