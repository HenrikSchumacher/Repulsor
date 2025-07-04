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
/// This should work for OpenBLAS. Don't forget to ussie the compiler flag `-lopenblas`.
    #include "../submodules/Tensors/OpenBLAS.hpp"
#endif

#include "../Repulsor.hpp"

using namespace Repulsor;
using namespace Tensors;
using namespace Tools;

// We have to toggle which domain dimensions and ambient dimensions shall be supported by runtime polymorphism before we load Repulsor.hpp
// You can activate everything you want, but compile times might increase substatially.
//using Int     = Int64;
using Int     = Int32;  // Used for the indices of the mesh triangles. 32 bits should suffice.
using LInt    = Int64;  // We need this 64 for the row pointers of big sparse matrices
using ExtInt  = Int64;

using Real    = Real64; // We need full precision here, definitely.
using SReal   = Real64; // We need full precision here, definitely.
using ExtReal = Real64;


int main(void)
{
    print("");
    print("###############################################################");
    print("###              Test program preconditioner               ####");
    print("###############################################################");
    print("");
    
    Profiler::Clear();

    int thread_count = 8;
    
    using MeshBase_T = SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>;
    using Mesh_T     = SimplicialMesh<2,3,Real,Int,LInt,SReal,ExtReal>;

    SimplicialMesh_Factory<MeshBase_T,2,2,3,3> mesh_factory;
    
    tic("Initializing mesh");
    
    std::filesystem::path this_file { __FILE__ };
    std::filesystem::path repo_dir  = this_file.parent_path().parent_path();
//    std::filesystem::path mesh_file = repo_dir / "Meshes" / "TorusMesh_00002400T.txt";
//    std::filesystem::path mesh_file = repo_dir / "Meshes" / "TorusMesh_00038400T.txt";
    
    std::filesystem::path mesh_file = "/Users/Henrik/github/BAEMM/Meshes/Spot_00093696T.txt";
    
//    std::filesystem::path mesh_file = "/Users/Henrik/github/BAEMM/Meshes/Spot_00023424T.txt";
    
    std::unique_ptr<MeshBase_T> M_ptr = mesh_factory.Make_FromFile(
        mesh_file, thread_count
    );
    
    auto & M = *M_ptr;
    
    // Some quite decent settings for 2-dimensional surfaces.
    M.cluster_tree_settings.split_threshold                        =  2;
    M.cluster_tree_settings.thread_count                           =  0; // take as many threads as there are used by SimplicialMesh M
    M.block_cluster_tree_settings.far_field_separation_parameter   =  0.25;
    M.adaptivity_settings.theta                                    = 10.0;
    
    M.SetH1SolverID( 1 ); // Try to use AMD reordering if available.
    
    // You can adapt the weights of the two terms in the H^1 metric.
    // Typically, you should set H1StiffnessWeight to 1 and H1MassWeight to a relatively small value. However, if the latter is too small, the solver might become instable.
    M.SetH1StiffnessWeight( 1 );
    M.SetH1MassWeight( 0.1 );

    TOOLS_DUMP(M.ThreadCount());

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

    TOOLS_DUMP(en);

    
    print("");

//    constexpr Int NRHS = 3;
//    const Int print_rows = 4;
    const Int ld   = 6;
    const Int nrhs = 3;
    
    const Int ldB = ld;
    const Int ldX = ld;
    const Int ldY = ld;
    const Int ldZ = ld;
    
    Tensor2<ExtReal,Int> B ( M.VertexCount(), ldB );
    Tensor2<ExtReal,Int> X ( M.VertexCount(), ldX );
    Tensor2<ExtReal,Int> Y ( M.VertexCount(), ldY );
    Tensor2<ExtReal,Int> Z ( M.VertexCount(), ldZ );
    
    B.SetZero();
    X.SetZero();
    Y.SetZero();
    
    TOOLS_DUMP( B.Dim(0) );
    TOOLS_DUMP( B.Dim(1) );
    
    print("");
    print("Experiment 1: Solve A.X = B for X and check the perconditioner norm of B - A.X ");
    print("");
    
    tic("tpe.Differential");
    tpe.Differential( M, Real(1), Real(0), &B[0][0], ldB );
    toc("tpe.Differential");
    
    
    print("");
//    print("tpe.Differential( M, Real(0.5), Real(0), &B[0][nrhs], ldB );");

    tic("tpe.Differential");
    tpe.Differential( M, Real(0.5), Real(0), &B[0][nrhs], ldB );
    toc("tpe.Differential");
    
//    TOOLS_DUMP( ArrayToString( B.data(), {print_rows,ld}, 8 ) );
    print("");
    

    // Computing X = alpha * A^{-1} . B + beta * X using CG method.
    const Int  max_iter  = 50;
    const Real tolerance = 0.00001;
    const Real alpha = 1;
    const Real beta  = 0;
    

    tic("Solving for gradient");
    tpm.Solve( M, 
        alpha, &B[0][0   ], ldB,
        beta,  &X[0][0   ], ldX,
        nrhs, max_iter, tolerance
    );
    toc("Solving for gradient");
    
    TOOLS_DUMP( B.CountNaNs() );
    TOOLS_DUMP( X.CountNaNs() );
    
    TOOLS_DUMP(tpm.Solver_IterationCount());
    TOOLS_DUMP(tpm.Solver_RelativeResiduals());
    
    print("");
    
    tic("Solving for gradient");
    tpm.Solve( M,
        alpha, &B[0][nrhs], ldB,
        beta,  &X[0][nrhs], ldX,
        nrhs, max_iter, tolerance
    );
    toc("Solving for gradient");
    
    TOOLS_DUMP( B.CountNaNs() );
    TOOLS_DUMP( X.CountNaNs() );
    
    TOOLS_DUMP(tpm.Solver_IterationCount());
    TOOLS_DUMP(tpm.Solver_RelativeResiduals());
    
    // Checking the result.
    
    const Real alpha_inv = Real(1) / alpha;

    Y = B;

    tpm.MultiplyMetric( M,
        alpha_inv, &X[0][0   ], ldX,
        Real(-1),  &Y[0][0   ], ldY,
        nrhs
    );
    tpm.MultiplyMetric( M,
        alpha_inv, &X[0][nrhs], ldX,
        Real(-1),  &Y[0][nrhs], ldY,
        nrhs
    );
    
    Tiny::Vector<2,Real,Int> B_prec_norms {
        tpm.PreconditionerNorm( M, &B[0][0   ], ldB, nrhs ),
        tpm.PreconditionerNorm( M, &B[0][nrhs], ldB, nrhs )
    };
    Tiny::Vector<2,Real,Int> abs_prec_errors {
        tpm.PreconditionerNorm( M, &Y[0][0   ], ldB, nrhs ),
        tpm.PreconditionerNorm( M, &Y[0][nrhs], ldB, nrhs )
    };
    
    Tiny::Vector<2,Real,Int> rel_prec_errors {
        abs_prec_errors[0] / B_prec_norms[0],
        abs_prec_errors[1] / B_prec_norms[1]
    };
    
    TOOLS_DUMP(rel_prec_errors);
    
    print("");
    print("Experiment 1: Done.");
    print("");
    
    
    print("");
    print("Experiment 2: Compute B = A.Z, then solve A.X = B for X and check the metric norm of Z - X");
    print("");
    
    M.VertexCoordinates().Write( &X[0][0   ], ldX );
    M.VertexCoordinates().Write( &X[0][nrhs], ldX );
    
    Tiny::Vector<2,Real,Int> X_true_metric_norms {
        tpm.MetricNorm( M, &X[0][0   ], ldX, nrhs ),
        tpm.MetricNorm( M, &X[0][nrhs], ldX, nrhs )
    };
    
    // B = (1/alpha) * A.X
    tpm.MultiplyMetric( M,
        alpha_inv, &X[0][0   ], ldX,
        Real(0),   &B[0][0   ], ldB,
        nrhs
    );
    tpm.MultiplyMetric( M,
        alpha_inv, &X[0][nrhs], ldX,
        Real(0),   &B[0][nrhs], ldB,
        nrhs
    );
    
    tic("Computing X = alpha * A^{-1}.B - X_true.");
    tpm.Solve( M,
        alpha,    &B[0][0   ], ldB,
        Real(-1), &X[0][0   ], ldX,
        nrhs, max_iter, tolerance
    );
    toc("Computing X = alpha * A^{-1}.B - X_true.");
    
    TOOLS_DUMP(tpm.Solver_IterationCount());
    TOOLS_DUMP(tpm.Solver_RelativeResiduals());
    
    tic("Computing X = alpha * A^{-1}.B - X_true.");
    tpm.Solve( M,
        alpha,    &B[0][nrhs], ldB,
        Real(-1), &X[0][nrhs], ldX,
        nrhs, max_iter, tolerance
    );
    toc("Computing X = alpha * A^{-1}.B - X_true.");
    
    TOOLS_DUMP(tpm.Solver_IterationCount());
    TOOLS_DUMP(tpm.Solver_RelativeResiduals());
    
    Tiny::Vector<2,Real,Int> abs_metric_errors {
        tpm.MetricNorm( M, &X[0][0   ], ldX, nrhs ),
        tpm.MetricNorm( M, &X[0][nrhs], ldX, nrhs )
    };
    Tiny::Vector<2,Real,Int> rel_metric_errors {
        abs_metric_errors[0] / X_true_metric_norms[0],
        abs_metric_errors[1] / X_true_metric_norms[1]
    };
    
    TOOLS_DUMP(rel_metric_errors);

    print("");
    print("Experiment 2: Done.");
    print("");
    
    return 0;
}
