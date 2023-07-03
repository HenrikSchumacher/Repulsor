#include <iostream>

#include <sys/types.h>
#include <pwd.h>



//#define REMESHER_VERBATIM

#define TOOLS_ENABLE_PROFILER // enable profiler
//#define TOOLS_DEBUG

#define LAPACK_DISABLE_NAN_CHECK
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>
//#include <cblas.h>
//#include <lapacke.h>

#include "../Repulsor.hpp"
#include "../Tensors/MyBLAS.hpp"
#include "../Tensors/Sparse.hpp"
#include "../Tensors/ConjugateGradient.hpp"
#include "../Tensors/GMRES.hpp"

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


int main(int argc, const char * argv[])
{
    print("");
    print("###############################################################");
    print("###              Test program preconditioner               ####");
    print("###############################################################");
    print("");
    
    
    const char * homedir = getenv("HOME");

    if( homedir == nullptr)
    {
        homedir = getpwuid(getuid())->pw_dir;
    }
    std::string home ( homedir );
    
    Profiler::Clear( home );

    
    int thread_count = 8;

    std::string path = home + "/github/BAEMM/Meshes/";
    std::string name = "TorusMesh_00038400T.txt";
//    std::string name = "Spot_00005856T.txt";
    
    
    using MeshBase_T = SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>;
    using Mesh_T     = SimplicialMesh<2,3,Real,Int,LInt,SReal,ExtReal>;

    const Real q = 6;
    const Real p = 12;
    
    const Real s = (p - 2) / q;
    
    SimplicialMesh_Factory<MeshBase_T,2,2,3,3> mesh_factory;
    
    TangentPointEnergy0<Mesh_T>       tpe        (q,p);
    TangentPointMetric0<Mesh_T>       tpm        (q,p);
    PseudoLaplacian    <Mesh_T,false> pseudo_lap (2-s);
    
    tic("Initializing mesh");
    std::unique_ptr<MeshBase_T> M_ptr = mesh_factory.Make_FromFile( path + name, thread_count );
    
    auto & M = *M_ptr;

    dump(M.ThreadCount());

    toc("Initializing mesh");

    // Some quite decent settings for 2-dimensional surfaces.
    M.cluster_tree_settings.split_threshold                        =  2;
    M.cluster_tree_settings.thread_count                           =  0; // take as many threads as there are used by SimplicialMesh M
    M.block_cluster_tree_settings.far_field_separation_parameter   =  0.25;
    M.adaptivity_settings.theta                                    = 10.0;

    
    tic("Creating ClusterTree");
    M.GetClusterTree();           // Not necessary. Will automatically called by al routines that require it.
    toc("Creating ClusterTree");

    tic("Creating BlockClusterTree");
    M.GetBlockClusterTree();      // Not necessary. Will automatically called by al routines that require it.
    toc("Creating BlockClusterTree");

    print("");

    
    double en;

    tic("tpe.Energy(M)");
    en = tpe.Value(M);
    toc("tpe.Energy(M)");

    dump(en);

    
    print("");

    
    
    constexpr Int NRHS = 3 * 1;
    
    Tensor2<Real,Int> B_buffer  ( M.VertexCount(), NRHS );
    Tensor2<Real,Int> X_buffer  ( M.VertexCount(), NRHS );
    Tensor2<Real,Int> Y_buffer  ( M.VertexCount(), NRHS );
    Tensor2<Real,Int> Z_buffer  ( M.VertexCount(), NRHS );
    
    mut<Real> B  = B_buffer.data();
    mut<Real> X  = X_buffer.data();
    mut<Real> Y  = Y_buffer.data();
    mut<Real> Z  = Z_buffer.data();
    
    tic("tpe.Differential(M)");
    tpe.Differential(M, B );
    toc("tpe.Differential(M)");

    print("");
    
    // The operator for the metric.
    auto A = [&]( ptr<Real> X, mut<Real> Y )
    {
        tpm.MultiplyMetric( M, Scalar::One<Real>, X, Scalar::Zero<Real>, Y, NRHS );
    };

    auto H1_metric = M.H1Metric(1,1);
    
    Permutation<Int> perm (
        M.NestedDissectionOrdering().data(), M.VertexCount(), Inverse::False, thread_count
    );
    
    Sparse::CholeskyDecomposition<Real,Int,LInt> S (
        H1_metric.Outer().data(), H1_metric.Inner().data(), std::move(perm)
    );

    S.SymbolicFactorization();

    S.NumericFactorization( H1_metric.Values().data(), Scalar::Zero<Real> );
    
    // The operator for the preconditioner.
    auto P = [&]( ptr<Real> X, mut<Real> Y )
    {
        S.Solve<Parallel>( X, Y, NRHS );
        pseudo_lap.MultiplyMetric( M, Scalar::One<Real>, Y, Scalar::Zero<Real>, Z, NRHS );
        S.Solve<Parallel>( Z, Y, NRHS );
    };

    print("");
    
    tic("P (first)");
    P(B,X);
    toc("P (first)");
    
    tic("A (first)");
    A(X,Y);
    toc("A (first)");
    
    print("");
    
    tic("P (second)");
    P(B,X);
    toc("P (second)");
    
    tic("A (second)");
    A(X,Y);
    toc("A (second)");
    
    print("");
    
    const Int max_iter = 30;
    
    ConjugateGradient<NRHS,Real,Int> CG ( M.VertexCount(), max_iter, thread_count );
    
    tic("CG");
    CG( A, P, B, NRHS, X, NRHS, 0.0001 );
    toc("CG");
    
    print("");
    
    dump(CG.IterationCount());
    dump(CG.RelativeResiduals());
    
    print("");
    
    A(X,Y);    
    Y_buffer -= B_buffer;
    
    dump( Y_buffer.MaxNorm() / B_buffer.MaxNorm() );
    
    return 0;
}
