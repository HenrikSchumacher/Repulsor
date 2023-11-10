#include <iostream>

#include <sys/types.h>
#include <pwd.h>

//#define REMESHER_VERBATIM

#define TOOLS_ENABLE_PROFILER // enable profiler
//#define TOOLS_DEBUG

#define LAPACK_DISABLE_NAN_CHECK
#define ACCELERATE_NEW_LAPACK
#include <Accelerate/Accelerate.h>

//#define LAPACK_DISABLE_NAN_CHECK
//#include <cblas.h>
//#include <lapack.h>

#include "../Repulsor.hpp"
#include "../submodules/Tensors/MyBLAS.hpp"
//#include "../submodules/Tensors/Sparse.hpp"
//#include "../submodules/Tensors/ConjugateGradient.hpp"
//#include "../submodules/Tensors/GMRES.hpp"

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

    SimplicialMesh_Factory<MeshBase_T,2,2,3,3> mesh_factory;
    
    tic("Initializing mesh");
    std::unique_ptr<MeshBase_T> M_ptr = mesh_factory.Make_FromFile( path + name, thread_count );
    
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
    M.GetClusterTree();           // Not necessary. Will automatically called by all routines that require it.
    toc("Creating ClusterTree");

    tic("Creating BlockClusterTree");
    M.GetBlockClusterTree();      // Not necessary. Will automatically called by all routines that require it.
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
    
//    X.Read(B.data());

    X.SetZero();
    
    tic("tpe.Differential");
    tpe.Differential( M, B.data() );
    toc("tpe.Differential");
    
    
    tpm.MultiplyMetric( M, 1., B.data(), 0, Y.data(), NRHS );
    
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
