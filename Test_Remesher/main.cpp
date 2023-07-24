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

using namespace Repulsor;
using namespace Tensors;
using namespace Tools;


// We have to toggle which domain dimensions and ambient dimensions shall be supported by runtime polymorphism before we load Repulsor.hpp
// You can activate everything you want, but compile times might increase substatially.
using Int     = Int32;
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
//    std::string name = "TorusMesh_00038400T.txt";
    std::string name = "Spot_00005856T.txt";
    
    
    using MeshBase_T = SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>;
    using Mesh_T     = SimplicialMesh<2,3,Real,Int,LInt,SReal,ExtReal>;
    
    SimplicialMesh_Factory<MeshBase_T,2,2,3,3> mesh_factory;
    
    
    tic("Initializing mesh");
    std::unique_ptr<MeshBase_T> M_ptr = mesh_factory.Make_FromFile( path + name, thread_count );
    
    auto & M = *M_ptr;

    dump(M.ThreadCount());

    toc("Initializing mesh");

    print("");

    
    print("");
    print("Testing remesher.");
    print("");
    
    std::unique_ptr<Mesh_T::RemesherBase_T> R = M_ptr->CreateRemesher();
    
    Tensor1<Real,Int> squared_edge_lengths = R->SquaredEdgeLengths();
    
    Sort( squared_edge_lengths.begin(), squared_edge_lengths.end() );
    
    const Real lower_bound = std::sqrt( squared_edge_lengths[Int(Real(R->EdgeCount()-1) * Real(0.1))] );
    const Real upper_bound = std::sqrt( squared_edge_lengths[Int(Real(R->EdgeCount()-1) * Real(0.9))]  );
    
    valprint("Minimum edge length            ", std::sqrt(squared_edge_lengths.First()) );
    valprint("Lower bound used for collapsing", lower_bound );
    valprint("Upper bound used for splitting ", upper_bound );
    valprint("Maximum edge length            ", std::sqrt(squared_edge_lengths.Last()) );
    
    const Int  collapse_split_iter = 4;

    print("");
    
    
    //Apply collapse_split_iter rounds of edge collapse and edge split.
    tic("UnifyEdgeLengths");
    R->UnifyEdgeLengths( lower_bound, upper_bound, collapse_split_iter );
    toc("UnifyEdgeLengths");
    
    print("");
    
    const Int flip_iter = 3;
    
    print("");
    
    //Apply flip_max_iter rounds of edge flips to improve triangle shapes.
    tic("DelaunayFlip");
    R->DelaunayFlip( flip_iter );
    toc("DelaunayFlip");
    
    const Int smooth_iter = 10;
    
    print("");
    
    // Apply smooth_iter iterations of tangential smoothening to further improve mesh regularity.
    tic("TangentialSmoothing");
    R->TangentialSmoothing( smooth_iter );
    toc("TangentialSmoothing");
    
    std::unique_ptr<MeshBase_T> N_ptr = mesh_factory.Make(
        R->VertexCoordinates().data(), R->VertexCount(),  R->AmbDim(),   false,
        R->Simplices().data(),         R->SimplexCount(), R->DomDim()+1, false,
        thread_count
    );
    
    R->SelfCheck();
    
    R->Compress();
    
    squared_edge_lengths = R->SquaredEdgeLengths();
    
    valprint("Minimum edge length            ", std::sqrt(squared_edge_lengths.Min()) );
    valprint("Maximum edge length            ", std::sqrt(squared_edge_lengths.Max()) );
    
    return 0;
}
