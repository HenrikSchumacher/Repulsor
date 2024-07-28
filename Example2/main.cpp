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
    #include "../submodules/Tensors/Accelerate.hpp"
#else
/// This should work for OpenBLAS.
    #include "../submodules/Tensors/OpenBLAS.hpp"
#endif

#include "../Repulsor.hpp"

using namespace Tools;

using Int     = int;
using LInt    = std::size_t;

using Real    = double;

int main(void)
{
    /// Set up profiler to write to `~/Tools_Profile.tsv` and `~/Tools_Log.txt`
    Profiler::Clear();
    
    using Mesh_T   = Repulsor::SimplicialMeshBase<Real,Int,LInt>;
    using Energy_T = Repulsor::EnergyBase<Mesh_T>;
    using Metric_T = Repulsor::MetricBase<Mesh_T>;

    const Int thread_count = 8;
    
    /// These turned out to be good choices for the parameters of the tangent-point energy for 2-dimensional surfaces.
    const Real q = 6;
    const Real p = 12;
    
    Repulsor::SimplicialMesh_Factory<Mesh_T,2,2,3,3> mesh_factory;
    Repulsor::TangentPointEnergy0_Factory<Mesh_T,2,2,3,3> TPE_factory;
    Repulsor::TangentPointMetric0_Factory<Mesh_T,2,2,3,3> TPM_factory;
    
    std::unique_ptr<Energy_T> tpe_ptr = TPE_factory.Make(2,3,q,p);
    std::unique_ptr<Metric_T> tpm_ptr = TPM_factory.Make(2,3,q,p);
    
    const auto & tpe = *tpe_ptr;
    const auto & tpm = *tpm_ptr;
    
    tic("Initializing mesh");
    
    std::filesystem::path this_file { __FILE__ };
    std::filesystem::path repo_dir  = this_file.parent_path().parent_path();
    std::filesystem::path mesh_file = repo_dir / "Meshes" / "TorusMesh_00038400T.txt";
    
    std::unique_ptr<Mesh_T> M_ptr = mesh_factory.Make_FromFile( mesh_file, thread_count );
    
    auto & M = *M_ptr;

    toc("Initializing mesh");
    
    print("");
    
    /// Some quite decent settings for 2-dimensional surfaces.
    M.cluster_tree_settings.split_threshold                        =  2;
    /// Take as many threads as there are used by SimplicialMesh M
    M.cluster_tree_settings.thread_count                           =  0;
    /// Set the separation parameters.
    M.block_cluster_tree_settings.far_field_separation_parameter   =  0.25;
    M.adaptivity_settings.theta                                    = 10.0;

    
    Real en;
    Mesh_T::CotangentVector_T diff ( M.VertexCount(), M.AmbDim() );
    
    tic("Compute tangent-point energy and derivative");
//    en = tpe.Value(M);
    en = tpe.Differential(M, diff.data());
    toc("Compute tangent-point energy and derivative");

    dump(en);

    
    /// This is how you can get the gradient of the tangent-point energy:
    const Int  max_iter  = 100;
    const Real relative_tolerance = 0.00001;
    
    Mesh_T::TangentVector_T gradient ( M.VertexCount(), M.AmbDim() );
    
    tic("Solving for gradient");
    /// We have to solve M.AmbDim() equations at once.
    tpm.Solve(M,
        Real(1), diff.data(),     M.AmbDim(),
        Real(0), gradient.data(), M.AmbDim(),
        M.AmbDim(),
        max_iter, relative_tolerance
    );
    toc("Solving for gradient");
    
    dump(tpm.CG_IterationCount());
    dump(tpm.CG_RelativeResiduals());
    
    
    Mesh_T::TangentVector_T downward_gradient = gradient;
    
    downward_gradient*= -1.;
    
    print("");
    
    /// Given a displacement field like `downward_gradient` along the manifold `M`, you can use `MaximumSafeStepSize` to compute an (almost) maximal collision-free step size in the interval [0,1].
    tic("MaximumSafeStepSize");
    const Real t = M.MaximumSafeStepSize(downward_gradient.data(), 1.);
    
    toc("MaximumSafeStepSize");
    
    valprint("Collision-free step size",t);
    
    print("");
    
    
    print("");
    
    tic("SemiStaticUpdate");
    /// X = M.VertexCoordinates() + t * downward_gradient.
    
    Mesh_T::TangentVector_T X = downward_gradient;
    X *= t;
    X += M.VertexCoordinates();

    
    /// Update mesh without rebuilding the bounding volume hierarchy.
    /// Use this for line search!
    M.SemiStaticUpdate( X.data() );
    
    toc("SemiStaticUpdate");
    
    print("");
    
    /// Let's see whether the energy has decreased...
    
    Real en_new;

    tic("Compute tangent-point energy of updated mesh");
    en_new = tpe.Value(M);
    toc("Compute tangent-point energy of updated mesh");

    dump(en_new);
    
    dump(en_new < en);
    

    return 0;
}
