#pragma once

namespace Repulsor
{
    template<typename Real_, typename Int_>
    class SimplicialRemesherBase
    {
        
        ASSERT_REAL(Real_);
        ASSERT_INT(Int_);
        
    public:
        
        using Real    = Real_;
        using Int     = Int_;
        
        using Vertex_T   = Int;
        using Edge_T     = Int;
        using Simplex_T  = Int;

        SimplicialRemesherBase() = default;

        virtual ~SimplicialRemesherBase() = default;
        
        virtual Int VertexCount() const = 0;

        virtual Int EdgeCount() const = 0;
        
        virtual Int SimplexCount() const = 0;
        
        virtual void Compress() = 0;

        virtual void LoadMesh(
            cptr<Real> V_coords_ ,  const Int vertex_count_,
            cptr<Int>  simplices_,  const Int simplex_count_,
            cptr<Real> V_data_,     const Int V_data_dim_,
            const Int thread_count_ = 1
        ) = 0;
        
//        virtual void LoadMesh( cref<MeshBase_T> M ) = 0;
//
//        virtual void LoadMesh( cref<MeshBase_T> M, cref<Tensor2<Real,Int>> V_data_ ) = 0;

//        virtual void LoadMesh( cref<MeshBase_T> M, cptr<Real> V_data_, const Int data_dim ) = 0;
//        
//        virtual void LoadMesh_External( cref<MeshBase_T> M, cptr<ExtReal> V_data_, const Int data_dim ) = 0;
        
//        virtual std::unique_ptr<MeshBase_T> CreateMesh() = 0;

        virtual cref<Tensor2<Real,Int>> VertexCoordinates() = 0;
        
        virtual cref<Tensor2<Real,Int>> VertexData() = 0;
        
        virtual cref<Tensor2<Int,Int>> Simplices() = 0;
        
        virtual Int DomDim() const = 0;
        
        virtual Int AmbDim() const = 0;
        
        virtual Int SplitEdges( cptr<Edge_T> e_list, const Int n ) = 0;
        
        virtual Int SplitEdges( cref<std::vector<Int>> e_list ) = 0;
        
        virtual Int CollapseEdges( cptr<Edge_T> e_list, const Int n ) = 0;
        
        virtual Int CollapseEdges( cref<std::vector<Int>> e_list ) = 0;

        
        
        virtual Int FlipEdges( cptr<Edge_T> e_list, const Int n, const bool check_Delaunay = false ) = 0;
        
        virtual Int FlipEdges( cref<std::vector<Int>> e_list, const bool check_Delaunay = false ) = 0;

        virtual bool UnifyEdgeLengths(
            const Real collapse_threshold,
            const Real split_threshold,
            const Int  max_iter = 100
        ) = 0;
        
        virtual Int DelaunayFlip( const Int max_iter = 100 ) = 0;
        
        virtual void SelfCheck() = 0;
        
    public:
        
        virtual std::string ClassName() const
        {
            return std::string("SimplicialRemesherBase<")+TypeName<Real>+","+TypeName<Int>+">";
        }
        
    }; // class SimplicialRemesherBase
    
} // namespace Repulsor
