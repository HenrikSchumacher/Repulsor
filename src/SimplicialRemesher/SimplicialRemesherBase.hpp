#pragma once

namespace Repulsor
{
    template<typename Real_, typename Int_, typename SReal_, typename ExtReal_>
    class SimplicialRemesherBase
    {
        
    public:
        
        using Real    = Real_;
        using Int     = Int_;
        using SReal   = SReal_;
        using ExtReal = ExtReal_;
        
        using Vertex_T   = Int;
        using Edge_T     = Int;
        using Simplex_T  = Int;
        
        using MeshBase_T = SimplicialMeshBase<Real,Int,SReal,ExtReal>;

        SimplicialRemesherBase() = default;

        virtual ~SimplicialRemesherBase() = default;
        
        virtual Int VertexCount() const = 0;

        virtual Int EdgeCount() const = 0;
        
        virtual Int SimplexCount() const = 0;
        
        virtual void Compress() = 0;

        virtual void LoadMesh( const MeshBase_T & M ) = 0;
        
        virtual void LoadMesh( const MeshBase_T & M, const Tensor2<Real,Int> & V_data_ ) = 0;

        virtual void LoadMesh( const MeshBase_T & M, const Real * V_data_, const Int data_dim ) = 0;
        
        virtual void LoadMesh_External( const MeshBase_T & M, const ExtReal * V_data_, const Int data_dim ) = 0;
        
        virtual std::unique_ptr<MeshBase_T> CreateMesh() = 0;

        virtual Tensor2<Real,Int> VertexData() = 0;

        virtual Int SplitEdges( const Edge_T * const e_list, const Int n ) = 0;
        
        virtual Int CollapseEdges( const Edge_T * const e_list, const Int n ) = 0;
        
        virtual Int FlipEdges( const Edge_T * const e_list, const Int n, const bool check_Delaunay = false ) = 0;

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
            return "SimplicialRemesherBase<"+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    }; // class SimplicialRemesherBase
    
} // namespace Repulsor
