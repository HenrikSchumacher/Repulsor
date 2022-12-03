#pragma once

#define CLASS SimplicialRemesherBase

namespace Repulsor
{
    template<typename Real_, typename Int_, typename SReal_, typename ExtReal_>
    class CLASS
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

        CLASS() = default;

        virtual ~CLASS() = default;

//        virtual void CleanseVertices() = 0;

//        virtual void CleanseEdges() = 0;

//        virtual void CleanseSimplices() = 0;

//        virtual void CleanseAll() = 0;
        
        virtual Int VertexCount() const = 0;

        virtual Int EdgeCount() const = 0;
        
        virtual Int SimplexCount() const = 0;
        
        virtual void Compress() = 0;

        virtual void LoadFromMesh( const MeshBase_T & M ) = 0;
        
        virtual void LoadFromMesh( const MeshBase_T & M, const Tensor2<Real,Int> & u ) = 0;
        
        virtual std::unique_ptr<MeshBase_T> CreateMesh() = 0;

        virtual Tensor2<Real,Int> & VertexData() = 0;
        
//        virtual Int CollapseEdge( const Vertex_T v_0, const Vertex_T v_1 ) = 0;
//
//        virtual Int CollapseEdge( const Edge_T e ) = 0;
//
//        virtual Int SplitEdge( const Vertex_T v_0, const Vertex_T v_1 ) = 0;
//
//        virtual Int SplitEdge( const Edge_T e ) = 0;


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
            return TO_STD_STRING(CLASS)+"<"+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
    };
    
} // namespace Repulsor

#undef CLASS
