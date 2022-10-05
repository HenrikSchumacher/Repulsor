#pragma once

#define CLASS SimplicialRemesherBase

namespace Repulsion
{
    template<typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS
    {
        
    public:
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

//        virtual void UpdateFromMesh( const MeshBase_T & M ) = 0;
        
        virtual std::unique_ptr<MeshBase_T> CreateMesh() = 0;

        virtual Int CollapseEdge( const Vertex_T v_0, const Vertex_T v_1 ) = 0;
        
        virtual Int CollapseEdge( const Edge_T e ) = 0;
        
        virtual Int SplitEdge( const Vertex_T v_0, const Vertex_T v_1 ) = 0;
        
        virtual Int SplitEdge( const Edge_T e ) = 0;
        
        virtual bool UnifyEdgeLengths( const Real collapse_threshold, const Real split_threshold ) = 0;
        
        virtual void SelfCheck() = 0;
        
    public:
        virtual CLASS & DownCast() = 0;
        
        virtual const CLASS & DownCast() const = 0;
        
        virtual std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
    };
}

#undef CLASS
