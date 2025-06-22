#pragma once

namespace Repulsor
{
    

    template<typename Real_, typename Int_, typename ExtReal_, typename ExtInt_>
    class SimplicialRemesherBase
    {
        static_assert(FloatQ<Real_>,"");
        static_assert(FloatQ<ExtReal_>,"");
        
        static_assert(IntQ<Int_>,"");
        static_assert(IntQ<ExtInt_>,"");
        
    public:
        
        using Real       = Real_;
        using Int        = Int_;
        using ExtReal    = ExtReal_;
        using ExtInt     = ExtInt_;
        
        using Vertex_T   = Int;
        using Edge_T     = Int;
        using Simplex_T  = Int;

        // Default constructor
        SimplicialRemesherBase() = default;
        // Destructor
        virtual ~SimplicialRemesherBase() = default;
        // Copy constructor
        SimplicialRemesherBase( const SimplicialRemesherBase & other ) = default;
        // Copy assignment operator
        SimplicialRemesherBase & operator=( const SimplicialRemesherBase & other ) = default;
        // Move constructor
        SimplicialRemesherBase( SimplicialRemesherBase && other ) = default;
        // Move assignment operator
        SimplicialRemesherBase & operator=( SimplicialRemesherBase && other ) = default;
        
        virtual Int VertexCount() const = 0;

        virtual Int EdgeCount() const = 0;
        
        virtual Int SimplexCount() const = 0;
        
        virtual void Compress() = 0;

        virtual void LoadMesh(
            cptr<ExtReal> vertex_coords_, const Int vertex_count_,    const bool vertex_coords_ColMajorQ,
            cptr<ExtInt>  simplices_,     const Int simplex_count_,   const bool simplices_ColMajorQ,
            cptr<ExtReal> vertex_data_,   const Int vertex_data_dim_, const bool vertex_data_ColMajorQ,
            const Int thread_count_ = 1
        ) = 0;
        
//        virtual void LoadMesh( cref<MeshBase_T> M ) = 0;
//
//        virtual void LoadMesh( cref<MeshBase_T> M, cref<Tensor2<Real,Int>> V_data_ ) = 0;

//        virtual void LoadMesh( cref<MeshBase_T> M, cptr<Real> V_data_, const Int data_dim ) = 0;
//        
//        virtual void LoadMesh_External( cref<MeshBase_T> M, cptr<ExtReal> V_data_, const Int data_dim ) = 0;
        
//        virtual std::unique_ptr<MeshBase_T> CreateMesh() = 0;

        virtual cref<Tensor2<Real,Int>> VertexCoordinates() const = 0;
        
        virtual cref<Tensor2<Real,Int>> VertexData() const = 0;
        
        virtual cref<Tensor3<Real,Int>> VertexErrorQuadrics() const = 0;
        
        virtual cref<Tensor3<Real,Int>> SimplexErrorQuadrics() const = 0;
        
        virtual cref<Tensor2<Int,Int>> Simplices() const = 0;
        
        virtual Tensor1<Real,Int> SquaredEdgeLengths() = 0;
        
        virtual Int DomDim() const = 0;
        
        virtual Int AmbDim() const = 0;
        
        virtual Int DataDim() const = 0;
        
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
        
        virtual void TangentialSmoothing( const Int max_iter = 1 ) = 0;
        
        virtual void SelfCheck() = 0;
        
    public:

        static std::string className()
        {
            return std::string("SimplicialRemesherBase")+"<"+TypeName<Real>+","+TypeName<Int>+">";
        }
        
        virtual std::string ClassName() const
        {
            return className();
        }
        
    }; // class SimplicialRemesherBase
    
} // namespace Repulsor
