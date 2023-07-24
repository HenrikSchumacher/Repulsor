#pragma once

//#define REMESHER_DEBUG
//#define REMESHER_VERBATIM

#include <functional>


namespace Repulsor
{
    // A hash function used to hash a pair of any kind
    struct PairHasher
    {
        template <class T1, class T2>
        Size_T operator()( const std::pair<T1, T2> & p ) const
        {
            // See https://stackoverflow.com/a/23860042/8248900
            
            Size_T hash_1 = std::hash<T1>{}(p.first);
            Size_T hash_2 = std::hash<T2>{}(p.second);
     
            hash_1 ^= hash_2 + 0x9e3779b9 + (hash_1<<6) + (hash_1>>2);
            return hash_1;
        }
    };
 
    template<int DOM_DIM, int AMB_DIM, typename Real_, typename Int_, typename ExtReal_, typename ExtInt_>
    class SimplicialRemesher : public SimplicialRemesherBase<Real_,Int_,ExtReal_,ExtInt_>
    {
    public:
        
        using Base_T     = SimplicialRemesherBase<Real_,Int_,ExtReal_,ExtInt_>;
        
        using Real       = typename Base_T::Real;
        using Int        = typename Base_T::Int;
        using ExtReal    = typename Base_T::ExtReal;
        using ExtInt     = typename Base_T::ExtInt;
        
        using Vertex_T   = typename Base_T::Vertex_T;
        using Edge_T     = typename Base_T::Edge_T;
        using Simplex_T  = typename Base_T::Simplex_T;
        
//        using MeshBase_T = typename Base_T::MeshBase_T;
        
        using Pair_T             = std::pair<Vertex_T,Vertex_T>;
        
        using EdgeContainer_T    = Tensor2<Vertex_T,Int>;
        using SimplexContainer_T = Tensor2<Vertex_T,Int>;
        
        using VertexList_T       = SortedList<Vertex_T,Int>;
        using SimplexList_T      = SortedList<Simplex_T,Int>;
        
        using Vector_T           = Tiny::Vector<AMB_DIM,        Real,Int>;
        using Matrix_T           = Tiny::Matrix<AMB_DIM,AMB_DIM,Real,Int>;

        using Quadric_T          = Tiny::Matrix<AMB_DIM+1,AMB_DIM+1,Real,Int>;
        
//        using BoolContainer_T = std::vector<bool>;
        using BoolContainer_T = Tensor1<bool,Int>;
        
        
    protected:
        
        static constexpr Int zero = 0;
        static constexpr Int one  = 1;
        static constexpr Int two  = 2;
        
        static constexpr Int S_vertex_count = DOM_DIM+1;
        static constexpr Int S_edge_count   = ((DOM_DIM+1)*DOM_DIM)/2;
        
        static constexpr Int V_max_simplex_valence = (DOM_DIM == 1) ? 2 : (DOM_DIM == 2) ? 9 : 100;
        
        const Real sqrt_eps = std::sqrt(Scalar::eps<Real>);
        
        
        Tiny::Vector<S_edge_count,Int,Int>    Tri_i;
        Tiny::Vector<S_edge_count,Int,Int>    Tri_j;
        Tiny::Matrix<AMB_DIM,AMB_DIM,Int,Int> Lin_k;
        
        Tensor2<Real,Int> V_coords; // vertex coordinates
        Tensor2<Real,Int> V_data;   // extra data per vertex which we will be transformed
        Tensor3<Real,Int> V_quadrics; // vertex error quadrics
        Tensor1<Real,Int> V_charges; // vertex charges
        Tensor1<Int, Int> V_lookup;
        std::vector<SimplexList_T> V_parent_simplices;
        BoolContainer_T V_active;
        BoolContainer_T V_modified;

        EdgeContainer_T edges;
        BoolContainer_T E_active;
        std::vector<SimplexList_T> E_parent_simplices;
    
        SimplexContainer_T simplices;
        Tensor3<Real,Int>  S_quadrics; // simplex error quadrics
        
        BoolContainer_T S_active;

        VertexList_T V_0_neighbors;
        VertexList_T V_1_neighbors;
        VertexList_T E_neighbors;
        VertexList_T E_opp_vertices;
        
        std::unordered_map<Pair_T,Int,PairHasher> edge_lookup;
        
        Vertex_T opp_buffer     [DOM_DIM+1];
        Vertex_T simplex_buffer [DOM_DIM+1];
        
        Int     vertex_count;
        Int max_vertex_count;
        
        Int     edge_count;
        Int max_edge_count;
        
        Int     simplex_count;
        Int max_simplex_count;
        
        Int thread_count = 1;
        
        bool compressedQ = false;
        bool with_data = false;
        
    public:
        
        void Init()
        {
            Int k = 0;
            
            for( Int i = 0; i < DOM_DIM+1; ++i )
            {
                for( Int j = i+1; j < DOM_DIM+1; ++j )
                {
                    Tri_i[k] = i;
                    Tri_j[k] = j;
                    Lin_k[i][j] = Lin_k[j][i] = k;
                    ++k;
                }
            }
        }
        
        SimplicialRemesher() = default;
        
//        explicit SimplicialRemesher( const Mesh_T & M )
//        {
//            Init();
//
//            this->LoadMesh(M);
//        }
//
//        explicit SimplicialRemesher( cref<Mesh_T> M, cref<Tensor2<Real,Int>> u )
//        {
//            Init();
//
//            this->LoadMesh(M,u);
//        }
        

        template<typename ExtReal, typename ExtInt>
        SimplicialRemesher(
            cptr<ExtReal> V_coords_ ,
            const Int vertex_count_,
            const bool vertex_coords_ColMajorQ,
            cptr< ExtInt> simplices_,
            const Int simplex_count_,
            const bool simplices_ColMajorQ,
            const Int thread_count_ = 1
        )
        {
            Init();
            
            Real * null = nullptr;
            
            this->LoadMesh(
                V_coords_ ,  vertex_count_,  vertex_coords_ColMajorQ,
                simplices_,  simplex_count_, simplices_ColMajorQ,
                null,        Int(0),         false,
                thread_count_
           );
        }
        
        SimplicialRemesher(
            cptr<ExtReal> vertex_coords_ ,
            const Int vertex_count_,
            const bool    vertex_coords_ColMajorQ,
            cptr<ExtInt>  simplices_,
            const Int     simplex_count_,
            const bool    simplices_ColMajorQ,
            cptr<ExtReal> vertex_data_,
            const Int     vertex_data_dim_,
            const bool    vertex_data_ColMajorQ,
            const Int     thread_count_ = 1
        )
        {
            Init();
            
            this->LoadMesh(
                vertex_coords_ , vertex_count_,    vertex_coords_ColMajorQ,
                simplices_,      simplex_count_,   simplices_ColMajorQ,
                vertex_data_,    vertex_data_dim_, vertex_data_ColMajorQ,
                thread_count_
           );
        }
        
        virtual ~SimplicialRemesher() = default;
        
        
        virtual Int VertexCount() const override
        {
            return vertex_count;
        }

        virtual Int EdgeCount() const override
        {
            return edge_count;
        }
        
        virtual Int SimplexCount() const override
        {
            return simplex_count;
        }
    
    public:
        
        void LoadMesh(
            cptr<ExtReal> vertex_coords_ ,
            const Int     vertex_count_,
            const bool    vertex_coords_ColMajorQ,
            cptr<ExtInt>  simplices_,
            const Int     simplex_count_,
            const bool    simplices_ColMajorQ,
            cptr<ExtReal> vertex_data_,
            const Int     vertex_data_dim_,
            const bool    vertex_data_ColMajorQ,
            const Int     thread_count_ = 1
        ) override
        {
            ptic(className()+"::LoadMesh");

            vertex_count  = vertex_count_;
            edge_count    = 0;
            simplex_count = simplex_count_;
            thread_count  = thread_count_;
            
            with_data = ( vertex_data_ != nullptr ) && ( vertex_data_dim_ > zero );
            
            ptic("Allocations");
//            max_vertex_count  = vertex_count + S_vertex_count * simplex_count;
//            max_edge_count    = S_edge_count * simplex_count + S_edge_count * simplex_count;
//            max_simplex_count = simplex_count + simplex_count;
            
            max_vertex_count  = vertex_count;
            max_edge_count    = simplex_count * S_edge_count;
            max_simplex_count = simplex_count;
            
            V_coords           = Tensor2<Real,Int>  ( max_vertex_count, AMB_DIM );
            V_charges          = Tensor1<Real,Int>  ( max_vertex_count, Scalar::One<Real> );
            V_modified         = BoolContainer_T    ( max_vertex_count, false );
            V_active           = BoolContainer_T    ( max_vertex_count );
            
            if( with_data )
            {
                V_data = Tensor2<Real,Int>  ( max_vertex_count, vertex_data_dim_ );
            }
            
            edges              = EdgeContainer_T    ( max_edge_count, 2, -1 );
            E_active           = BoolContainer_T    ( max_edge_count, false );
            
            simplices          = SimplexContainer_T ( max_simplex_count, DOM_DIM+1, -1 );
            S_active           = BoolContainer_T    ( max_simplex_count );
            
            std::fill( &V_active  [0            ], &V_active  [    vertex_count ], true  );
            std::fill( &V_active  [vertex_count ], &V_active  [max_vertex_count ], false );
            std::fill( &E_active  [0],             &E_active  [max_edge_count   ], false );
            std::fill( &S_active  [0            ], &S_active  [    simplex_count], true  );
            std::fill( &S_active  [simplex_count], &S_active  [max_simplex_count], false );
            
            V_parent_simplices = std::vector<SimplexList_T> ( max_vertex_count );
            E_parent_simplices = std::vector<SimplexList_T> ( max_edge_count );
            
            ptoc("Allocations");

            ptic("Copy");

            if( vertex_coords_ColMajorQ )
            {
                V_coords.ReadTransposed(vertex_coords_);
            }
            else
            {
                V_coords.Read(vertex_coords_);
            }
            
            if( simplices_ColMajorQ )
            {
                simplices.ReadTransposed(simplices_);
            }
            else
            {
                simplices.Read(simplices_);
            }
            
            if( with_data )
            {
                if( vertex_data_ColMajorQ )
                {
                    V_data.ReadTransposed( vertex_data_ );
                }
                else
                {
                    V_data.Read( vertex_data_ );
                }
            }
            
            ptoc("Copy");

            
            ptic("ComputeSimplexConnectivity");
            for( Simplex_T s = 0; s < simplex_count; ++s )
            {
              ComputeSimplexConnectivity(s);
            }
            ptoc("ComputeSimplexConnectivity");
            
            compressedQ = true;
            
            // ComputeErrorQuadrics();
            
            ptoc(className()+"::LoadMesh");
            
        } // LoadMesh
        
        
    public:
        
        virtual void Compress() override
        {
            // Packs all data in V_coords, edges, simplices to the beginning of the arrays.
            
            ptic(className()+"::Compress");
            
            if( compressedQ )
            {
                ptoc(className()+"::Compress");
                return;
            }
            
            const Int old_vertex_count = vertex_count;
            
            V_lookup.Resize(old_vertex_count);
            
            vertex_count = 0;
            
            for( Int v = 0; v < old_vertex_count; ++v )
            {
                if( V_active[v] && (V_parent_simplices[v].Size() > 0) )
                {
                    copy_buffer<AMB_DIM>( V_coords.data(v), V_coords.data(vertex_count) );

                    if( with_data )
                    {
                        copy_buffer<VarSize,Sequential>( V_data.data(v), V_data.data(vertex_count), V_data.Dimension(1) );
                    }
                    
//                    copy_buffer<(AMB_DIM+1)*(AMB_DIM+1)>( V_quadrics.data(v), V_quadrics.data(v_count) );
                    
                    V_charges[v]        = Scalar::Zero<Real>;
                    V_active[v]         = false;
                    V_modified[v]       = false;
                    
                    V_charges[vertex_count]  = Scalar::One<Real>;
                    V_active[vertex_count]   = true;
                    V_modified[vertex_count] = false;
                    
                    V_lookup[v] = vertex_count;
                    
                    V_parent_simplices[vertex_count].Clear();
                    
                    ++vertex_count;
                }
            }
            
//            // Should be superfluous.
//            for( Int v = old_vertex_count; v < max_vertex_count; ++v )
//            {
//                V_parent_simplices[v].Clear();
//            }
            
            
            std::fill( &E_active  [0],             &E_active  [max_edge_count   ], false );
            
            for( Int e = 0; e < max_edge_count; ++e )
            {
                E_parent_simplices[e].Clear();
            }
            
            const Int old_edge_count = edge_count;
            
            for( Int e = 0; e < old_edge_count; ++e )
            {
                E_parent_simplices[e].Clear();
            }
            
//            // Should be superfluous.
//            for( Int e = old_edge_count; e < max_edge_count; ++e )
//            {
//                E_parent_simplices[e].Clear();
//            }
            
            edge_lookup = std::unordered_map<Pair_T,Int,PairHasher>();
            
            edge_count = 0;
            
//
//            for( Int e = 0; e < old_edge_count; ++e )
//            {
//                if( E_active[e] && (E_parent_simplices[e].Size() > 0) )
//                {
//                    const Vertex_T v_0 = V_lookup[edges[e][0]];
//                    const Vertex_T v_1 = V_lookup[edges[e][1]];
//
//                    if( (v_0 < zero) || (v_1 < zero) )
//                    {
//                        eprint(className()+"::Compress: invalid vertex found in edge "+ToString(e)+".");
//                    }
//
//                    if( v_0 == v_1 )
//                    {
//                        eprint(className()+"::Compress: edge "+ToString(e)+" is degenerate.");
//                    }
//
//                    CreateEdge( v_0, v_1 );
////
////                    edges[edge_count][0] = std::min(v_0,v_1);
////                    edges[edge_count][1] = std::max(v_0,v_1);
////
////                    E_active[e]       = false;
////                    E_active[edge_count] = true;
//
////                    ++edge_count;
//                }
//            }
            
            const Int old_simplex_count = simplex_count;
            
            simplex_count = 0;
            
            for( Int s = 0; s < old_simplex_count; ++s )
            {
                if( S_active[s] )
                {
                    for( Int i = 0; i < S_vertex_count; ++i )
                    {
                        const Vertex_T v = V_lookup[simplices[s][i]];
                        
                        if( v < zero )
                        {
                            eprint(className()+"::Compress: invalid vertex "+ToString(v)+" found in simplex "+ToString(s)+".");
                        }
                        
                        simplices[simplex_count][i] = v;
                    }
                    
                    S_active[s]              = false;
                    S_active[simplex_count]  = true;
                    
                    ComputeSimplexConnectivity(simplex_count);
                    
                    ++simplex_count;
                }
            }
            
            compressedQ = true;

            ptoc(className()+"::Compress");
        } // Compress
        
//        virtual std::unique_ptr<MeshBase_T> CreateMesh() override
//        {
//            Compress();
//            return std::make_unique<Mesh_T>(
//                V_coords.data(), vertex_count,
//                simplices.data(), simplex_count,
//                thread_count
//            );
//        } //CreateMesh

#include "SimplicialRemesher/Vertices.hpp"
#include "SimplicialRemesher/Edges.hpp"
#include "SimplicialRemesher/Simplices.hpp"
        
#include "SimplicialRemesher/Checks.hpp"
#include "SimplicialRemesher/SplitEdges.hpp"
#include "SimplicialRemesher/CollapseEdges.hpp"
#include "SimplicialRemesher/FlipEdges.hpp"

//#include "SimplicialRemesher/ErrorQuadrics.hpp"
        
#include "SimplicialRemesher/UnifyEdgeLengths.hpp"
#include "SimplicialRemesher/TangentialSmoothing.hpp"
        
//##############################################################################################
//      SelfCheck
//##############################################################################################
            
        
        virtual void SelfCheck() override
        {
            // Check whether all vertices of each undeleted simplex are undeleted.
            ptic(className()+"::SelfCheck");
            
            
            Int deleted_parent_simplex_found_in_vertex = 0;
            
            Int deleted_vertices_found_in_edges = 0;
            Int deleted_parent_simplices_found_in_edges = 0;
            Int edges_agnostic_of_their_parent_simplices = 0;
            Int simplices_agnostic_of_their_child_edges = 0;
            
            Int deleted_vertices_found_in_simplices = 0;
            
            
            Int simplices_agnostic_of_their_child_vertices = 0;
            Int vertices_agnostic_of_their_parent_simplices = 0;
            
            Int invalid_edges_found_in_simplices = 0;
            Int deleted_edges_found_in_simplices = 0;
            
            Int simplices_with_duplicates = 0;
                            
            
            print("Checking vertices...");
            for( Int v = 0; v < vertex_count; ++v )
            {
                if( V_active[v] )
                {
                    if( V_parent_simplices[v].Size() > V_max_simplex_valence )
                    {
                        wprint(className()+"::SelfCheck: Vertex "+ToString(v)+" has simplex valence "+Tools::ToString(V_parent_simplices[v].Size())+". Allowed maximum is "+Tools::ToString(V_max_simplex_valence)+".");
                    }
                    
                    for( Simplex_T s : V_parent_simplices[v] )
                    {
                        if( !S_active[s] )
                        {
                            if( deleted_parent_simplex_found_in_vertex==0 )
                            {
                                eprint(className()+"::SelfCheck: Vertex "+ToString(v)+" has deleted parent simplex "+ToString(s)+".");
                            }
                            ++deleted_parent_simplex_found_in_vertex;
                        }
                        else
                        {
                            bool contained = false;
                            for( Int i = 0; i < S_vertex_count; ++i )
                            {
                                Vertex_T w = simplices(s,i);
                                
                                if( v == w )
                                {
                                    contained = true;
                                }
                            }
                            
                            if( !contained )
                            {
                                if( simplices_agnostic_of_their_child_vertices==0 )
                                {
                                    eprint(className()+"::SelfCheck: simplex "+ToString(s)+" does not know about its child vertex "+ToString(v)+".");
                                    
                                    print(
                                          "V_parent_simplices[v] = " + V_parent_simplices[v].ToString()
                                    );
                                }
                                ++simplices_agnostic_of_their_child_vertices;
                            }
                        }
                    }
                }
                
            } // for( Int v = 0; v < vertex_count; ++v )
            print("Done with checking vertices.");
            
            
            print("Checking edges...");
            for( Edge_T e = 0; e < edge_count; ++e )
            {
                if( E_active[e] )
                {
                    if( edges(e,0) == edges(e,1) )
                    {
                        eprint(className()+"::SelfCheck: Edge "+ToString(e)+" = { "+ToString(edges(e,0))+","+ToString(edges(e,0))+"} is degenerate.");
                    }
                    
                    for( Int i = 0; i < 2; ++i )
                    {
                        Vertex_T v = edges(e,i);
                        
                        if( !V_active[v] )
                        {
                            if( deleted_vertices_found_in_edges==0 )
                            {
                                eprint(className()+"::SelfCheck: Edge "+ToString(e)+" contains deleted vertex "+ToString(v)+".");
                            }
                            ++deleted_vertices_found_in_edges;
                        }
                    }
                    
                    for( Simplex_T s : E_parent_simplices[e] )
                    {
                        if( !S_active[s] )
                        {
                            if( deleted_parent_simplices_found_in_edges==0 )
                            {
                                eprint(className()+"::SelfCheck: Edge "+ToString(e)+" has deleted parent simplex "+ToString(s)+".");
                            }
                            ++deleted_parent_simplices_found_in_edges;
                        }
                        else
                        {
                            bool contained = false;
                            for( Int k = 0; k < S_edge_count; ++k )
                            {
                                Edge_T f = SimplexFindEdge(s, k);
                                if( e == f )
                                {
                                    contained = true;
                                }
                            }
                            
                            if( !contained )
                            {
                                if( simplices_agnostic_of_their_child_edges==0 )
                                {
                                    eprint(className()+"::SelfCheck: Simplex "+ToString(s)+" cannot find its child edge "+ToString(e)+" = { "+ToString(edges(e,0))+", "+ToString(edges(e,0))+" }.");
                                }
                                ++simplices_agnostic_of_their_child_edges;
                            }
                            
                        }
                    }
                }
                
            } // for( Edge_T e = 0; e < edge_count; ++e )
            print("Done with checking edges.");
            
            print("Checking simplices...");
            for( Simplex_T s = 0; s < simplex_count; ++s )
            {
                if( S_active[s] )
                {
                    {
                        bool duplicate_freeQ = true;
                        
                        for( Int i = 0; i < S_vertex_count; ++i )
                        {
                            for( Int j = i+1; j < S_vertex_count; ++j )
                            {
                                duplicate_freeQ = duplicate_freeQ && (simplices[s][i] != simplices[s][j] );
                            }
                        }
                        
                        if( !duplicate_freeQ )
                        {
                            if( simplices_with_duplicates==0 )
                            {
                                eprint(className()+"::SelfCheck: simplex "+ToString(s)+" containes duplicate vertices.");
                            }
                            ++simplices_with_duplicates;
                        }
                    }
                    
                    for( Int i = 0; i < S_vertex_count; ++i )
                    {
                        const Vertex_T v = simplices(s,i);
                        
                        if( !V_active[v] )
                        {
                            if( deleted_vertices_found_in_simplices==0 )
                            {
                                eprint(className()+"::SelfCheck: Deleted vertex "+ToString(v)+" found in simplex "+ToString(s)+".");
                            }
                            ++deleted_vertices_found_in_simplices;
                        }
                        else
                        {
                            bool contained = (V_parent_simplices[v].Find(s) >= 0);
//
                            if( !contained )
                            {
                                if( vertices_agnostic_of_their_parent_simplices==0 )
                                {
                                    eprint(className()+"::SelfCheck: vertex "+ToString(v)+" does not know about its parent simplex "+ToString(s));
                                }
                                ++vertices_agnostic_of_their_parent_simplices;
                            }
                        }
                    }
                    
                    
                    for( Int j = 0; j < S_vertex_count; ++j )
                    {
                        const Edge_T e = SimplexFindEdge(s,j);
                        
                        if( e < 0 )
                        {
                            if( invalid_edges_found_in_simplices==0 )
                            {
                                eprint(className()+"::SelfCheck: Invalid edge "+ToString(e)+" found in simplex "+ToString(s)+".");
                            }
                            ++invalid_edges_found_in_simplices;
                        }
                        else if( !E_active[e] )
                        {
                            if( deleted_edges_found_in_simplices==0 )
                            {
                                eprint(className()+"::SelfCheck: Deleted edge "+ToString(e)+" found in simplex "+ToString(s)+".");
                            }
                            ++deleted_edges_found_in_simplices;
                    
                        }
                        else
                        {
                            bool contained = (E_parent_simplices[e].Find(s) >= 0);
                            
                            if( !contained )
                            {
                                if( edges_agnostic_of_their_parent_simplices==0 )
                                {
                                    eprint(className()+"::SelfCheck: edge "+ToString(e)+" does not know about its parent simplex "+ToString(s)+".");
                                }
                                ++edges_agnostic_of_their_parent_simplices;
                            }
                        }
                    }
                    

                }
            }
            print("Done with checking simplices.");
            
            print("SelfCheck finished.");
            
//            dump(invalid_edges_found_in_simplices);
//            dump(deleted_edges_found_in_simplices);
//            dump(edges_agnostic_of_their_parent_simplices);
//
//            dump(deleted_parent_simplex_found_in_vertex);
//
//            dump(deleted_vertices_found_in_edges);
//            dump(deleted_parent_simplices_found_in_edges);
//            dump(edges_agnostic_of_their_parent_simplices);
//            dump(simplices_agnostic_of_their_child_edges);
//
//            dump(deleted_vertices_found_in_simplices);
//
//            dump(simplices_agnostic_of_their_child_vertices);
//            dump(vertices_agnostic_of_their_parent_simplices);
//
//            dump(invalid_edges_found_in_simplices);
//            dump(deleted_edges_found_in_simplices);
//            dump(simplices_with_duplicates);
            
            ptoc(className()+"::SelfCheck");
        }
        
    public:
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        

        virtual cref<Tensor2<Real,Int>> VertexCoordinates() const override
        {
            return V_coords;
        }

        virtual cref<Tensor2<Real,Int>> VertexData() const override
        {
            return V_data;
        }
        
        virtual cref<Tensor3<Real,Int>> VertexErrorQuadrics() const override
        {
            return V_quadrics;
        }
        
        virtual cref<Tensor3<Real,Int>> SimplexErrorQuadrics() const override
        {
            return S_quadrics;
        }
        
        virtual cref<Tensor2<Int,Int>> Simplices() const override
        {
            return simplices;
        }
        
        
        virtual Tensor1<Real,Int> SquaredEdgeLengths() override
        {
            Compress();
            
            Tensor1<Real,Int> result ( EdgeCount() );
            
            
            ParallelDo(
                [&,this]( const Edge_T e )
                {
                    result[e] = SquaredEdgeLength(e);
                },
                edge_count, thread_count
            );

            return result;
        }
        
        virtual Int DomDim() const override
        {
            return DOM_DIM;
        }
        
        virtual Int AmbDim() const override
        {
            return AMB_DIM;
        }
        
        virtual Int DataDim() const override
        {
            return V_data.Dimension(1);
        }
        
    private:
        
        static std::string className()
        {
            return "SimplicialRemesher<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+">";
        }
        
    };
    
        
} // namespace Repulsor
