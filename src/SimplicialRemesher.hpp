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
 
    template<int DOM_DIM, int AMB_DIM, typename Real_, typename Int_>
    class SimplicialRemesher : public SimplicialRemesherBase<Real_,Int_>
    {
    public:
        
        using Base_T  = SimplicialRemesherBase<Real_,Int_>;
        
        using Real    = typename Base_T::Real;
        using Int     = typename Base_T::Int;
        
        using Vertex_T   = typename Base_T::Vertex_T;
        using Edge_T     = typename Base_T::Edge_T;
        using Simplex_T  = typename Base_T::Simplex_T;
        
//        using MeshBase_T = typename Base_T::MeshBase_T;
        
        using Pair_T             = std::pair<Vertex_T,Vertex_T>;
        
        using EdgeContainer_T    = Tensor2<Vertex_T,Int>;
        using SimplexContainer_T = Tensor2<Vertex_T,Int>;
        
        using VertexList_T       = SortedList<Vertex_T,Int>;
        using SimplexList_T      = SortedList<Simplex_T,Int>;
        
        using Vector_T           = Tensors::Tiny::Vector<AMB_DIM,Real,Int>;

//        using BoolContainer_T = std::vector<bool>;
        using BoolContainer_T = Tensor1<bool,Int>;
        
        
    protected:
        
        static constexpr Int zero = 0;
        static constexpr Int one  = 1;
        static constexpr Int two  = 2;
        
        static constexpr Int S_vertex_count = DOM_DIM+1;
        static constexpr Int S_edge_count   = ((DOM_DIM+1)*DOM_DIM)/2;
        
        static constexpr Int V_max_simplex_valence = (DOM_DIM == 1) ? 2 : (DOM_DIM == 2) ? 9 : 100;
        
        Tensor2<Real,Int> V_coords; // vertex coordinates
        Tensor2<Real,Int> V_data;   // extra data per vertex which we will be transformed
        Tensor1<Int, Int> V_lookup;
        std::vector<SimplexList_T> V_parent_simplices;
        BoolContainer_T V_active;
        BoolContainer_T V_modified;

        EdgeContainer_T edges;
        BoolContainer_T E_active;
        std::vector<SimplexList_T> E_parent_simplices;
        
        SimplexContainer_T simplices;
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
        
        bool compressed = false;
        bool with_data = false;
        
    public:
        
        SimplicialRemesher() = default;
        
//        explicit SimplicialRemesher( const Mesh_T & M )
//        {
//            Init();
//
//            this->LoadMesh(M);
//        }
//
//        explicit SimplicialRemesher( const Mesh_T & M, const Tensor2<Real,Int> & u )
//        {
//            Init();
//
//            this->LoadMesh(M,u);
//        }
        

        template<typename Real_1, typename Int_1>
        SimplicialRemesher(
            ptr<Real_1> V_coords_ ,  const Int vertex_count_,
            ptr< Int_1> simplices_,   const Int simplex_count_,
            const Int thread_count_ = 1
        )
        {
            Real * null = nullptr;
            
            this->LoadMesh(
                V_coords_ ,  vertex_count_,
                simplices_,  simplex_count_,
                null,        Int(0),
                thread_count_
           );
        }
        
        template<typename Real_1, typename Int_1, typename Real_2>
        SimplicialRemesher(
            ptr<Real_1> V_coords_ ,  const Int vertex_count_,
            ptr< Int_1> simplices_,  const Int simplex_count_,
            ptr<Real_2> V_data_,     const Int V_data_dim_,
            const Int thread_count_ = 1
        )
        {
            this->LoadMesh(
                V_coords_ ,  vertex_count_,
                simplices_,  simplex_count_,
                V_data_,     V_data_dim_,
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
            ptr<Real> V_coords_ ,  const Int vertex_count_,
            ptr<Int>  simplices_,  const Int simplex_count_,
            ptr<Real> V_data_,     const Int V_data_dim_,
            const Int thread_count_ = 1
        ) override
        {
            ptic(className()+"::LoadMesh");

            vertex_count  = vertex_count_;
            edge_count    = 0;
            simplex_count = simplex_count_;
            thread_count  = thread_count_;
            
            with_data = ( V_data_ != nullptr ) && ( V_data_dim_ > zero );
            
            ptic("Allocations");
//            max_vertex_count  = vertex_count + S_vertex_count * simplex_count;
//            max_edge_count    = S_edge_count * simplex_count + S_edge_count * simplex_count;
//            max_simplex_count = simplex_count + simplex_count;
            
            max_vertex_count  = vertex_count;
            max_edge_count    = simplex_count * S_edge_count;
            max_simplex_count = simplex_count;
            
            V_parent_simplices = std::vector<SimplexList_T> ( max_vertex_count );
            V_coords           = Tensor2<Real,Int>  ( max_vertex_count, AMB_DIM );
            V_modified         = BoolContainer_T    ( max_vertex_count, false );
            V_active           = BoolContainer_T    ( max_vertex_count );
            std::fill( &V_active[0           ], &V_active[    vertex_count], true  );
            std::fill( &V_active[vertex_count], &V_active[max_vertex_count], false );
            
            if( with_data )
            {
                V_data = Tensor2<Real,Int>  ( max_vertex_count, V_data_dim_ );
            }
            
            E_parent_simplices = std::vector<SimplexList_T> ( max_edge_count );
            edges              = EdgeContainer_T    ( max_edge_count, 2, -1 );
            E_active           = BoolContainer_T    ( max_edge_count, false );
            
            
            simplices          = SimplexContainer_T ( max_simplex_count, DOM_DIM+1, -1 );
            S_active           = BoolContainer_T    ( max_simplex_count );
            std::fill( &S_active[0            ], &S_active[    simplex_count], true  );
            std::fill( &S_active[simplex_count], &S_active[max_simplex_count], false );
            
            ptoc("Allocations");

            ptic("Copy");

            copy_buffer( V_coords_,  V_coords.data(),  vertex_count * AMB_DIM );
            copy_buffer( simplices_, simplices.data(), simplex_count * S_vertex_count );
            
            if( with_data )
            {
               copy_buffer( V_data_,  V_data.data(), vertex_count * V_data_dim_ );
            }
            ptoc("Copy");

            ptic("ComputeConnectivity");
            // Must be run single-threaded
            for( Simplex_T s = 0; s < simplex_count; ++s )
            {
                ComputeSimplexConnectivity(s);
            }
            ptoc("ComputeConnectivity");
            compressed = true;

            ptoc(className()+"::LoadMesh");
            
        } // LoadMesh
        
    public:
        
        virtual void Compress() override
        {
            // Packs all data in V_coords, edges, simplices to the beginning of the arrays.
            
            ptic(className()+"::Compress");
            
            if( compressed )
            {
                ptoc(className()+"::Compress");
                return;
            }
            
            V_lookup.Resize(vertex_count);
            
            Int v_count = 0;
            
            for( Int v = 0; v < vertex_count; ++v )
            {
                if( V_active[v] && (V_parent_simplices[v].Size() > 0) )
                {
                    copy_buffer<AMB_DIM,Sequential>( V_coords.data(v), V_coords.data(v_count) );

                    if( with_data )
                    {
                        copy_buffer<VarSize,Sequential>( V_data.data(v), V_data.data(v_count), V_data.Dimension(1) );
                    }
                    
                    V_active[v]         = false;
                    V_modified[v]       = false;
                    V_active[v_count]   = true;
                    V_modified[v_count] = false;
                    
                    V_lookup[v] = v_count;
                    ++v_count;
                }
            }
            
            Int e_count = 0;
            
            for( Int e = 0; e < edge_count; ++e )
            {
                if( E_active[e] && (E_parent_simplices[e].Size() > 0) )
                {
                    const Vertex_T v_0 = V_lookup[edges(e,0)];
                    const Vertex_T v_1 = V_lookup[edges(e,1)];
                    
                    if( (v_0 < zero) || (v_1 < zero) )
                    {
                        eprint(className()+"::Compress: invalid vertex found in edge "+ToString(e)+".");
                    }
                    
                    if( v_0 <= v_1 )
                    {
                        edges(e_count,0) = v_0;
                        edges(e_count,1) = v_1;
                    }
                    else
                    {
                        edges(e_count,0) = v_1;
                        edges(e_count,1) = v_0;
                    }

                    E_active[e]       = false;
                    E_active[e_count] = true;
                    
                    ++e_count;
                }
            }
            
            Int s_count = 0;
            
            for( Int s = 0; s < simplex_count; ++s )
            {
                if( S_active[s] )
                {
                    for( Int i = 0; i < S_vertex_count; ++i )
                    {
                        const Vertex_T v = V_lookup[simplices(s,i)];
                        
                        if( v < zero )
                        {
                            eprint(className()+"::Compress: invalid vertex "+ToString(v)+" found in simplex "+ToString(s)+".");
                        }
                        
                        simplices(s_count,i) = v;
                    }
                    
                    S_active[s]         = false;
                    S_active[s_count]   = true;
                    
                    ++s_count;
                }
            }
            
            vertex_count  = v_count;
            edge_count    = e_count;
            simplex_count = s_count;
            
            compressed = true;
            
            // TODO: set ?_parent_simplices etc. back onto a valid state.
            
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

#include "SimplicialRemesher/UnifyEdgeLengths.hpp"
  
        
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
                    for( Int i = 0; i < S_vertex_count; ++i )
                    {
                        Vertex_T v = simplices(s,i);
                        
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
                        Edge_T e = SimplexFindEdge(s,j);
                        
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
            
            ptoc(className()+"::SelfCheck");
        }
        
    public:
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        
        
        virtual const Tensor2<Real,Int> & VertexCoordinates() override
        {
            return V_coords;
        }
        
//        virtual Tensor2<Real,Int> VertexData() override
//        {
//            Compress();
//            return Tensor2<Real,Int>( V_data.data(), vertex_count, V_data.Dimension(1) );
//        }
        
        virtual const Tensor2<Real,Int> & VertexData() override
        {
            return V_data;
        }
        
        virtual const Tensor2<Int,Int> & Simplices() override
        {
            return simplices;
        }
        
        virtual Int DomDim() const override
        {
            return DOM_DIM;
        }
        
        virtual Int AmbDim() const override
        {
            return AMB_DIM;
        }
        
    private:
        
        static std::string className()
        {
            return "SimplicialRemesher<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+">";
        }
        
    };
    
        
} // namespace Repulsor