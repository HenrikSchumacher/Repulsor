#pragma once

#define CLASS SimplicialRemesher
#define BASE SimplicialRemesherBase<Real,Int,SReal,ExtReal>

//#define REMESHER_DEBUG
//#define REMESHER_VERBATIM

#include <functional>


#define I(x) static_cast<Int>(x)
#define R(x) static_cast<Real>(x)

namespace Repulsor
{
    
//    template <class T>
//    struct PointerPairHasher {
//        std::size_t operator()( const std::pair<T *, T *> & p) const
//        {
//            // See https://stackoverflow.com/a/23860042/8248900
//            std::size_t hash_1 = reinterpret_cast<std::size_t>( p.first );
//            std::size_t hash_2 = reinterpret_cast<std::size_t>( p.second );
//            hash_1 ^= hash_2 + 0x9e3779b9 + (hash_1<<6) + (hash_1>>2);
//            return hash_1;
//        }
//    };

    
    // A hash function used to hash a pair of any kind
    struct PairHasher {
        template <class T1, class T2>
        size_t operator()( const std::pair<T1, T2> & p ) const
        {
            size_t hash_1 = std::hash<T1>{}(p.first);
            size_t hash_2 = std::hash<T2>{}(p.second);
     
            hash_1 ^= hash_2 + 0x9e3779b9 + (hash_1<<6) + (hash_1>>2);
            return hash_1;
        }
    };
 
    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
        using Vertex_T   = Int;
        using Edge_T     = Int;
        using Simplex_T  = Int;
        
        using Pair_T     = std::pair<Vertex_T,Vertex_T>;
        
        using EdgeContainer_T    = Tensor2<Vertex_T,Int>;
        using SimplexContainer_T = Tensor2<Vertex_T,Int>;
        
        using MeshBase_T = SimplicialMeshBase<                Real,Int,SReal,ExtReal>;
        using Mesh_T     = SimplicialMesh    <DOM_DIM,AMB_DIM,Real,Int,SReal,ExtReal>;

//        using BoolContainer_T = std::vector<bool>;
        using BoolContainer_T = Tensor1<bool,Int>;
        
    protected:
        
        static constexpr Int S_vertex_count = DOM_DIM+1;
        static constexpr Int S_edge_count   = ((DOM_DIM+1)*DOM_DIM)/2;
        
        static constexpr Int V_max_simplex_valence = (DOM_DIM == 1) ? 2 : (DOM_DIM == 2) ? 9 : 100;
        
        Int tri_i [S_edge_count]     = {};
        Int tri_j [S_edge_count]     = {};
        Int lin_k [AMB_DIM][AMB_DIM] = {};
        
        Tensor2<Real,Int> V_coords;
        Tensor1<Int,Int>  V_lookup;
        
        EdgeContainer_T edges;
        SimplexContainer_T simplices;
        
        std::vector<SortedList<Simplex_T,Int>> V_parent_simplices;
        std::vector<SortedList<Simplex_T,Int>> E_parent_simplices;
                
        BoolContainer_T V_active;
        BoolContainer_T V_modified;

        BoolContainer_T E_active;

        BoolContainer_T S_active;

        SortedList<Vertex_T,Int> V_0_neighbors;
        SortedList<Vertex_T,Int> V_1_neighbors;
        SortedList<Vertex_T,Int> E_neighbors;
        SortedList<Vertex_T,Int> E_opp_vertices;
        
        std::unordered_map<Pair_T, Int, PairHasher> edge_lookup;

        std::array<Vertex_T,DOM_DIM+1> opp_buffer;
        std::array<Vertex_T,DOM_DIM+1> simplex_buffer;
        
        Int     vertex_count;
        Int max_vertex_count;
        
        Int     edge_count;
        Int max_edge_count;
        
        Int     simplex_count;
        Int max_simplex_count;
        
        Int thread_count = 1;
        
        bool compressed = false;
        
    public:
        
        CLASS() = default;
        
        explicit CLASS( const Mesh_T & M )
        {
            Init();
            
            updateFromMesh(M);
        }
        
        
        virtual ~SimplicialRemesher() = default;
        
        //DONE.
        void Init()
        {
            Int k = 0;
            
            for( Int i = 0; i < DOM_DIM+1; ++i )
            {
                for( Int j = i+1; j < DOM_DIM+1; ++j )
                {
                    tri_i[k] = i;
                    tri_j[k] = j;
                    lin_k[i][j] = lin_k[j][i] = k;
                    ++k;
                }
            }
        }
        
        //DONE.
        void CheckInit() const
        {
            print(ClassName()+"::CheckInit");
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                for( Int j = i+1; j < AMB_DIM; ++j )
                {
                    print("{ "+ToString(i)+"," +ToString(j)+" } -> " + ToString(lin_k[i][j]) );
                }
            }

            for( Int k = 0; k < S_vertex_count; ++k )
            {
                print( ToString(k) +" -> { "+ToString(tri_i[k])+"," +ToString(tri_j[k])+" }");
            }
        }
        
        
        //DONE.
        virtual Int VertexCount() const override
        {
            return vertex_count;
        }

        //DONE.
        virtual Int EdgeCount() const override
        {
            return edge_count;
        }
        
        //DONE.
        virtual Int SimplexCount() const override
        {
            return simplex_count;
        }
        
//        virtual void UpdateFromMesh( const MeshBase_T & M ) override
//        {
//            ptic(className()+"::UpdateFromMesh");
//
//            updateFromMesh( M );
//
//            ptoc(className()+"::UpdateFromMesh");
//        }
    
    private:
        
        //TODO: Check this!
        void updateFromMesh( const MeshBase_T & M )
        {
            ptic(className()+"::updateFromMesh");
            vertex_count  = M.VertexCount();
            edge_count    = 0;
            simplex_count = M.SimplexCount();
            thread_count  = M.ThreadCount();
            
            ptic("Allocations");
            max_vertex_count  = vertex_count + S_vertex_count * simplex_count;
            max_edge_count    = S_edge_count * simplex_count + S_edge_count * simplex_count;
            max_simplex_count = simplex_count + simplex_count;
            
            V_coords   = Tensor2<Real,Int>  ( max_vertex_count, AMB_DIM, 0 );
            V_active   = BoolContainer_T    ( max_vertex_count, false );
            V_modified = BoolContainer_T    ( max_vertex_count, false );
            V_parent_simplices = std::vector<SortedList<Simplex_T, Int>> ( max_vertex_count );
            
            edges      = EdgeContainer_T    ( max_edge_count, 2, -1 );
            E_active   = BoolContainer_T    ( max_edge_count, false );
            E_parent_simplices = std::vector<SortedList<Simplex_T, Int>> ( max_edge_count );
            
            simplices  = SimplexContainer_T ( max_simplex_count, DOM_DIM+1, -1 );
            S_active   = BoolContainer_T    ( max_simplex_count, false );
            
            std::fill( &V_active[0], &V_active[vertex_count],  true );
            std::fill( &S_active[0], &S_active[simplex_count], true );
            ptoc("Allocations");
            
            ptic("Copy");
            M.VertexCoordinates().Write( V_coords.data() );
            M.Simplices().Write( simplices.data() );
            ptoc("Copy");
            
            ptic("ComputeConnectivity");
            // Must be run single-threaded
            for( Simplex_T s = 0; s < simplex_count; ++s )
            {
                Simplex_ComputeConnectivity(s);
            }
            ptoc("ComputeConnectivity");
            compressed = true;
            
            ptoc(className()+"::updateFromMesh");
            
        } // updateFromMesh
        
    public:
        
        virtual void Compress() override
        {
            // Packs all data in V_coords, edges, simplices to the beginning of the arrays.
            
//            SelfCheck();
            
            ptic(className()+"::Compress");
            
            if( compressed )
            {
                ptoc(className()+"::Compress");
                return;
            }
            
            V_lookup.Resize(vertex_count);
            
            V_lookup.Fill(-3);
            
            Int v_count = 0;
            
            for( Int v = 0; v < vertex_count; ++v )
            {
                if( V_active[v] && (V_parent_simplices[v].Size() > 0) )
                {
                   for( Int k = 0; k < AMB_DIM; ++k )
                    {
                        V_coords(v_count,k) = V_coords(v,k);
                    }
                    
                    V_active[v]         = false;
                    V_modified[v]       = false;
                    V_active[v_count]   = true;
                    V_modified[v_count] = false;
                    
                    V_lookup[v] = v_count;
                    v_count++;
                }
            }
            
            
//            V_coords.Resize(v_count,AMB_DIM);
            
            Int e_count = 0;
            
            for( Int e = 0; e < edge_count; ++e )
            {
                if( E_active[e] && (E_parent_simplices[e].Size() > 0) )
                {
                    const Vertex_T v_0 = V_lookup[edges(e,0)];
                    const Vertex_T v_1 = V_lookup[edges(e,1)];
                    
                    if( v_0 < 0  || v_1 < 0)
                    {
                        eprint("Bah!");
                    }
                    
                    if( !V_active[v_0] || !V_active[v_1] )
                    {
                        eprint("Buh!");
                    }
                    
                    edges(e_count,0) = v_0;
                    edges(e_count,1) = v_1;
                    
                    E_active[e]         = false;
                    E_active[e_count]   = true;
                    
                    e_count++;
                }
            }
            
            
            Int s_count = 0;
            
            for( Int s = 0; s < simplex_count; ++s )
            {
                if( S_active[s] )
                {
                    for( Int i = 0; i < S_vertex_count; ++ i )
                    {
                        const Vertex_T v = V_lookup[simplices(s,i)];
                        
                        if( v < 0 )
                        {
                            eprint("Bah!!!");
                            
                            valprint("v",v);
                            valprint("simplices(s,i)",simplices(s,i));
                            
                        }
                        
//                        if( !V_active[v] )
//                        {
//                            eprint("Buh!!!");
//                        }
                        
                        simplices(s_count,i) = v;
                    }
                    
                    S_active[s]         = false;

                    S_active[s_count]   = true;
                    
                    s_count++;
                }
            }
            
            vertex_count  = v_count;
            edge_count    = e_count;
            simplex_count = s_count;
            
            compressed = true;
            
            // TODO: set ?_parent_simplices etc. back onto a valid state.
            
//            SelfCheck();
            
            ptoc(className()+"::Compress");
        } // Compress

        
        virtual std::unique_ptr<MeshBase_T> CreateMesh() override
        {
            Compress();
            return std::make_unique<Mesh_T>(
                V_coords.data(), vertex_count,
                simplices.data(), simplex_count,
                thread_count
            );
        } //CreateMesh
        
//################################################################################################
//##    vertex related
//################################################################################################


        Vertex_T CreateVertex()
        {
            Vertex_T w = vertex_count++;
            
            if( vertex_count >= max_vertex_count)
            {
#ifdef REMESHER_VERBATIM
                print(className()+"::CreateVertex: Reassembling vertex array.");
#endif
                max_vertex_count *= I(2);
                
                V_coords.Resize( max_vertex_count, AMB_DIM );
                V_active.Resize( max_vertex_count );
                V_modified.Resize( max_vertex_count );
                
            }
            
            V_active[w] = true;
            V_modified[w] = true;
            
            vertex_count++;
            
            return w;
        }
        
        void Vertex_Deactivate( const Vertex_T v )
        {
#ifdef REMESHER_VERBATIM
            print(className()+"::Vertex_Deactivate("+ToString(v)+")");
#endif
            V_active[v] = false;
        }
        
        void Vertex_MarkAsModified( const Vertex_T v )
        {
#ifdef REMESHER_VERBATIM
            print(className()+"::Vertex_MarkAsModified("+ToString(v)+")");
#endif
            V_modified[v] = true;
        }
        
        void ComputeVertexPosition( const Vertex_T v_0, const Vertex_T v_1, const Vertex_T w )
        {
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                V_coords(w,k) = 0.5 * ( V_coords(v_0,k) + V_coords(v_1,k) );
            }
        }
        
        void Vertex_CollectNeighborVertices( const Vertex_T v, SortedList<Vertex_T,Int> & neighbors ) const
        {
            neighbors.Clear();
            
            for( Simplex_T s : V_parent_simplices[v] )
            {
                for( Int i = 0; i < S_vertex_count; ++ i )
                {
                    const Vertex_T w = simplices(s,i);
                    
                    if( v != w )
                    {
                        neighbors.Insert(w);
                    }
                }
            }
        }
        
//#####################################################################################################
//##    edge related
//#####################################################################################################
        
        Edge_T CreateEdge( const Vertex_T v_0, const Vertex_T v_1 )
        {
            Pair_T p = std::minmax(v_0,v_1);
            
            Edge_T e = edge_count++;
            
            if( e >= max_edge_count )
            {
#ifdef REMESHER_VERBATIM
                print(className()+"::CreateEdge: Reassembling edge array.");
#endif
                max_edge_count *= I(2);
                
                edges.Resize( max_edge_count, 2 );
                E_active.Resize( max_edge_count );
            }
        
            edge_lookup.insert( {p,e} );
            
            edges(e,0) = p.first;
            edges(e,1) = p.second;
            
            E_active[e] = true;
            
#ifdef REMESHER_DEBUG
            print(className()+"::CreateEdge: Created edge with ID = "+ToString(e)+" and vertices {"+ToString(v_0)+","+ToString(v_1)+"}.");
#endif
            return e;
        }
        
        //DONE.
        Edge_T RequireEdge( const Vertex_T v_0, const Vertex_T v_1 )
        {
            Edge_T e = FindEdge(v_0,v_1);
            
            if( e < I(0) )
            {
                e = CreateEdge(v_0,v_1);
                
#ifdef REMESHER_DEBUG
                print(className()+"RequireEdge: Created edge "+ToString(e)+" with vertices {"+ToString(v_0)+","+ToString(v_1)+"}.");
#endif
            }
            return e;
        }
        
        void DeleteEdge( const Edge_T e )
        {
#ifdef REMESHER_VERBATIM
            print(className()+"::Edge_Delete("+ToString(e)+")");
#endif
            E_active[e] = false;
            LookupErase(e);
        }
        
        //DONE.
        Edge_T FindEdge( const Vertex_T v_0, const Vertex_T v_1 )
        {
            const Pair_T p = std::minmax(v_0,v_1);
            
            return ( edge_lookup.count(p) > 0) ? edge_lookup[p] : -9;
        }
        
        //DONE.
        void LookupErase( const Edge_T e )
        {
            const Pair_T p = std::minmax( edges(e,0), edges(e,1) );
            edge_lookup.erase( p );
        }
        
        //DONE.
        void LookupInsert( const Edge_T e )
        {
            const Pair_T p = std::minmax( edges(e,0), edges(e,1) );
            edge_lookup.insert( {p,e} );
        }
        
        
        //DONE.
        void Edge_CollectOpposingVertices( const Edge_T e, SortedList<Vertex_T, Int> & opp_vertices ) const
        {
            opp_vertices.Clear();
            
            const Vertex_T v_0 = edges(e,0);
            const Vertex_T v_1 = edges(e,1);
            
            // Going through the simplices to find opposing vertices.
            for( Simplex_T s : E_parent_simplices[e] )
            {
                for( Int i = 0; i < S_vertex_count; ++i )
                {
                    const Vertex_T v = simplices(s,i);
                    
                    if( (v != v_0) && (v != v_1) )
                    {
                        opp_vertices.Insert(v);
                    }
                }
            }
        }
        
        //DONE.
        Real SquaredEdgeLength( Edge_T e ) const
        {
            Real L2 = 0.;
            
            Real const * restrict const V = V_coords.data(edges(e,0));
            Real const * restrict const W = V_coords.data(edges(e,1));
            
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                Real delta = V[k] -  W[k];
                L2 += delta * delta;
            }
            
            return L2;
        }
        
//######################################################################################################
//##    simplex related
//######################################################################################################

        //DONE.
        void Simplex_ComputeConnectivity( const Simplex_T s )
        {
//            S_active[s] = true;
            
            for( Int i = 0; i < S_vertex_count; ++i )
            {
                Vertex_T v = simplices(s,i);
                
                V_parent_simplices[v].Insert(s);
                
                for( Int j = i+1; j < S_vertex_count; ++j )
                {
                    Vertex_T w = simplices(s,j);
                    
                    Edge_T e = RequireEdge(v,w);
       
                    if( e < 0 )
                    {
                        print("Uuuuh!");
                    }
                    
                    E_parent_simplices[e].Insert(s);
                }
            }
        }
        
//        Simplex_T * CreateSimplex()
//        {
//            Simplex_T * S = new Simplex_T( SimplexCount(), simplex_buffer );
//            simplices.push_back(S);
//
//            for( Int i = 0; i < Simplex_T::VertexCount(); ++i )
//            {
//                Vertex_T * V = S->Vertex(i);
//                V->InsertParentSimplex(S);
//
//                for( Int j = i+1; j < Simplex_T::VertexCount(); ++j )
//                {
//                    Vertex_T * W = S->Vertex(j);
//
//                    Edge_T * E = RequireEdge(V, W);
//
//                    E->InsertParentSimplex(S);
//                }
//            }
//
//            return S;
//        }
        
        Simplex_T CreateSimplex( const Vertex_T * vertex_list )
        {
            Simplex_T s = simplex_count++;
            
            if( s >= max_simplex_count )
            {
#ifdef REMESHER_VERBATIM
                print(className()+"::CreateSimplex: Reassembling simplex array.");
#endif
                max_simplex_count *= 2;
                
                simplices.Resize(max_simplex_count,DOM_DIM+1);
                S_active.Resize(max_simplex_count);
            }
            
            S_active[s] = true;
            
            copy_buffer(vertex_list, simplices.data(s), S_vertex_count);
            
            Simplex_ComputeConnectivity(s);
            
            return s;
        }
        
        Edge_T Simplex_FindEdge( const Simplex_T s, const Int k )
        {
            const Pair_T p = std::minmax( simplices(s,tri_i[k]), simplices(s,tri_j[k]) );
            
            return ( edge_lookup.count(p) > 0) ? edge_lookup[p] : -13;
        }
        
        void Simplex_OppositeVertices(
            const Simplex_T s,
            const Vertex_T v_0,
            const Vertex_T v_1,
            Vertex_T * restrict opp_vertex_list
        ) const
        {
            Int counter = 0;

            for( Int i = 0; i < S_vertex_count; ++i )
            {
                const Vertex_T v = simplices(s,i);
                
                if( v != v_0 && v!=v_1 )
                {
                    opp_vertex_list[counter] = v;
                    ++counter;
                }
            }
        }
        
        void DeleteSimplex( const Simplex_T s )
        {
            S_active[s] = false;
            
            for( Vertex_T v = 0; v < S_vertex_count; ++v )
            {
                if( v>= 0 )
                {
                    V_parent_simplices[v].Drop(s);
                }
            }
        }

#include "SplitEdge.hpp"
#include "CollapseEdge.hpp"
        
        
            
//######################################################################################################
//##    Unify edge lengths
//######################################################################################################
        
        virtual bool UnifyEdgeLengths(
            const Real collapse_threshold,
            const Real split_threshold
        ) override
        {
            ptic(className()+"::UnifyEdgeLengths");
            
//            SelfCheck();
            
            if( split_threshold < collapse_threshold )
            {
                eprint(className()+"::UnifyEdgeLengths: split_threshold < collapse_threshold. Aborting");
                ptoc(className()+"::UnifyEdgeLengths");
                return 0;
            }
            
            const Real split    = split_threshold * split_threshold;
            const Real collapse = collapse_threshold * collapse_threshold;
            
            Int split_counter = 0;
            Int collapse_counter = 0;
            
            for( Int e = 0; e < edge_count; ++e )
            {
                if( !E_active[e] )
                {
#ifdef REMESHER_VERBATIM
                    wprint(className()+"::UnifyEdgeLengths: Skipping edge "+ToString(e)+" because it is inactive.");
#endif
                    continue;
                }
                
                const Real L2 = SquaredEdgeLength(e);
                
                if( L2 > split )
                {
                    const Int r = SplitEdge(e);
                    if( r >= 0 )
                    {
                        ++split_counter;
                    }
                    else
                    {
#ifdef REMESHER_VERBATIM
                        wprint(className()+"::UnifyEdgeLengths: SplitEdge failed to split edge "+ToString(e)+".");
#endif
                    }
                }
                else
                    if( L2 < collapse )
                {
                    const Int r = CollapseEdge(e);
                    if( r >= 0 )
                    {
                        ++collapse_counter;
                    }
                    else
                    {
#ifdef REMESHER_VERBATIM
                        wprint(className()+"::UnifyEdgeLengths: CollapseEdge failed to collapse edge "+ToString(e)+".");
#endif
                    }
                }
                
            } // for( Int e = 0; e < edge_count; ++e )
            
            valprint("split_counter",split_counter);
            valprint("collapse_counter",collapse_counter);
            
            ptoc(className()+"::UnifyEdgeLengths");
            
//            SelfCheck();
            
            return (split_counter>0) || (collapse_counter>0);
            
        } // UnifyEdgeLengths
  
        
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
                            
            
            print("Cycling over vertices...");
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
                            deleted_parent_simplex_found_in_vertex++;
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
//                                    S->PrintStats();
//                                    V->PrintStats();
//                                    V->PrintParentSimplices();
                                }
                                simplices_agnostic_of_their_child_vertices++;
                            }
                        }
                    }
                }
                
            } // for( Int v = 0; v < vertex_count; ++v )
            
            
            
            print("Cycling over edges...");
            for( Edge_T e = 0; e < edge_count; ++e )
            {
                if( E_active[e] )
                {
                    for( Int i = 0; i < 2; ++i )
                    {
                        Vertex_T v = edges(e,i);
                        
                        if( !V_active[v] )
                        {
                            if( deleted_vertices_found_in_edges==0 )
                            {
                                eprint(className()+"::SelfCheck: Edge "+ToString(e)+" contains deleted vertex "+ToString(v)+".");
                            }
                            deleted_vertices_found_in_edges++;
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
                            deleted_parent_simplices_found_in_edges++;
                        }
                        else
                        {
                            bool contained = false;
                            for( Int k = 0; k < S_edge_count; ++k )
                            {
                                Edge_T f = Simplex_FindEdge(s, k);
                                if( e == f )
                                {
                                    contained = true;
                                }
                            }
                            
                            if( !contained )
                            {
                                if( simplices_agnostic_of_their_child_edges==0 )
                                {
                                    eprint(className()+"::SelfCheck: Simplex "+ToString(s)+" cannot find its child edge "+ToString(e)+".");
                                    
//                                    S->PrintStats();
//                                    E->PrintStats();
//                                    E->PrintParentSimplices();
                                }
                                simplices_agnostic_of_their_child_edges++;
                            }
                            
                        }
                    }
                }
                
            } // for( Edge_T e = 0; e < edge_count; ++e )
            
            
            print("Cycling over simplices...");
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
                            deleted_vertices_found_in_simplices++;
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
//                                    V->PrintStats();
//                                    V->PrintParentSimplices();
//                                    S->PrintStats();
                                }
                                vertices_agnostic_of_their_parent_simplices++;
                            }
                        }
                    }
                    
                    
                    for( Int j = 0; j < S_vertex_count; ++j )
                    {
                        Edge_T e = Simplex_FindEdge(s,j);
                        
                        if( e < 0 )
                        {
                            if( invalid_edges_found_in_simplices==0 )
                            {
                                eprint(className()+"::SelfCheck: Invalid edge "+ToString(e)+" found in simplex "+ToString(s)+".");
                            }
                            invalid_edges_found_in_simplices++;
                        }
                        else if( !E_active[e] )
                        {
                            if( deleted_edges_found_in_simplices==0 )
                            {
                                eprint(className()+"::SelfCheck: Deleted edge "+ToString(e)+" found in simplex "+ToString(s)+".");
                            }
                            deleted_edges_found_in_simplices++;
                    
                        }
                        else
                        {
                            bool contained = (E_parent_simplices[e].Find(s) >= 0);
                            
                            if( !contained )
                            {
                                if( edges_agnostic_of_their_parent_simplices==0 )
                                {
                                    eprint(className()+"::SelfCheck: edge "+ToString(e)+" does not know about its parent simplex "+ToString(s)+".");
//                                    E->PrintStats();
//                                    E->PrintParentSimplices();
//                                    S->PrintStats();
                                }
                                edges_agnostic_of_their_parent_simplices++;
                            }
                        }
                    }
                }
            }
            
            
            ptoc(className()+"::SelfCheck");
        }
        
        
    protected:
        
        
                    
                               
    public:
        
        virtual CLASS & DownCast() override
        {
            return *this;
        }
        
        virtual const CLASS & DownCast() const override
        {
            return *this;
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        
    private:
        
        static std::string className()
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    };
    
        
} // namespace Repulsor

#undef I
#undef R
                            
#undef CLASS
#undef BASE
