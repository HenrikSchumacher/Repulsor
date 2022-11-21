#pragma once

namespace Repulsor
{
    
    template<int DOM_DIM, typename Int> class CombinatorialVertex;
    template<int DOM_DIM, typename Int> class CombinatorialEdge;
    template<int DOM_DIM, typename Int> class CombinatorialSimplex;
    
    template<int DOM_DIM, typename Int>
    class CombinatorialEdge
    {
    protected:
        
        using Vertex_T  = CombinatorialVertex <DOM_DIM,Int>;
        using Edge_T    = CombinatorialEdge   <DOM_DIM,Int>;
        using Simplex_T = CombinatorialSimplex<DOM_DIM,Int>;
        
        using Pair_T = std::pair<Int,Int>;
        
        Int identifier = -1;
        
        bool deleted = false;
        
        std::array<Vertex_T *,2> vertices = {nullptr,nullptr};
        SortedList<Simplex_T *,Int> parent_simplices;
        
    public:
        
        ~CombinatorialEdge() = default;
        
        CombinatorialEdge()
        :   deleted(true)
        {}
        
        CombinatorialEdge( const Int ID_, Vertex_T * & V_0, Vertex_T * & V_1 )
        :   identifier(ID_)
        ,   deleted(false)
        {
            vertices[0] = V_0;
            vertices[1] = V_1;
            
            switch (DOM_DIM)
            {
                case 1:
                {
                    parent_simplices.Reserve(1);
                    break;
                }
                case 2:
                {
                    parent_simplices.Reserve(2);
                    break;
                }
                case 3:
                {
                    parent_simplices.Reserve(6);
                    break;
                }
                default:
                {
                    break;
                }
            }
        }
        
        Int ID() const
        {
            return identifier;
        }
        
        void SetID( Int ID_ )
        {
            identifier = ID_;
        }
        
        bool DeletedQ() const
        {
            return deleted;
        }

        void MarkAsDeleted()
        {
            deleted = true;
        }
        
        Int ParentSimplexCount() const
        {
            return parent_simplices.Size();
        }
        
        Simplex_T * ParentSimplex( const Int k )
        {
            return parent_simplices[k];
        }
        
        void InsertParentSimplex( Simplex_T * S )
        {
            parent_simplices.Insert(S);
        }
        
        void DropParentSimplex( Simplex_T * S )
        {
            parent_simplices.Drop(S);
        }
        
        void ReplaceParentSimplex( Simplex_T * S, Simplex_T * T )
        {
            parent_simplices.Drop(S);
            parent_simplices.Insert(T);
        }
        
        Vertex_T * Vertex( const Int i )
        {
            // cppcheck-suppress CastIntegerToAddressAtReturn
            return vertices[i];
        }
        
        std::array<Vertex_T *,2> & Vertices()
        {
            return vertices;
        }
        
        void SetVertex( const Int i, Vertex_T * v )
        {
            vertices[i] = v;
        }
        
        SortedList<Simplex_T *,Int> & ParentSimplices()
        {
            return parent_simplices;
        }
        
        const SortedList<Simplex_T *,Int> & ParentSimplices() const
        {
            return parent_simplices;
        }
        
        void ReplaceVertex( const Vertex_T * V, Vertex_T * W )
        {
            if( vertices[0] == V )
            {
                vertices[0] = W;
            }
            if( vertices[1] == V )
            {
                vertices[1] = W;
            }
        }
        
        
        void CollectOpposingVertices( SortedList<Vertex_T *, Int> & opp_vertices ) const
        {
//            ptic(ClassName()+"::CollectOpposingVertices");
            
            opp_vertices.Clear();
            
            const Vertex_T * V_0 = vertices[0];
            const Vertex_T * V_1 = vertices[1];
            
            // Going through the simplices to find opposing vertices.
            for( Simplex_T * S : parent_simplices )
            {
//                print(S->ToString());
                for( Vertex_T * V : S->Vertices() )
                {
                    if( V != V_0 && V!=V_1 )
                    {
                        opp_vertices.Insert(V);
                    }
                }
            }
            
//            ptoc(ClassName()+"::CollectOpposingVertices");
        }
        
//        friend std::string ToString( const Edge_T * E )
//        {
//            if( E == nullptr )
//            {
//                return std::string( "nullptr" );
//            }
//            {
//                if( E->DeletedQ() )
//                {
//                    return std::string( "DEL" );
//                }
//                else
//                {
//                    return std::to_string(E->ID());
//                }
//            }
//        }
        

        
        std::string ToString() const
        {
            std::stringstream s;
            s << ClassName() << " = { ID = " << ID() <<", {";
            s << Tools::ToString(vertices[0]);
            s << ", ";
            s << Tools::ToString(vertices[1]);
            s << "} }";
            
            return s.str();
        }

        void PrintStats() const
        {
            print(ToString());
        }
        
        void PrintParentSimplices() const
        {
            std::stringstream s;
            s<< ClassName()+"::PrintParentSimplices() = { ID = "+Tools::ToString(ID())+", "+ParentSimplices().ToString()+" }.";
            print(s.str());
        }
        
    public:
        
        std::string ClassName() const
        {
            return "CombinatorialEdge<"+Tools::ToString(DOM_DIM)+","+TypeName<Int>::Get()+">";
        }
        
    }; // CombinatorialEdge
    
} // namespace Repulsor
