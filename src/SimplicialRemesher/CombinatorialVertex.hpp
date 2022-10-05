#pragma once

namespace Repulsor
{

    template<int DOM_DIM, typename Int> class CombinatorialVertex;
    template<int DOM_DIM, typename Int> class CombinatorialEdge;
    template<int DOM_DIM, typename Int> class CombinatorialSimplex;
    
    template<int DOM_DIM, typename Int>
    class CombinatorialVertex
    {
    public:
        
        using Vertex_T  = CombinatorialVertex <DOM_DIM,Int>;
        using Edge_T    = CombinatorialEdge   <DOM_DIM,Int>;
        using Simplex_T = CombinatorialSimplex<DOM_DIM,Int>;
        
        static constexpr Int max_simplex_valence = (DOM_DIM == 1) ? 3 : (DOM_DIM == 2) ? 9 : 100;
        
    protected:
        
        Int identifier = -1;
        
        bool deleted = false;
        
        bool modified = false;
        
        SortedList<Simplex_T *,Int> parent_simplices;
        
    public:
        
        ~CombinatorialVertex() = default;
        
//        CombinatorialVertex() : deleted(true)
//        {}
        
        explicit CombinatorialVertex( const Int ID_ ) : identifier(ID_), deleted(false)
        {
            switch (DOM_DIM)
            {
                case 1:
                {
                    parent_simplices.Reserve(2);
                    break;
                }
                case 2:
                {
                    parent_simplices.Reserve(8);
                    break;
                }
                case 3:
                {
                    parent_simplices.Reserve(24);
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
        
        bool ModifiedQ() const
        {
            return modified;
        }
        
        void MarkAsModified()
        {
            modified = true;
        }
        
        Int ParentSimplexCount() const
        {
            return parent_simplices.Size();
        }
        
        Simplex_T * ParentSimplex( const Int k )
        {
            return parent_simplices[k];
        }
        
        const SortedList<Simplex_T *,Int> & ParentSimplices() const
        {
            return parent_simplices;
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
        
        
        void CollectNeighborVertices( SortedList<Vertex_T *,Int> & neighbors ) const
        {
//            ptic(ClassName()+"NeighborVertices");
            
            neighbors.Clear();
            
            for( Simplex_T * S : parent_simplices )
            {
//                if( S == nullptr )
//                {
//                    eprint("!!!");
//                    PrintParentSimplices();
//                    break;
//                }
                
//                print(S->ToString());
                
                for( Vertex_T * V : S->Vertices() )
                {
//                    if( V == nullptr )
//                    {
//                        eprint("!!");
//                    }
//                    else
//                    {
//                        valprint("V->ID()",V->ID());
//                    }

                    if( V->ID() != this->ID() )
                    {
//                        print(V->ToString());
//                        Int i = neighbors.FindPosition(V);
//                        valprint("i",i);
                        
                        
                        neighbors.Insert(V);
                    }
                }
            
//                print(neighbors.ToString());
            }

//            ptoc(ClassName()+"NeighborVertices");
        }
        
//        friend std::string ToString( const Vertex_T * S )
//        {
//            if( S == nullptr )
//            {
//                return std::string( "NULL" );
//            }
//            {
//                if( S->DeletedQ() )
//                {
//                    return std::string( "DEL" );
//                }
//                else
//                {
//                    return std::to_string(S->ID());
//                }
//            }
//        }
        
        std::string ToString() const
        {
            return ClassName()+" = { ID = "+Tools::ToString(ID())+", {"+Tools::ToString(this)+"} }";
        }
        
        void PrintStats() const
        {
            print(ToString());
        }
        
        void PrintParentSimplices() const
        {
            std::stringstream s;
            s<< ClassName()+"::PrintParentSimplices() = { ID = "+Tools::ToString(ID())+", "+parent_simplices.ToString()+" }.";
            print(s.str());
        }
        
    public:
        
        std::string ClassName() const
        {
            return "CombinatorialVertex<"+Tools::ToString(DOM_DIM)+","+TypeName<Int>::Get()+">";
        }
        
    }; // CombinatorialVertex
    
} // namespace Repulsor
