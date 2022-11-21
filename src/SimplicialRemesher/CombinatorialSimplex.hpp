#pragma once

namespace Repulsor
{
    
    template<int DOM_DIM, typename Int> class CombinatorialVertex;
    template<int DOM_DIM, typename Int> class CombinatorialEdge;
    template<int DOM_DIM, typename Int> class CombinatorialSimplex;
    
    template<int DOM_DIM, typename Int>
    class CombinatorialSimplex
    {
    public:
        
        using Vertex_T  = CombinatorialVertex <DOM_DIM,Int>;
        using Edge_T    = CombinatorialEdge   <DOM_DIM,Int>;
        using Simplex_T = CombinatorialSimplex<DOM_DIM,Int>;
        
//        using Pair_T = std::pair<Vertex_T *,Vertex_T *>;
        
        using Pair_T = std::pair<Int,Int>;
        
        static constexpr Int edge_count = ((DOM_DIM+1)*DOM_DIM)/2;
        
        // List of {i,j} pairs for the edges, where i, j are simplex-local indices.
        static Int idx [edge_count];
        static Int jdx [edge_count];
        
    protected:

        std::array<Vertex_T *,DOM_DIM+1> vertices;
        
        Int identifier = -1;
        
        bool deleted = false;
        
    public:
        
        CombinatorialSimplex() = default;
        
        ~CombinatorialSimplex() = default;

        CombinatorialSimplex( const Int identifier_, std::array<Vertex_T *,DOM_DIM+1> & vlist )
        :   vertices(vlist)
        ,   identifier(identifier_)
        ,   deleted(false)
        {}

        CombinatorialSimplex( const Int identifier_, Vertex_T ** vlist )
        :   identifier(identifier_)
        ,   deleted(false)
        {
            Vertex_T ** restrict w = vlist;
            
            for( Int i = 0; i < DOM_DIM+1; ++i )
            {
                vertices[i] = w[i];
            }
        }
        
        CombinatorialSimplex( const Int identifier_, Vertex_T ** vlist, const Int * const index_list )
        :   identifier(identifier_)
        ,   deleted(false)
        {
//            ptic(ClassName()+"::CombinatorialSimplex");

            Vertex_T ** restrict w = vlist;
            const Int * restrict const s = index_list;
            
            for( Int i = 0; i < DOM_DIM+1; ++i )
            {
                vertices[i] = w[s[i]];
            }
//            ptoc(ClassName()+"::CombinatorialSimplex");
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
        
        static constexpr Int VertexCount()
        {
            return DOM_DIM+1;
        }
        
        Vertex_T * Vertex( const Int i )
        {
            // cppcheck-suppress CastIntegerToAddressAtReturn
            return vertices[i];
        }
        
        const Vertex_T * Vertex( const Int i ) const
        {
            // cppcheck-suppress CastIntegerToAddressAtReturn
            return vertices[i];
        }
        
        void SetVertex( const Int i, Vertex_T * v )
        {
            vertices[i] = v;
        }
        
        std::array<Vertex_T *,DOM_DIM+1> & Vertices()
        {
            return vertices;
        }
        
        const std::array<Vertex_T *,DOM_DIM+1> & Vertices() const
        {
            return vertices;
        }
        
        static constexpr Int EdgeCount()
        {
            return edge_count;
        }
        
        void ReplaceVertex( const Vertex_T * const  V, Vertex_T * const W )
        {
            for( Int i = 0; i < DOM_DIM+1; ++i )
            {
                if( vertices[i] == V )
                {
                    vertices[i] = W;
                }
            }
        }
        


        void OppositeVertices( Edge_T * E, std::array<Vertex_T *, DOM_DIM+1> & opp_vertex_list ) const
        {
            Int counter = 0;

            const Vertex_T * V_0 = E->Vertex(0);
            const Vertex_T * V_1 = E->Vertex(1);

            for( Vertex_T * V : vertices )
            {
                if( V != V_0 && V!=V_1 )
                {
                    opp_vertex_list[counter] = V;
                    ++counter;
                }
            }
        }
        
        
//        friend std::string ToString( const Simplex_T * S )
//        {
//            if( S == nullptr )
//            {
//                return std::string( "nullptr" );
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
            std::stringstream s;
            s << ClassName() <<" = { ID = " << ID() <<", {";
            
            s << Tools::ToString(vertices[0]);
            
            const Int n = static_cast<Int>(vertices.size());
            
            for( Int i = 1; i < n; ++i )
            {
                s << ",";
                
                //  cppcheck-suppress containerOutOfBounds
                s << Tools::ToString( vertices[i] );
            }
            
            s << "} }";
            
            return s.str();
        }
        
        void PrintStats() const
        {
            print(ToString());
        }
        
    public:
        
        std::string ClassName() const
        {
            return "CombinatorialSimplex<"+Tools::ToString(DOM_DIM)+","+TypeName<Int>::Get()+">";
        }
        
    }; // CombinatorialSimplex
    
    template<>
    int CombinatorialSimplex<0,int>::idx [0] = {};
    template<>
    int CombinatorialSimplex<0,int>::jdx [0] = {};
    
    template<>
    int CombinatorialSimplex<1,int>::idx [1] = {0};
    template<>
    int CombinatorialSimplex<1,int>::jdx [1] = {1};
    
    template<>
    int CombinatorialSimplex<2,int>::idx [3] = {0,0,1};
    template<>
    int CombinatorialSimplex<2,int>::jdx [3] = {1,2,2};
    
    template<>
    int CombinatorialSimplex<3,int>::idx [6] = {0,0,0,1,1,2};
    template<>
    int CombinatorialSimplex<3,int>::jdx [6] = {1,2,3,2,3,3};
    
    template<>
    int CombinatorialSimplex<4,int>::idx [10] = {0,0,0,0,1,1,1,2,3,4};
    template<>
    int CombinatorialSimplex<4,int>::jdx [10] = {1,2,3,4,2,3,4,3,4,4};
    
    
    
    template<>
    long long CombinatorialSimplex<0,long long>::idx [0] = {};
    template<>
    long long CombinatorialSimplex<0,long long>::jdx [0] = {};
    
    template<>
    long long CombinatorialSimplex<1,long long>::idx [1] = {0};
    template<>
    long long CombinatorialSimplex<1,long long>::jdx [1] = {1};
    
    template<>
    long long CombinatorialSimplex<2,long long>::idx [3] = {0,0,1};
    template<>
    long long CombinatorialSimplex<2,long long>::jdx [3] = {1,2,2};
    
    template<>
    long long CombinatorialSimplex<3,long long>::idx [6] = {0,0,0,1,1,2};
    template<>
    long long CombinatorialSimplex<3,long long>::jdx [6] = {1,2,3,2,3,3};
    
    template<>
    long long CombinatorialSimplex<4,long long>::idx [10] = {0,0,0,0,1,1,1,2,3,4};
    template<>
    long long CombinatorialSimplex<4,long long>::jdx [10] = {1,2,3,4,2,3,4,3,4,4};
    
    
//    template<typename Int>
//    Int CombinatorialSimplex<0,Int>::idx [0] = {};
//    template<typename Int>
//    Int CombinatorialSimplex<0,Int>::jdx [0] = {};
//
//    template<typename Int>
//    Int CombinatorialSimplex<1,Int>::idx [1] = {0};
//    template<typename Int>
//    Int CombinatorialSimplex<1,Int>::jdx [1] = {1};
//
//    template<typename Int>
//    Int CombinatorialSimplex<2,Int>::idx [3] = {0,0,1};
//    template<typename Int>
//    Int CombinatorialSimplex<2,Int>::jdx [3] = {1,2,2};
//
//    template<typename Int>
//    Int CombinatorialSimplex<3,Int>::idx [6] = {0,0,0,1,1,2};
//    template<typename Int>
//    Int CombinatorialSimplex<3,Int>::jdx [6] = {1,2,3,2,3,3};
//
//    template<typename Int>
//    Int CombinatorialSimplex<4,Int>::idx [10] = {0,0,0,0,1,1,1,2,3,4};
//    template<typename Int>
//    Int CombinatorialSimplex<4,Int>::jdx [10] = {1,2,3,4,2,3,4,3,4,4};
    
} // namespace Repulsor
