#pragma once

#define CLASS BlockSplitter
#define BASE  Traversor_Kernel<ClusterTree_T_>

namespace Repulsor
{
    template<typename ClusterTree_T_>
    class CLASS : public BASE
    {
    public:
        
        using ClusterTree_T = typename BASE::ClusterTree_T;
        using Real    = typename BASE::Real;
        using Int     = typename BASE::Int;
        using SReal   = typename BASE::SReal;
        using ExtReal = typename BASE::ExtReal;
        
    public:
        
        virtual ~CLASS() = default;
        
        CLASS(
            const ClusterTree_T & S,
            const ClusterTree_T & T,
            const Real far_theta2_,
            const Real near_theta2_
        )
        :   BASE( S, T )
        ,   S_C_proto_ptr  ( S.ClusterPrototype().Clone() )
        ,   T_C_proto_ptr  ( T.ClusterPrototype().Clone() )
        ,   S_P_proto_ptr  ( S.PrimitivePrototype().Clone() )
        ,   T_P_proto_ptr  ( T.PrimitivePrototype().Clone() )
        ,   S_C_serialized ( S.ClusterSerializedData().data() )
        ,   T_C_serialized ( T.ClusterSerializedData().data() )
        ,   S_P_serialized ( S.PrimitiveSerializedData().data() )
        ,   T_P_serialized ( T.PrimitiveSerializedData().data() )
        ,   A ( S.PrimitiveAdjacencyMatrix() )
        ,   far_theta2  ( far_theta2_ )
        ,   near_theta2 ( near_theta2_ )
        ,   intersection_theta2 ( near_theta2_ )
        {
            const Int expected = static_cast<Int>(10) * ( S.PrimitiveCount() + T.PrimitiveCount() );
            
            near_idx.reserve(expected);
            near_jdx.reserve(expected);
            far_idx.reserve(expected);
            far_jdx.reserve(expected);
        }
        
    public:
        
        std::vector<Int>    inter_idx;
        std::vector<Int>    inter_jdx;
        std::vector<Int> verynear_idx;
        std::vector<Int> verynear_jdx;
        std::vector<Int>     near_idx;
        std::vector<Int>     near_jdx;
        std::vector<Int>      far_idx;
        std::vector<Int>      far_jdx;

    protected:
        
        std::shared_ptr<typename ClusterTree_T::BoundingVolume_T> S_C_proto_ptr;
        std::shared_ptr<typename ClusterTree_T::BoundingVolume_T> T_C_proto_ptr;
        std::shared_ptr<typename ClusterTree_T::Primitive_T>      S_P_proto_ptr;
        std::shared_ptr<typename ClusterTree_T::Primitive_T>      T_P_proto_ptr;
        
        typename ClusterTree_T::BoundingVolume_T & S_C_proto = *S_C_proto_ptr;
        typename ClusterTree_T::BoundingVolume_T & T_C_proto = *T_C_proto_ptr;
        typename ClusterTree_T::Primitive_T      & S_P_proto = *S_P_proto_ptr;
        typename ClusterTree_T::Primitive_T      & T_P_proto = *T_P_proto_ptr;
        
        typename ClusterTree_T::GJK_T G;

        SReal * restrict const S_C_serialized;
        SReal * restrict const T_C_serialized;

        SReal * restrict const S_P_serialized;
        SReal * restrict const T_P_serialized;
        
        SparseBinaryMatrixCSR<Int> & A;
        
        Int C_i = 0;
        Int C_j = 0;
        Int P_i = 0;
        Int P_j = 0;
        
        const Real far_theta2;
        const Real near_theta2;
        const Real intersection_theta2;
        
    public:

        virtual force_inline void LoadClusterS( const Int i ) override
        {
            C_i = i;
            S_C_proto.SetPointer(S_C_serialized,i);
        }
        
        virtual force_inline void LoadClusterT( const Int j ) override
        {
            C_j = j;
            T_C_proto.SetPointer(T_C_serialized,j);
        }
        
        virtual force_inline void LoadPrimitiveS( const Int i ) override
        {
            P_i = i;
            S_P_proto.SetPointer(S_C_serialized,i);
        }
        virtual force_inline void LoadPrimitiveT( const Int j ) override
        {
            P_j = j;
            T_P_proto.SetPointer(T_C_serialized,j);
        }
        
        virtual force_inline bool IsAdmissable() override
        {
            return G.MultipoleAcceptanceCriterion( S_C_proto, T_C_proto, far_theta2 );
        }
        
        virtual force_inline bool ClusterScoreS() override
        {
            return S_C_proto.SquaredRadius();
        }
        virtual force_inline bool ClusterScoreT() override
        {
            return T_C_proto.SquaredRadius();
        }
        
        virtual force_inline void ComputeLeafDiagonal() override
        {
            near_idx.push_back(P_i);
            near_jdx.push_back(P_i);
        }
        
        virtual force_inline void ComputeLeaf() override
        {
            const bool neighbor_found = A.FindNonzeroPosition(P_i,P_j) >= 0;

            const bool admissable = neighbor_found || G.MultipoleAcceptanceCriterion( S_P_proto, T_P_proto, near_theta2);
            
            if( admissable )
            {
                near_idx.push_back(P_i);
                near_jdx.push_back(P_j);
            }
            else
            {
                const bool intersecting = (!G.SeparatedQ()) && G.IntersectingQ( S_P_proto, T_P_proto, intersection_theta2 );
                
                if( intersecting )
                {
                    inter_idx.push_back(P_i);
                    inter_jdx.push_back(P_j);
                }
                else
                {
                    verynear_idx.push_back(P_i);
                    verynear_jdx.push_back(P_j);
                }
            }
        }
        
        virtual force_inline void ComputeLeafSwapped() override
        {
            const bool neighbor_found = A.FindNonzeroPosition(P_i,P_j) >= 0;

            const bool admissable = neighbor_found || G.MultipoleAcceptanceCriterion( S_P_proto, T_P_proto, near_theta2);
            
            if( admissable )
            {
                near_idx.push_back(P_j);
                near_jdx.push_back(P_i);
            }
            else
            {
                const bool intersecting = (!G.SeparatedQ()) && G.IntersectingQ( S_P_proto, T_P_proto, intersection_theta2 );
                
                if( intersecting )
                {
                    inter_idx.push_back(P_j);
                    inter_jdx.push_back(P_i);
                }
                else
                {
                    verynear_idx.push_back(P_j);
                    verynear_jdx.push_back(P_i);
                }
            }
        }
        
        virtual force_inline void ComputeAdmissable() override
        {
            far_idx.push_back(C_i);
            far_jdx.push_back(C_j);
        }
        virtual force_inline void ComputeAdmissableSwapped() override
        {
            far_idx.push_back(C_j);
            far_jdx.push_back(C_i);
        }

        
    protected:
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        
    private:
      
        std::string className()
        {
            return TO_STD_STRING(CLASS) + "<"+this->tree_string+">";
        }
    };
    
} //namespace Repulsor

#undef BASE
#undef CLASS

