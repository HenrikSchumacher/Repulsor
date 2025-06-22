#pragma once

namespace Repulsor
{
    template<typename ClusterTree_T_>
    class BlockSplit_Kernel
    {
    public:
        
        using ClusterTree_T = ClusterTree_T_;
        using Real          = typename ClusterTree_T::Real;
        using SReal         = typename ClusterTree_T::SReal;
        using ExtReal       = typename ClusterTree_T::ExtReal;
        using Int           = typename ClusterTree_T::Int;
        using LInt          = typename ClusterTree_T::LInt;
        
    public:

        
        BlockSplit_Kernel(
            cref<ClusterTree_T> S, cref<ClusterTree_T> T,
            const Int thread_,
            const Real far_theta2_,
            const Real near_theta2_
        )
        :   inter_idx           ( 1 )
        ,   verynear_idx        ( 1 )
        ,   near_idx            ( 2 * ( S.PrimitiveCount() + T.PrimitiveCount() ) )
        ,   far_idx             ( 2 * ( S.PrimitiveCount() + T.PrimitiveCount() ) )
        ,   thread              ( thread_ )
        ,   S_C_proto           ( S.ClusterPrototype().Clone()              )
        ,   T_C_proto           ( T.ClusterPrototype().Clone()              )
        ,   S_P_proto           ( S.PrimitivePrototype().Clone()            )
        ,   T_P_proto           ( T.PrimitivePrototype().Clone()            )
        ,   S_C_serialized      ( S.ClusterSerialized().data()              )
        ,   T_C_serialized      ( T.ClusterSerialized().data()              )
        ,   S_P_serialized      ( S.PrimitiveSerialized().data()            )
        ,   T_P_serialized      ( T.PrimitiveSerialized().data()            )
        ,   A                   ( S.PrimitiveAdjacencyMatrix()              )
        ,   far_theta2          ( far_theta2_                               )
        ,   near_theta2         ( near_theta2_                              )
        ,   intersection_theta2 ( near_theta2_                              )
        {}
        
        // Default constructor
        BlockSplit_Kernel() = delete;
        // Destructor
        ~BlockSplit_Kernel() = default;
        // Copy constructor
        BlockSplit_Kernel( const BlockSplit_Kernel & other )
        :   inter_idx           ( other.inter_idx.Capacity()    )
        ,   verynear_idx        ( other.verynear_idx.Capacity() )
        ,   near_idx            ( other.near_idx.Capacity()     )
        ,   far_idx             ( other.far_idx.Capacity()      )
        ,   thread              ( other.thread                  )
        ,   S_C_proto           ( other.S_C_proto->Clone()      ) // !!!
        ,   T_C_proto           ( other.T_C_proto->Clone()      ) // !!!
        ,   S_P_proto           ( other.S_P_proto->Clone()      ) // !!!
        ,   T_P_proto           ( other.T_P_proto->Clone()      ) // !!!
        ,   S_C_serialized      ( other.S_C_serialized          )
        ,   T_C_serialized      ( other.T_C_serialized          )
        ,   S_P_serialized      ( other.S_P_serialized          )
        ,   T_P_serialized      ( other.T_P_serialized          )
        ,   A                   ( other.A                       )
        ,   far_theta2          ( other.far_theta2              )
        ,   near_theta2         ( other.near_theta2             )
        ,   intersection_theta2 ( other.near_theta2             )
        {}
        
        // Swap function
        friend void swap( BlockSplit_Kernel & X, BlockSplit_Kernel & Y ) noexcept
        {
            // see https://stackoverflow.com/questions/5695548/public-friend-swap-member-function for details
            using std::swap;
            
            swap( X.inter_idx,              Y.inter_idx             );
            swap( X.verynear_idx,           Y.verynear_idx          );
            swap( X.near_idx,               Y.near_idx              );
            swap( X.far_idx,                Y.far_idx               );
            swap( X.thread,                 Y.thread                );
            swap( X.S_C_proto,              Y.S_C_proto             );
            swap( X.T_C_proto,              Y.T_C_proto             );
            swap( X.T_P_proto,              Y.T_P_proto             );
            swap( X.S_C_serialized,         Y.S_C_serialized        );
            swap( X.T_C_serialized,         Y.T_C_serialized        );
            swap( X.S_P_serialized,         Y.S_P_serialized        );
            swap( X.T_P_serialized,         Y.T_P_serialized        );
            swap( X.A,                      Y.A                     );
            swap( X.far_theta2,             Y.far_theta2            );
            swap( X.near_theta2,            Y.near_theta2           );
            swap( X.intersection_theta2,    Y.intersection_theta2   );
        }

        // Copy assignment
        BlockSplit_Kernel & operator=( BlockSplit_Kernel other ) noexcept
        {
            swap(*this, other);
            return *this;
        }
        
        // Let's better not talk about the move constructor and move assignment here.
        // Since this kernel is put into a std::vector, some move semantics are necessary. We let the compiler let them generate for us. The data in this kernel are mostly pointers and references, so we we are happy if copy semantics are used throughout.
        
    public:
        
        PairAggregator<Int,Int,LInt> inter_idx;
        PairAggregator<Int,Int,LInt> verynear_idx;
        PairAggregator<Int,Int,LInt> near_idx;
        PairAggregator<Int,Int,LInt> far_idx;

    protected:
        
        Int thread;
        
        std::shared_ptr<typename ClusterTree_T::BoundingVolume_T> S_C_proto;
        std::shared_ptr<typename ClusterTree_T::BoundingVolume_T> T_C_proto;
        std::shared_ptr<typename ClusterTree_T::Primitive_T>      S_P_proto;
        std::shared_ptr<typename ClusterTree_T::Primitive_T>      T_P_proto;
        
        GJK<ClusterTree_T::AMB_DIM,GJK_Real,Int> G;
        
        mptr<SReal> S_C_serialized = nullptr;
        mptr<SReal> T_C_serialized = nullptr;

        mptr<SReal> S_P_serialized = nullptr;
        mptr<SReal> T_P_serialized = nullptr;
        
        cref<typename ClusterTree_T::SparseBinaryMatrix_T> A;
        
        Int C_i = -1;
        Int C_j = -1;
        Int P_i = -1;
        Int P_j = -1;
        
        const Real far_theta2;
        const Real near_theta2;
        const Real intersection_theta2;
        
    public:

        TOOLS_FORCE_INLINE void LoadClusterS( const Int i )
        {
            C_i = i;
            S_C_proto->SetPointer( S_C_serialized, i );
        }
        
        TOOLS_FORCE_INLINE void LoadClusterT( const Int j )
        {
            C_j = j;
            T_C_proto->SetPointer( T_C_serialized, j );
        }
        
        TOOLS_FORCE_INLINE void LoadPrimitiveS( const Int i )
        {
            P_i = i;
            S_P_proto->SetPointer( S_P_serialized, i );
        }
        TOOLS_FORCE_INLINE void LoadPrimitiveT( const Int j )
        {
            P_j = j;
            T_P_proto->SetPointer( T_P_serialized, j );
        }
        
        TOOLS_FORCE_INLINE bool AdmissableQ()
        {
            const bool result = G.MultipoleAcceptanceCriterion( *S_C_proto, *T_C_proto, far_theta2 );
            
            return result;
        }
        
        TOOLS_FORCE_INLINE SReal ClusterScoreS()
        {
            return S_C_proto->SquaredRadius();
        }
        
        TOOLS_FORCE_INLINE SReal ClusterScoreT()
        {
            return T_C_proto->SquaredRadius();
        }
        
        TOOLS_FORCE_INLINE void ComputeLeafDiagonal()
        {
//            near_idx.Push(P_i,P_i);
        }
        
        TOOLS_FORCE_INLINE void ComputeLeaf()
        {
            const bool neighbor_found = A.NonzeroPositionQ(P_i,P_j);
            
            const bool admissableQ = neighbor_found || G.MultipoleAcceptanceCriterion(
                    *S_P_proto, *T_P_proto, near_theta2
                );
            
            if( admissableQ )
            {
                near_idx.Push(P_i,P_j);
            }
            else
            {
                const bool intersectingQ = (!G.SeparatedQ()) && G.IntersectingQ(
                    *S_P_proto, *T_P_proto, intersection_theta2, true
                );
                
                if( intersectingQ )
                {
                    inter_idx.Push(P_i,P_j);
                }
                else
                {
                    verynear_idx.Push(P_i,P_j);
                }
            }
        }
        
        TOOLS_FORCE_INLINE void ComputeLeafSwapped()
        {
            const bool neighbor_found = A.NonzeroPositionQ(P_i,P_j);
            
            const bool admissableQ = neighbor_found || G.MultipoleAcceptanceCriterion(
                *S_P_proto, *T_P_proto, near_theta2
            );
            
            if( admissableQ )
            {
                near_idx.Push(P_j,P_i);
            }
            else
            {
                const bool intersectingQ = (!G.SeparatedQ()) && G.IntersectingQ(
                    *S_P_proto, *T_P_proto, intersection_theta2, true
                );
                
                if( intersectingQ )
                {
                    inter_idx.Push(P_j,P_i);
                }
                else
                {
                    verynear_idx.Push(P_j,P_i);
                }
            }
        }
        
        TOOLS_FORCE_INLINE void ComputeAdmissable()
        {
            far_idx.Push(C_i,C_j);
        }
        TOOLS_FORCE_INLINE void ComputeAdmissableSwapped()
        {
            far_idx.Push(C_j,C_i);
        }
        
    private:
      
        static std::string className()
        {
            return std::string("BlockSplit_Kernel<...>");
        }
        
    public:
        
        std::string ClassName() const
        {
            return className();
        }
        
    }; // class BlockSplit_Kernel
    
} // namespace Repulsor

