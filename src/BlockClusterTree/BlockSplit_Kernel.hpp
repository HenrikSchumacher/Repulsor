#pragma once

#define CLASS BlockSplit_Kernel
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
        
        CLASS() = default;
        
        virtual ~CLASS() = default;
        
        CLASS(
            const ClusterTree_T & S,
            const ClusterTree_T & T,
            const Real far_theta2_,
            const Real near_theta2_
        )
        :   BASE( S, T )
        ,   S_C_proto           ( S.ClusterPrototype().Clone()          )
        ,   T_C_proto           ( T.ClusterPrototype().Clone()          )
        ,   S_P_proto           ( S.PrimitivePrototype().Clone()        )
        ,   T_P_proto           ( T.PrimitivePrototype().Clone()        )
        ,   S_C_serialized      ( S.ClusterSerializedData().data()      )
        ,   T_C_serialized      ( T.ClusterSerializedData().data()      )
        ,   S_P_serialized      ( S.PrimitiveSerializedData().data()    )
        ,   T_P_serialized      ( T.PrimitiveSerializedData().data()    )
        ,   A                   ( S.PrimitiveAdjacencyMatrix()          )
        ,   far_theta2          ( far_theta2_                           )
        ,   near_theta2         ( near_theta2_                          )
        ,   intersection_theta2 ( near_theta2_                          )
        {}
        
        CLASS( const CLASS & other )
        : BASE( other )
        ,   inter_idx           ( other.inter_idx                       )
        ,   inter_jdx           ( other.inter_jdx                       )
        ,   verynear_idx        ( other.verynear_idx                    )
        ,   verynear_jdx        ( other.verynear_jdx                    )
        ,   near_idx            ( other.near_idx                        )
        ,   near_jdx            ( other.near_jdx                        )
        ,   far_idx             ( other.far_idx                         )
        ,   far_jdx             ( other.far_jdx                         )
        ,   S_C_proto           ( other.S_C_proto->Clone()              )
        ,   T_C_proto           ( other.T_C_proto->Clone()              )
        ,   S_P_proto           ( other.S_P_proto->Clone()              )
        ,   T_P_proto           ( other.T_P_proto->Clone()              )
        ,   S_C_serialized      ( other.S_C_serialized                  )
        ,   T_C_serialized      ( other.T_C_serialized                  )
        ,   S_P_serialized      ( other.S_P_serialized                  )
        ,   T_P_serialized      ( other.T_P_serialized                  )
        ,   A                   ( other.A                               )
        ,   far_theta2          ( other.far_theta2                      )
        ,   near_theta2         ( other.near_theta2                     )
        ,   intersection_theta2 ( other.near_theta2                     )
        {}
        
        friend void swap(CLASS &A, CLASS &B)
        {
            // see https://stackoverflow.com/questions/5695548/public-friend-swap-member-function for details
            using std::swap;

            swap( A.tree_string,            B.tree_string           );
            swap( A.inter_idx,              B.inter_idx             );
            swap( A.inter_jdx,              B.inter_jdx             );
            swap( A.verynear_idx,           B.verynear_idx          );
            swap( A.verynear_jdx,           B.verynear_jdx          );
            swap( A.near_idx,               B.near_idx              );
            swap( A.near_jdx,               B.near_jdx              );
            swap( A.far_idx,                B.far_idx               );
            swap( A.far_jdx,                B.far_jdx               );
            swap( A.S_C_proto,              B.S_C_proto             );
            swap( A.T_C_proto,              B.T_C_proto             );
            swap( A.T_P_proto,              B.T_P_proto             );
            swap( A.S_C_serialized,         B.S_C_serialized        );
            swap( A.T_C_serialized,         B.T_C_serialized        );
            swap( A.S_P_serialized,         B.S_P_serialized        );
            swap( A.T_P_serialized,         B.T_P_serialized        );
            swap( A.A,                      B.A                     );
            swap( A.far_theta2,             B.far_theta2            );
            swap( A.near_theta2,            B.near_theta2           );
            swap( A.intersection_theta2,    B.intersection_theta2   );
        }

        // Copy assignment
        CLASS & operator=(CLASS other) // Intentionally no pass-by-reference here!
        {
            swap(*this, other);

            return *this;
        }
                
//        // Move constructor
//        CLASS( CLASS && other ) noexcept
//        {
//            swap(*this, other);
//        }
//
//        // Move assignment operator
//        CLASS & operator=( CLASS && other ) noexcept
//        {
//            if( this != &other )
//            {
//                swap( *this, other );
//            }
//            return *this;
//        }
        
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
        
        std::shared_ptr<typename ClusterTree_T::BoundingVolume_T> S_C_proto;
        std::shared_ptr<typename ClusterTree_T::BoundingVolume_T> T_C_proto;
        std::shared_ptr<typename ClusterTree_T::Primitive_T>      S_P_proto;
        std::shared_ptr<typename ClusterTree_T::Primitive_T>      T_P_proto;
        
        GJK_Algorithm<ClusterTree_T::AMB_DIM,GJK_Real,Int> G;
        
        SReal * restrict const S_C_serialized;
        SReal * restrict const T_C_serialized;

        SReal * restrict const S_P_serialized;
        SReal * restrict const T_P_serialized;
        
        const SparseBinaryMatrixCSR<Int> & A;
        
        Int C_i = -1;
        Int C_j = -1;
        Int P_i = -1;
        Int P_j = -1;
        
        const Real far_theta2;
        const Real near_theta2;
        const Real intersection_theta2;
        
    public:

        virtual force_inline void LoadClusterS( const Int i ) override
        {
            C_i = i;
            S_C_proto->SetPointer(S_C_serialized,i);
        }
        
        virtual force_inline void LoadClusterT( const Int j ) override
        {
            C_j = j;
            T_C_proto->SetPointer(T_C_serialized,j);
        }
        
        virtual force_inline void LoadPrimitiveS( const Int i ) override
        {
            P_i = i;
            S_P_proto->SetPointer(S_P_serialized,i);
        }
        virtual force_inline void LoadPrimitiveT( const Int j ) override
        {
            P_j = j;
            T_P_proto->SetPointer(T_P_serialized,j);
        }
        
        virtual force_inline bool IsAdmissable() override
        {
            return G.MultipoleAcceptanceCriterion( *S_C_proto, *T_C_proto, far_theta2 );
        }
        
        virtual force_inline SReal ClusterScoreS() override
        {
            return S_C_proto-> SquaredRadius();
        }
        virtual force_inline SReal ClusterScoreT() override
        {
            return T_C_proto->SquaredRadius();
        }
        
        virtual force_inline void ComputeLeafDiagonal() override
        {
            near_idx.push_back(P_i);
            near_jdx.push_back(P_i);
        }
        
        virtual force_inline void ComputeLeaf() override
        {
            const bool neighbor_found = A.FindNonzeroPosition(P_i,P_j) >= 0;

            const bool admissable = neighbor_found || G.MultipoleAcceptanceCriterion(
                    *S_P_proto, *T_P_proto, near_theta2
            );
            
            if( admissable )
            {
                near_idx.push_back(P_i);
                near_jdx.push_back(P_j);
            }
            else
            {
                const bool intersecting = (!G.SeparatedQ()) && G.IntersectingQ(
                    *S_P_proto, *T_P_proto, intersection_theta2
                );
                
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

            const bool admissable = neighbor_found || G.MultipoleAcceptanceCriterion(
                *S_P_proto, *T_P_proto, near_theta2
            );
            
            if( admissable )
            {
                near_idx.push_back(P_j);
                near_jdx.push_back(P_i);
            }
            else
            {
                const bool intersecting = (!G.SeparatedQ()) && G.IntersectingQ(
                    *S_P_proto, *T_P_proto, intersection_theta2
                );
                
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

        
    public:
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        
    private:
      
        std::string className() const
        {
            return TO_STD_STRING(CLASS) + "<"+this->tree_string+">";
        }
    };
    
} //namespace Repulsor

#undef BASE
#undef CLASS

