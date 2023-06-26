#pragma once

namespace Repulsor
{
    template<typename ClusterTree_T_>
    class Collision_Kernel : public ClusterTreePairTraversor_Kernel<ClusterTree_T_>
    {
    private:
        
        using Base_T = ClusterTreePairTraversor_Kernel<ClusterTree_T_>;
        
    public:
        
        using ClusterTree_T = ClusterTree_T_;
        using Real    = typename Base_T::Real;
        using Int     = typename Base_T::Int;
        using SReal   = typename Base_T::SReal;
        using ExtReal = typename Base_T::ExtReal;
        using LInt    = typename Base_T::LInt;
        
    public:
        
        Collision_Kernel() = delete;
        
        ~Collision_Kernel() = default;
        
        Collision_Kernel(
            const ClusterTree_T & S,
            const ClusterTree_T & T,
            const SReal t_init_,
            const SReal TOL
        )
        :   Base_T( S, T )
        ,   S_C_ser     ( S.ClusterSerialized()              )
        ,   S_C_up_ser  ( S.ClusterUpdatedSerialized()       )
        ,   T_C_ser     ( T.ClusterSerialized()              )
        ,   T_C_up_ser  ( T.ClusterUpdatedSerialized()       )
        ,   S_P_ser     ( S.PrimitiveSerialized()            )
        ,   S_P_v_ser   ( S.PrimitiveVelocitiesSerialized()  )
        ,   T_P_ser     ( T.PrimitiveSerialized()            )
        ,   T_P_v_ser   ( T.PrimitiveVelocitiesSerialized()  )
        ,   A           ( S.PrimitiveAdjacencyMatrix()       )
        ,   C           ( S.MovingPrimitivePrototype(), T.MovingPrimitivePrototype(), TOL )
        ,   t_init      ( t_init_ )
        ,   t_max       ( t_init_ )
        {}
        

        Collision_Kernel( const Collision_Kernel & other )
        :   Base_T( other )
        ,   S_C_ser     ( other.S_C_ser     )
        ,   S_C_up_ser  ( other.S_C_up_ser  )
        ,   T_C_ser     ( other.T_C_ser     )
        ,   T_C_up_ser  ( other.T_C_up_ser  )
        ,   S_P_ser     ( other.S_P_ser     )
        ,   S_P_v_ser   ( other.S_P_v_ser   )
        ,   T_P_ser     ( other.T_P_ser     )
        ,   T_P_v_ser   ( other.T_P_v_ser   )
        ,   A           ( other.A           )
        ,   C           ( other.C           )
        ,   t_init      ( other.t_init      )
        ,   t_max       ( other.t_max       )
        ,   triples     ( other.triples     )
        {}
        
        friend void swap(Collision_Kernel &A, Collision_Kernel &B)
        {
            // see https://stackoverflow.com/questions/5695548/public-friend-swap-member-function for details
            using std::swap;

            swap( A.S_C_ser,        B.S_C_ser       );
            swap( A.S_C_up_ser,     B.S_C_up_ser    );
            swap( A.T_C_ser,        B.T_C_ser       );
            swap( A.T_C_up_ser,     B.T_C_up_ser    );
            swap( A.S_P_ser,        B.S_P_ser       );
            swap( A.S_P_v_ser,      B.S_P_v_ser     );
            swap( A.T_P_ser,        B.T_P_ser       );
            swap( A.T_P_v_ser,      B.T_P_v_ser     );
            swap( A.A,              B.A             );
            swap( A.C,              B.C             );
            swap( A.t_init,         B.t_init        );
            swap( A.t_max,          B.t_max         );
            swap( A.triples,        B.triples       );
        }

        // Copy assignment
        Collision_Kernel & operator=(Collision_Kernel other) // Intentionally no pass-by-reference here!
        {
            swap(*this, other);

            return *this;
        }
                
//        // Move constructor
//        Collision_Kernel( Collision_Kernel && other ) noexcept
//        {
//            swap(*this, other);
//        }
//
//        // Move assignment operator
//        Collision_Kernel & operator=( Collision_Kernel && other ) noexcept
//        {
//            if( this != &other )
//            {
//                swap( *this, other );
//            }
//            return *this;
//        }
        
    public:
        
    protected:
        
        const Tensor2<SReal,Int> & S_C_ser;
        const Tensor2<SReal,Int> & S_C_up_ser;
        const Tensor2<SReal,Int> & T_C_ser;
        const Tensor2<SReal,Int> & T_C_up_ser;

        const Tensor2<SReal,Int> & S_P_ser;
        const Tensor2<SReal,Int> & S_P_v_ser;
        const Tensor2<SReal,Int> & T_P_ser;
        const Tensor2<SReal,Int> & T_P_v_ser;
        
        const typename ClusterTree_T::SparseBinaryMatrix_T & A;
        
        CollisionFinder<ClusterTree_T::AMB_DIM,Real,Int,SReal> C;
        
        const SReal t_init;
        
        SReal a = static_cast<SReal>(0);
        SReal b = static_cast<SReal>(1);
        
        SReal t_max;
        
        Int C_i = -1;
        Int C_j = -1;
        Int P_i = -1;
        Int P_j = -1;
        
    public:
        
        TripleAggregator<Int,Int,SReal,LInt> triples;
        
    public:

        force_inline void LoadClusterS( const Int i )
        {
            C_i = i;
        }
        
        force_inline void LoadClusterT( const Int j )
        {
            C_j = j;
        }
        
        force_inline void LoadPrimitiveS( const Int i )
        {
            P_i = i;
            
        }
        force_inline void LoadPrimitiveT( const Int j )
        {
            P_j = j;
        }
        
        force_inline bool AdmissableQ()
        {
            Compute_AABB_CollisionTimeInterval<ClusterTree_T::AMB_DIM, SReal, Int> (
                S_C_ser.data(C_i), S_C_up_ser.data(C_i),
                T_C_ser.data(C_j), T_C_up_ser.data(C_j),
                a, b
            );
            
            return a >= static_cast<SReal>(1);
        }
        
        force_inline SReal ClusterScoreS()
        {
            return S_C_ser(C_i,0);
        }
        force_inline SReal ClusterScoreT()
        {
            return T_C_ser(C_j,0);
        }
        
        force_inline void ComputeLeafDiagonal()
        {}
        
        force_inline void ComputeLeaf()
        {
            if( A.FindNonzeroPosition(P_i,P_j).found )
            {
                const SReal t = C.FindMaximumSafeStepSize(
                    S_P_ser.data(P_i), S_P_v_ser.data(P_i),
                    T_P_ser.data(P_j), T_P_v_ser.data(P_j),
                    t_init, false
                );
                
                if( t < t_init )
                {
                    triples.Push( P_i, P_j, t );
                    
                    t_max = std::min( t_max, t );
                }
            }
        }
        
        force_inline void ComputeLeafSwapped()
        {
            if( A.FindNonzeroPosition(P_i,P_j).found )
            {
                const SReal t = C.FindMaximumSafeStepSize(
                    S_P_ser.data(P_i), S_P_v_ser.data(P_i),
                    T_P_ser.data(P_j), T_P_v_ser.data(P_j),
                    t_init, false
                );
                
                if( t < t_init )
                {
                    triples.Push( P_j, P_i, t );
                    
                    t_max = std::min( t_max, t );
                }
            }
        }
        
        force_inline void ComputeAdmissable()
        {}
        
        force_inline void ComputeAdmissableSwapped()
        {}
        

    public:
      
        std::string ClassName() const
        {
            return "Collision_Kernel<"+this->tree_string+">";
        }
        
    }; // class Collision_Kernel
    
} //namespace Repulsor
