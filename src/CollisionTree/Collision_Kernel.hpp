#pragma once

namespace Repulsor
{
    template<typename ClusterTree_T_>
    class Collision_Kernel
    {
    public:
        
        using ClusterTree_T = ClusterTree_T_;
        using Real    = typename ClusterTree_T::Real;
        using Int     = typename ClusterTree_T::Int;
        using SReal   = typename ClusterTree_T::SReal;
        using ExtReal = typename ClusterTree_T::ExtReal;
        using LInt    = typename ClusterTree_T::LInt;
        
    public:
        
        Collision_Kernel() = delete;
        
        ~Collision_Kernel() = default;
        
        Collision_Kernel(
            cref<ClusterTree_T> S,
            cref<ClusterTree_T> T,
            const Int   thread_,
            const SReal t_init_,
            const SReal TOL
        )
        :   thread      ( thread_                            )
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
        

//        Collision_Kernel( cref<Collision_Kernel> other )
//        :   thread      ( other.thread      )
//        ,   S_C_ser     ( other.S_C_ser     )
//        ,   S_C_up_ser  ( other.S_C_up_ser  )
//        ,   T_C_ser     ( other.T_C_ser     )
//        ,   T_C_up_ser  ( other.T_C_up_ser  )
//        ,   S_P_ser     ( other.S_P_ser     )
//        ,   S_P_v_ser   ( other.S_P_v_ser   )
//        ,   T_P_ser     ( other.T_P_ser     )
//        ,   T_P_v_ser   ( other.T_P_v_ser   )
//        ,   A           ( other.A           )
//        ,   C           ( other.C           )
//        ,   t_init      ( other.t_init      )
//        ,   t_max       ( other.t_max       )
//        ,   triples     ( other.triples     )
//        {}
//        
//        friend void swap( Collision_Kernel & X, Collision_Kernel & Y)
//        {
//            // see https://stackoverflow.com/questions/5695548/public-friend-swap-member-function for details
//            using std::swap;
//
//            swap( X.thread,         Y.thread        );
//            swap( X.S_C_ser,        Y.S_C_ser       );
//            swap( X.S_C_up_ser,     Y.S_C_up_ser    );
//            swap( X.T_C_ser,        Y.T_C_ser       );
//            swap( X.T_C_up_ser,     Y.T_C_up_ser    );
//            swap( X.S_P_ser,        Y.S_P_ser       );
//            swap( X.S_P_v_ser,      Y.S_P_v_ser     );
//            swap( X.T_P_ser,        Y.T_P_ser       );
//            swap( X.T_P_v_ser,      Y.T_P_v_ser     );
//            swap( X.A,              Y.A             );
//            swap( X.C,              Y.C             );
//            swap( X.t_init,         Y.t_init        );
//            swap( X.t_max,          Y.t_max         );
//            swap( X.triples,        Y.triples       );
//        }
//
//        // Copy assignment
//        Collision_Kernel & operator=(Collision_Kernel other) // Intentionally no pass-by-reference here!
//        {
//            swap(*this, other);
//
//            return *this;
//        }
                
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
        
    protected:
        
        const Int thread;
        
        cref<Tensor2<SReal,Int>> S_C_ser;
        cref<Tensor2<SReal,Int>> S_C_up_ser;
        cref<Tensor2<SReal,Int>> T_C_ser;
        cref<Tensor2<SReal,Int>> T_C_up_ser;

        cref<Tensor2<SReal,Int>> S_P_ser;
        cref<Tensor2<SReal,Int>> S_P_v_ser;
        cref<Tensor2<SReal,Int>> T_P_ser;
        cref<Tensor2<SReal,Int>> T_P_v_ser;
        
        cref<typename ClusterTree_T::SparseBinaryMatrix_T> A;
        
        CollisionFinder<ClusterTree_T::AMB_DIM,Real,Int,SReal> C;
        
        const SReal t_init;
        
        SReal a = Scalar::Zero<SReal>;
        SReal b = Scalar::One<SReal>;
        
        SReal t_max;
        
        Int C_i = -1;
        Int C_j = -1;
        Int P_i = -1;
        Int P_j = -1;
        
    public:
        
        TripleAggregator<Int,Int,SReal,LInt> triples;
        
    public:

        TOOLS_FORCE_INLINE void LoadClusterS( const Int i )
        {
            C_i = i;
        }
        
        TOOLS_FORCE_INLINE void LoadClusterT( const Int j )
        {
            C_j = j;
        }
        
        TOOLS_FORCE_INLINE void LoadPrimitiveS( const Int i )
        {
            P_i = i;
            
        }
        TOOLS_FORCE_INLINE void LoadPrimitiveT( const Int j )
        {
            P_j = j;
        }
        
        TOOLS_FORCE_INLINE bool AdmissableQ()
        {
            Compute_AABB_CollisionTimeInterval<ClusterTree_T::AMB_DIM, SReal, Int> (
                S_C_ser.data(C_i), S_C_up_ser.data(C_i),
                T_C_ser.data(C_j), T_C_up_ser.data(C_j),
                a, b
            );
            
            return a >= Scalar::One<SReal>;
        }
        
        TOOLS_FORCE_INLINE SReal ClusterScoreS()
        {
            return S_C_ser(C_i,0);
        }
        TOOLS_FORCE_INLINE SReal ClusterScoreT()
        {
            return T_C_ser(C_j,0);
        }
        
        TOOLS_FORCE_INLINE void ComputeLeafDiagonal()
        {}
        
        TOOLS_FORCE_INLINE void ComputeLeaf()
        {
            if( A.NonzeroPositionQ(P_i,P_j) )
            {
                const SReal t = C.FindMaximumSafeStepSize(
                    S_P_ser.data(P_i), S_P_v_ser.data(P_i),
                    T_P_ser.data(P_j), T_P_v_ser.data(P_j),
                    t_init, false
                );
                
                if( t < t_init )
                {
                    triples.Push( P_i, P_j, t );
                    
                    t_max = Min( t_max, t );
                }
            }
        }
        
        TOOLS_FORCE_INLINE void ComputeLeafSwapped()
        {
            if( A.NonzeroPositionQ(P_i,P_j) )
            {
                const SReal t = C.FindMaximumSafeStepSize(
                    S_P_ser.data(P_i), S_P_v_ser.data(P_i),
                    T_P_ser.data(P_j), T_P_v_ser.data(P_j),
                    t_init, false
                );
                
                if( t < t_init )
                {
                    triples.Push( P_j, P_i, t );
                    
                    t_max = Min( t_max, t );
                }
            }
        }
        
        TOOLS_FORCE_INLINE void ComputeAdmissable()
        {}
        
        TOOLS_FORCE_INLINE void ComputeAdmissableSwapped()
        {}
        

    private:
      
        static std::string className()
        {
            return std::string("Collision_Kernel<...>");
        }
        
    public:
        
        std::string ClassName() const
        {
            return className();
        }
        
    }; // class Collision_Kernel
    
} //namespace Repulsor
