#pragma once

namespace Repulsor
{
    template<typename ClusterTree_T_>
    class MaximumSafeStepSize_Kernel
    {
        
    public:
        
        using ClusterTree_T = ClusterTree_T_;
        using Real    = typename ClusterTree_T::Real;
        using Int     = typename ClusterTree_T::Int;
        using SReal   = typename ClusterTree_T::SReal;
        using ExtReal = typename ClusterTree_T::ExtReal;
        
    public:
        
        MaximumSafeStepSize_Kernel() = delete;
        
        ~MaximumSafeStepSize_Kernel() = default;
        
        MaximumSafeStepSize_Kernel(
            cref<ClusterTree_T> S,
            cref<ClusterTree_T> T,
            const Int thread_,
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
        

        MaximumSafeStepSize_Kernel( cref<MaximumSafeStepSize_Kernel> other, const Int thread_ )
        :   thread      ( other.thread      )
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
        {}
        
        friend void swap( mref<MaximumSafeStepSize_Kernel> X, mref<MaximumSafeStepSize_Kernel> Y )
        {
            // see https://stackoverflow.com/questions/5695548/public-friend-swap-member-function for details
            using std::swap;

            swap( X.thread,         Y.thread        );
            swap( X.S_C_ser,        Y.S_C_ser       );
            swap( X.S_C_up_ser,     Y.S_C_up_ser    );
            swap( X.T_C_ser,        Y.T_C_ser       );
            swap( X.T_C_up_ser,     Y.T_C_up_ser    );
            swap( X.S_P_ser,        Y.S_P_ser       );
            swap( X.S_P_v_ser,      Y.S_P_v_ser     );
            swap( X.T_P_ser,        Y.T_P_ser       );
            swap( X.T_P_v_ser,      Y.T_P_v_ser     );
            swap( X.A,              Y.A             );
            swap( X.C,              Y.C             );
            swap( X.t_init,         Y.t_init        );
            swap( X.t_max,          Y.t_max         );
        }

        // Copy assignment
        MaximumSafeStepSize_Kernel & operator=(MaximumSafeStepSize_Kernel other) // Intentionally no pass-by-reference here!
        {
            swap(*this, other);

            return *this;
        }
                
//        // Move constructor
//        MaximumSafeStepSize_Kernel( MaximumSafeStepSize_Kernel && other ) noexcept
//        {
//            swap(*this, other);
//        }
//
//        // Move assignment operator
//        MaximumSafeStepSize_Kernel & operator=( MaximumSafeStepSize_Kernel && other ) noexcept
//        {
//            if( this != &other )
//            {
//                swap( *this, other );
//            }
//            return *this;
//        }
        
    public:
        
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
        
        SReal t_max;
        
        Int C_i = -1;
        Int C_j = -1;
        Int P_i = -1;
        Int P_j = -1;
        
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
            SReal a = Scalar::Zero<SReal>;
            SReal b = Scalar::One<SReal>;
            
            Compute_AABB_CollisionTimeInterval<ClusterTree_T::AMB_DIM, SReal, Int> (
                S_C_ser.data(C_i), S_C_up_ser.data(C_i),
                T_C_ser.data(C_j), T_C_up_ser.data(C_j),
                a, b
            );
            
            return (a * t_init) >= t_max;
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
            if( !A.FindNonzeroPosition(P_i,P_j).found )
            {
                t_max = C.FindMaximumSafeStepSize(
                    S_P_ser.data(P_i), S_P_v_ser.data(P_i),
                    T_P_ser.data(P_j), T_P_v_ser.data(P_j),
                    t_max, false
                );
            }
        }
        
        force_inline void ComputeLeafSwapped()
        {
            if( !A.FindNonzeroPosition(P_i,P_j).found )
            {
                t_max = C.FindMaximumSafeStepSize(
                    S_P_ser.data(P_i), S_P_v_ser.data(P_i),
                    T_P_ser.data(P_j), T_P_v_ser.data(P_j),
                    t_max, false
                );
            }
        }
        
        force_inline void ComputeAdmissable()
        {}
        
        force_inline void ComputeAdmissableSwapped()
        {}

        SReal MaxTime() const
        {
            return t_max;
        }
        
    public:
      
        std::string ClassName() const
        {
            return "MaximumSafeStepSize_Kernel<...>";
        }
        
    }; // class MaximumSafeStepSize_Kernel
 
} //namespace Repulsor


