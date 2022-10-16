#pragma once

namespace Repulsor
{
#define CLASS Collision_Kernel
#define BASE  Traversor_Kernel<ClusterTree_T_>
    
    template<typename ClusterTree_T_>
    class CLASS : public BASE
    {
    public:
        
        using ClusterTree_T = ClusterTree_T_;
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
            const SReal t_init_
        )
        :   BASE( S, T )
        ,   S_C_ser     ( S.ClusterSerializedData()              )
        ,   S_C_up_ser  ( S.ClusterUpdatedSerializedData()       )
        ,   T_C_ser     ( T.ClusterSerializedData()              )
        ,   T_C_up_ser  ( T.ClusterUpdatedSerializedData()       )
        ,   S_P_ser     ( S.PrimitiveSerializedData()            )
        ,   S_P_v_ser   ( S.PrimitiveVelocitiesSerializedData()  )
        ,   T_P_ser     ( T.PrimitiveSerializedData()            )
        ,   T_P_v_ser   ( T.PrimitiveVelocitiesSerializedData()  )
        ,   A           ( S.PrimitiveAdjacencyMatrix()           )
        ,   C           ( S.MovingPrimitivePrototype(), T.MovingPrimitivePrototype() )
        ,   t_init      ( t_init_ )
        ,   t_max       ( t_init_ )
        {}
        

        CLASS( const CLASS & other )
        : BASE( other )
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
        ,   collision_i ( other.collision_i )
        ,   collision_j ( other.collision_j )
        ,   collision_t ( other.collision_t )
        {}
        
        friend void swap(CLASS &A, CLASS &B)
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
            swap( A.collision_i,    B.collision_i   );
            swap( A.collision_j,    B.collision_j   );
            swap( A.collision_t,    B.collision_t   );
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
        
    protected:
        
        const Tensor2<SReal,Int> & S_C_ser;
        const Tensor2<SReal,Int> & S_C_up_ser;
        const Tensor2<SReal,Int> & T_C_ser;
        const Tensor2<SReal,Int> & T_C_up_ser;

        const Tensor2<SReal,Int> & S_P_ser;
        const Tensor2<SReal,Int> & S_P_v_ser;
        const Tensor2<SReal,Int> & T_P_ser;
        const Tensor2<SReal,Int> & T_P_v_ser;
        
        const SparseBinaryMatrixCSR<Int> & A;
        
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
        
        std::vector<Int>   collision_i;
        std::vector<Int>   collision_j;
        std::vector<SReal> collision_t;
        
    public:

        virtual force_inline void LoadClusterS( const Int i ) override
        {
            C_i = i;
        }
        
        virtual force_inline void LoadClusterT( const Int j ) override
        {
            C_j = j;
        }
        
        virtual force_inline void LoadPrimitiveS( const Int i ) override
        {
            P_i = i;
            
        }
        virtual force_inline void LoadPrimitiveT( const Int j ) override
        {
            P_j = j;
        }
        
        virtual force_inline bool IsAdmissable() override
        {
            Compute_AABB_CollisionTimeInterval<ClusterTree_T::AMB_DIM, SReal, Int> (
                S_C_ser.data(C_i), S_C_up_ser.data(C_i),
                T_C_ser.data(C_j), T_C_up_ser.data(C_j),
                a, b
            );
            
            return a >= static_cast<SReal>(1);
        }
        
        virtual force_inline SReal ClusterScoreS() override
        {
            return S_C_ser(C_i,0);
        }
        virtual force_inline SReal ClusterScoreT() override
        {
            return T_C_ser(C_j,0);
        }
        
        virtual force_inline void ComputeLeafDiagonal() override
        {}
        
        virtual force_inline void ComputeLeaf() override
        {
            if( A.FindNonzeroPosition(P_i,P_j) < static_cast<int>(0) )
            {
                const SReal t = C.FindMaximumSafeStepSize(
                    S_P_ser.data(P_i), S_P_v_ser.data(P_i),
                    T_P_ser.data(P_j), T_P_v_ser.data(P_j),
                    t_init, false
                );
                
                if( t < t_init )
                {
                    collision_i.push_back(P_i);
                    collision_j.push_back(P_j);
                    collision_t.push_back(t);
                    
                    t_max = std::min( t_max, t );
                }
            }
        }
        
        virtual force_inline void ComputeLeafSwapped() override
        {
            if( A.FindNonzeroPosition(P_i,P_j) < static_cast<int>(0) )
            {
                const SReal t = C.FindMaximumSafeStepSize(
                    S_P_ser.data(P_i), S_P_v_ser.data(P_i),
                    T_P_ser.data(P_j), T_P_v_ser.data(P_j),
                    t_init, false
                );
                
                if( t < t_init )
                {
                    collision_i.push_back(P_j);
                    collision_j.push_back(P_i);
                    collision_t.push_back(t);
                    
                    t_max = std::min( t_max, t );
                }
            }
        }
        
        virtual force_inline void ComputeAdmissable() override
        {}
        
        virtual force_inline void ComputeAdmissableSwapped() override
        {}
        
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
 
#undef BASE
#undef CLASS
    
} //namespace Repulsor
