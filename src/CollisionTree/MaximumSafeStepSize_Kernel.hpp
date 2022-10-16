#pragma once

namespace Repulsor
{
#define CLASS MaximumSafeStepSize_Kernel
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
        ,   S_C_ser     ( S.ClusterSerialized()              )
        ,   S_C_up_ser  ( S.ClusterUpdatedSerialized()       )
        ,   T_C_ser     ( T.ClusterSerialized()              )
        ,   T_C_up_ser  ( T.ClusterUpdatedSerialized()       )
        ,   S_P_ser     ( S.PrimitiveSerialized()            )
        ,   S_P_v_ser   ( S.PrimitiveVelocitiesSerialized()  )
        ,   T_P_ser     ( T.PrimitiveSerialized()            )
        ,   T_P_v_ser   ( T.PrimitiveVelocitiesSerialized()  )
        ,   A           ( S.PrimitiveAdjacencyMatrix()       )
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
        
        SReal t_max;
        
        Int C_i = -1;
        Int C_j = -1;
        Int P_i = -1;
        Int P_j = -1;
        
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
            SReal a = static_cast<SReal>(0);
            SReal b = static_cast<SReal>(1);
            
            Compute_AABB_CollisionTimeInterval<ClusterTree_T::AMB_DIM, SReal, Int> (
                S_C_ser.data(C_i), S_C_up_ser.data(C_i),
                T_C_ser.data(C_j), T_C_up_ser.data(C_j),
                a, b
            );
            
            return (a * t_init) >= t_max;
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
                t_max = C.FindMaximumSafeStepSize(
                    S_P_ser.data(P_i), S_P_v_ser.data(P_i),
                    T_P_ser.data(P_j), T_P_v_ser.data(P_j),
                    t_max, false
                );
            }
        }
        
        virtual force_inline void ComputeLeafSwapped() override
        {
            if( A.FindNonzeroPosition(P_i,P_j) < static_cast<int>(0) )
            {
                t_max = C.FindMaximumSafeStepSize(
                    S_P_ser.data(P_i), S_P_v_ser.data(P_i),
                    T_P_ser.data(P_j), T_P_v_ser.data(P_j),
                    t_max, false
                );
            }
        }
        
        virtual force_inline void ComputeAdmissable() override
        {}
        
        virtual force_inline void ComputeAdmissableSwapped() override
        {}

        force_inline SReal MaxTime() const
        {
            return t_max;
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
 
#undef BASE
#undef CLASS
    
} //namespace Repulsor


