#pragma once

#define CLASS FMM_Kernel_NearField
#define BASE  FMM_Kernel<ClusterTree_T_,is_symmetric_,energy_flag,diff_flag,hess_flag,metric_flag>

namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename ClusterTree_T_,
        bool is_symmetric_,
        bool energy_flag, bool diff_flag, bool hess_flag, bool metric_flag
    >
    class CLASS : public BASE
    {
    public:
        
        using ClusterTree_T = ClusterTree_T_;
        
        using Real    = typename ClusterTree_T::Real;
        using Int     = typename ClusterTree_T::Int;
        using SReal   = typename ClusterTree_T::SReal;
        using ExtReal = typename ClusterTree_T::ExtReal;
        
        using BASE::AMB_DIM;
        using BASE::PROJ_DIM;
        using BASE::S_Tree;
        using BASE::T_Tree;
        using BASE::symmetry_factor;
        

        static constexpr Int S_DOM_DIM = S_DOM_DIM_;
        static constexpr Int T_DOM_DIM = T_DOM_DIM_;
        
        static constexpr Int S_COORD_DIM = (S_DOM_DIM+1)*AMB_DIM;
        static constexpr Int T_COORD_DIM = (T_DOM_DIM+1)*AMB_DIM;
        static constexpr Int S_DATA_DIM  = 1 + S_COORD_DIM + PROJ_DIM;
        static constexpr Int T_DATA_DIM  = 1 + T_COORD_DIM + PROJ_DIM;
        
        using S_TREE_T = SimplexHierarchy<S_DOM_DIM,AMB_DIM,GJK_Real,Int,SReal>;
        using T_TREE_T = SimplexHierarchy<T_DOM_DIM,AMB_DIM,GJK_Real,Int,SReal>;
        
        using GJK_T    = GJK_Algorithm<AMB_DIM,GJK_Real,Int>;
        
    protected:
        
        mutable Real sum = static_cast<Real>(0);
        
        mutable Real a   = static_cast<Real>(0);
        mutable Real b   = static_cast<Real>(0);
        
#ifdef NearField_S_Copy
        mutable Real x [AMB_DIM] = {};
        mutable Real p [PROJ_DIM] = {};
        mutable Real x_buffer [S_DATA_DIM] = {};
#else
        mutable Real x [AMB_DIM] = {};
        mutable Real const * restrict p         = nullptr;
        mutable Real const * restrict x_buffer  = nullptr;
#endif
        
#ifdef NearField_T_Copy
        mutable Real y [AMB_DIM] = {};
        mutable Real q [PROJ_DIM] = {};
        mutable Real y_buffer [T_DATA_DIM] = {};
#else
        mutable Real y [AMB_DIM] = {};
        mutable Real const * restrict q         = nullptr;
        mutable Real const * restrict y_buffer  = nullptr;
#endif
        
        mutable Real DX [S_DATA_DIM] = {};
        mutable Real DY [T_DATA_DIM] = {};
        
        const  Real * restrict const S_data     = nullptr;
               Real * restrict const S_D_data   = nullptr;
        const SReal * restrict const S_ser      = nullptr;

        
        const  Real * restrict const T_data     = nullptr;
               Real * restrict const T_D_data   = nullptr;
        const SReal * restrict const T_ser      = nullptr;

        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
        using BASE::S_ID;
        using BASE::T_ID;
        
        static const constexpr Real S_scale = static_cast<Real>(1)/static_cast<Real>(S_DOM_DIM+1);
        static const constexpr Real T_scale = static_cast<Real>(1)/static_cast<Real>(T_DOM_DIM+1);

        mutable S_TREE_T S;
        mutable T_TREE_T T;
        
        const SReal * restrict const lambda = S.Center();
        const SReal * restrict const mu     = T.Center();
        
        mutable GJK_T gjk;
        
    public:
        
        CLASS() = default;
        
        CLASS( const ClusterTree_T & S_, const ClusterTree_T & T_ )
        :   BASE( S_, T_ )
        ,   S_data      ( S_Tree.PrimitiveNearFieldData().data()                            )
        ,   S_D_data    ( S_Tree.ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   S_ser       ( S_Tree.PrimitiveSerialized().data()                               )
        ,   T_data      ( T_Tree.PrimitiveNearFieldData().data()                            )
        ,   T_D_data    ( T_Tree.ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   T_ser       ( T_Tree.PrimitiveSerialized().data()                               )
        {
            if( S_Tree.PrimitiveNearFieldData().Dimension(1) != S_DATA_DIM )
            {
                eprint(className()+" Constructor: S_Tree.PrimitiveNearFieldData().Dimension(1) != S_DATA_DIM");
            }
            
            if( S_Tree.ThreadPrimitiveDNearFieldData().Dimension(2) != S_DATA_DIM )
            {
                eprint(className()+" Constructor: S_Tree.ThreadPrimitiveDNearFieldData().Dimension(2) != S_DATA_DIM");
            }
            
            if( S_Tree.PrimitiveSerialized().Dimension(1) != S.SimplexPrototype().Size() )
            {
                eprint(className()+" Constructor: S_Tree.PrimitiveSerialized().Dimension(1) != S.SimplexPrototype().Size()");
            }
               
            if( T_Tree.PrimitiveNearFieldData().Dimension(1) != T_DATA_DIM)
            {
               eprint(className()+" Constructor: T_Tree.PrimitiveNearFieldData().Dimension(1) != T_DATA_DIM ");
            }
            
            if( T_Tree.ThreadPrimitiveDNearFieldData().Dimension(2) != T_DATA_DIM )
            {
                eprint(className()+" Constructor: T_Tree.ThreadPrimitiveDNearFieldData().Dimension(2) != T_DATA_DIM ");
            }
            
            if( T_Tree.PrimitiveSerialized().Dimension(1) != T.SimplexPrototype().Size() )
            {
                eprint(className()+" Constructor: T_Tree.PrimitiveSerialized().Dimension(1) != T.SimplexPrototype().Size()");
            }
        }
        
        CLASS( const CLASS & other )
        :   BASE( other )
        ,   S_data      ( other.S_data                      )
        ,   S_D_data    ( other.S_Tree.ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   S_ser       ( other.S_ser                       )
        ,   T_data      ( other.T_data                      )
        ,   T_D_data    ( other.T_Tree.ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   T_ser       ( other.T_ser                       )
        {}
        
        virtual ~CLASS() = default;

    public:
        
        virtual void LoadS( const Int i ) override
        {
            S_ID = i;
            const Real * const restrict X = &S_data[S_DATA_DIM * S_ID];
            
            a = X[0];
            
#ifdef NearField_S_Copy
            copy_buffer( &X[1+S_COORD_DIM], &p[0], PROJ_DIM );
#else
            p = &X[1+S_COORD_DIM];
#endif
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                x[k] = S_scale * X[1 + AMB_DIM * 0 + k];
                
                for( Int kk = 1; kk < S_DOM_DIM+1; ++kk )
                {
                    x[k] += S_scale * X[1 + AMB_DIM * kk + k];
                }
            }
            
            if constexpr ( diff_flag )
            {
                zerofy_buffer( &DX[0], S_DATA_DIM );
            }
        }
        
        virtual void LoadT( const Int j ) override
        {
            T_ID = j;
            const Real * const restrict Y = &T_data[T_DATA_DIM * T_ID];
            
            b = Y[0];
            
#ifdef NearField_S_Copy
            copy_buffer( &Y[1+T_COORD_DIM], &q[0], PROJ_DIM );
#else
            q = &Y[1+S_COORD_DIM];
#endif
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                y[k] = T_scale * Y[1 + AMB_DIM * 0 + k];
                
                for( Int kk = 1; kk < T_DOM_DIM+1; ++kk )
                {
                    y[k] += T_scale * Y[1 + AMB_DIM * kk + k];
                }
            }
            
            if constexpr ( diff_flag )
            {
                zerofy_buffer( &DY[0], T_DATA_DIM );
            }
        }

//        void PrefetchS( const Int i ) const
//        {}
        
        virtual void PrefetchT( const Int j ) const override
        {
            prefetch_range<T_DATA_DIM,0,0>( &T_data[T_DATA_DIM * j] );
            
            if( diff_flag )
            {
                prefetch_range<T_DATA_DIM,1,0>( &T_D_data[T_DATA_DIM * j] );
            }
        }
        
        virtual Real compute() override = 0;
            
        virtual void WriteS() override
        {
            if constexpr (diff_flag )
            {
                Real * restrict const to = &S_D_data[S_DATA_DIM * S_ID];
                
                for( Int k = 0; k < S_DATA_DIM; ++k )
                {
                    to[k] += symmetry_factor * DX[k];
                }
            }
        }
        
        virtual void WriteT() override
        {
            if constexpr (diff_flag )
            {
                Real * restrict const to = &T_D_data[T_DATA_DIM * T_ID];
                
                for( Int k = 0; k < T_DATA_DIM; ++k )
                {
                    to[k] += symmetry_factor * DY[k];
                }
            }
        }
        
    public:
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        
    private:
        
        std::string className() const
        {
            return TO_STD_STRING(CLASS)+"<"+S_Tree.ClassName()+">";
        }

    };

} // namespace Repulsor

#undef BASE
#undef CLASS
