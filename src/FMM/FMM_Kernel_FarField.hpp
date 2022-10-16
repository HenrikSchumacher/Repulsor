#pragma once

#define CLASS FMM_Kernel_FarField
#define BASE  FMM_Kernel<ClusterTree_T_,is_symmetric_,energy_flag,diff_flag,hess_flag,metric_flag>

namespace Repulsor
{
    template<
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
        using BASE::S;
        using BASE::T;
        using BASE::symmetry_factor;
        using BASE::is_symmetric;
        
        static constexpr Int S_DATA_DIM  = 1 + AMB_DIM + PROJ_DIM;
        static constexpr Int T_DATA_DIM  = 1 + AMB_DIM + PROJ_DIM;
        
    protected:
        
        mutable Real sum = static_cast<Real>(0);
        
        mutable Real a   = static_cast<Real>(0);
        mutable Real b   = static_cast<Real>(0);
        
#ifdef FarField_S_Copy
        mutable Real x [AMB_DIM] = {};
        mutable Real P [PROJ_DIM] = {};
#else
        mutable Real const * restrict x = nullptr;
        mutable Real const * restrict P = nullptr;
#endif
        
#ifdef FarField_T_Copy
        mutable Real y [AMB_DIM] = {};
        mutable Real Q [PROJ_DIM] = {};
#else
        mutable Real const * restrict y = nullptr;
        mutable Real const * restrict Q = nullptr;
#endif
        
        mutable Real DX [S_DATA_DIM] = {};
        mutable Real DY [T_DATA_DIM] = {};
        
        const  Real * restrict const S_data     = nullptr;
               Real * restrict const S_D_data   = nullptr;

        const  Real * restrict const T_data     = nullptr;
               Real * restrict const T_D_data   = nullptr;

        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;

        using BASE::S_ID;
        using BASE::T_ID;
                
    public:
        
        CLASS() = default;
        
        CLASS( const ClusterTree_T & S_, const ClusterTree_T & T_ )
        :   BASE( S_, T_ )
        ,   S_data      ( S.ClusterFarFieldData().data()                            )
        ,   S_D_data    ( S.ThreadClusterDFarFieldData().data(omp_get_thread_num()) )
        ,   T_data      ( T.ClusterFarFieldData().data()                            )
        ,   T_D_data    ( T.ThreadClusterDFarFieldData().data(omp_get_thread_num()) )
        {
            if( S.ClusterFarFieldData().Dimension(1) != S_DATA_DIM )
            {
                eprint(className()+" Constructor: S.ClusterFarFieldData().Dimension(1) != S_DATA_DIM");
            }
            
            if( S.ThreadClusterDFarFieldData().Dimension(2) != S_DATA_DIM )
            {
                eprint(className()+" Constructor: S.ThreadClusterDFarFieldData().Dimension(2) != S_DATA_DIM");
            }
               
            if( T.ClusterFarFieldData().Dimension(1) != T_DATA_DIM)
            {
               eprint(className()+" Constructor: T.ClusterFarFieldData().Dimension(1) != T_DATA_DIM ");
            }
            
            if( T.ThreadClusterDFarFieldData().Dimension(2) != T_DATA_DIM )
            {
                eprint(className()+" Constructor: T.ThreadClusterDFarFieldData().Dimension(2) != T_DATA_DIM ");
            }
        }
        
        CLASS( const CLASS & other )
        :   BASE( other )
        ,   S_data      ( other.S_data                      )
        ,   S_D_data    ( other.S.ThreadClusterDFarFieldData().data(omp_get_thread_num()) )
        ,   T_data      ( other.T_data                      )
        ,   T_D_data    ( other.T.ThreadClusterDFarFieldData().data(omp_get_thread_num()) )
        {}
        
        virtual ~CLASS() = default;

    public:
        
        virtual void LoadS( const Int i ) override
        {
            S_ID = i;
            const Real * const restrict X = &S_data[S_DATA_DIM * S_ID];

            a = X[0];
#ifdef FarField_S_Copy
            copy_buffer( &X[1],         &x[0], AMB_DIM  );
            copy_buffer( &X[1+AMB_DIM], &P[0], PROJ_DIM );
#else
            x = &X[1];
            P = &X[1 + AMB_DIM];
#endif
            
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
#ifdef FarField_T_Copy
            copy_buffer( &Y[1],         &y[0], AMB_DIM  );
            copy_buffer( &Y[1+AMB_DIM], &Q[0], PROJ_DIM );
#else
            y = &Y[1];
            Q = &Y[1 + AMB_DIM];
#endif
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
        
        virtual Real Compute() override
        {
           return symmetry_factor * compute();
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
            return TO_STD_STRING(CLASS)+"<"
            + S.ClassName() + ","
            + ToString(is_symmetric) + ","
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(hess_flag) + ","
            + ToString(metric_flag) + ","
            ">";
        }

    };

} // namespace Repulsor

#undef BASE
#undef CLASS

