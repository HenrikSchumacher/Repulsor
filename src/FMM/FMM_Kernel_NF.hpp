#pragma once

#define CLASS FMM_Kernel_NF
#define BASE  FMM_Kernel<ClusterTree_T_,is_symmetric_,energy_flag_,diff_flag_,metric_flag_>

namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename ClusterTree_T_,
        bool is_symmetric_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_
    >
    class CLASS : public BASE
    {
    public:
        
        using ClusterTree_T = ClusterTree_T_;
        
        using Real    = typename ClusterTree_T::Real;
        using Int     = typename ClusterTree_T::Int;
        using SReal   = typename ClusterTree_T::SReal;
        using ExtReal = typename ClusterTree_T::ExtReal;
        using typename BASE::Configurator_T;
        
        using BASE::AMB_DIM;
        using BASE::PROJ_DIM;
        using BASE::S;
        using BASE::T;
        using BASE::symmetry_factor;
        using BASE::is_symmetric;
        using BASE::energy_flag;
        using BASE::diff_flag;
        using BASE::metric_flag;

        static constexpr Int S_DOM_DIM = S_DOM_DIM_;
        static constexpr Int T_DOM_DIM = T_DOM_DIM_;
        
        static constexpr Int S_COORD_DIM = (S_DOM_DIM+1)*AMB_DIM;
        static constexpr Int T_COORD_DIM = (T_DOM_DIM+1)*AMB_DIM;
        static constexpr Int S_DATA_DIM  = 1 + S_COORD_DIM + PROJ_DIM;
        static constexpr Int T_DATA_DIM  = 1 + T_COORD_DIM + PROJ_DIM;
        
    protected:
        
        mutable Real sum = static_cast<Real>(0);
        
        mutable Real a   = static_cast<Real>(0);
        mutable Real b   = static_cast<Real>(0);
        
        mutable Real x [AMB_DIM] = {};
#ifdef NearField_S_Copy
        mutable Real P [PROJ_DIM] = {};
#else
        mutable Real const * restrict P = nullptr;
#endif
        
        mutable Real y [AMB_DIM] = {};
#ifdef NearField_T_Copy
        mutable Real Q [PROJ_DIM] = {};
#else
        mutable Real const * restrict Q = nullptr;
#endif
        
        mutable Real DX [S_DATA_DIM] = {};
        mutable Real DY [T_DATA_DIM] = {};
        
        const  Real * restrict const S_data    = nullptr;
               Real * restrict const S_D_data  = nullptr;
               Real * restrict const S_diag    = nullptr;
        
        const  Real * restrict const T_data    = nullptr;
               Real * restrict const T_D_data  = nullptr;
               Real * restrict const T_diag    = nullptr;


        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;

        using BASE::i_global;
        using BASE::j_global;
        using BASE::k_global;
        
        using BASE::metric_data;
        
        static const constexpr Real S_scale = static_cast<Real>(1)/static_cast<Real>(S_DOM_DIM+1);
        static const constexpr Real T_scale = static_cast<Real>(1)/static_cast<Real>(T_DOM_DIM+1);
        
        using BASE::loadS;
        using BASE::loadT;
        using BASE::compute;
        using BASE::writeT;
        using BASE::writeS;
        
    public:
        
        CLASS() = default;
        
        // To be used for configuration of kernel.
        CLASS( Configurator_T & conf )
        :   BASE        ( conf                                                         )
        ,   S_data      ( S.PrimitiveNearFieldData().data()                            )
        ,   S_D_data    ( S.ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   S_diag      ( S.NF_Accumulator().data(               omp_get_thread_num()) )
        ,   T_data      ( T.PrimitiveNearFieldData().data()                            )
        ,   T_D_data    ( T.ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   T_diag      ( T.NF_Accumulator().data(               omp_get_thread_num()) )
        {
            if( S.PrimitiveNearFieldData().Dimension(1) != S_DATA_DIM )
            {
                eprint(className()+" Constructor: S.PrimitiveNearFieldData().Dimension(1) != S_DATA_DIM");
            }
            
            if constexpr ( diff_flag )
            {
                if( S.ThreadPrimitiveDNearFieldData().Dimension(2) != S_DATA_DIM )
                {
                    eprint(className()+" Constructor: S.ThreadPrimitiveDNearFieldData().Dimension(2) != S_DATA_DIM");
                }
            }
            
            if( T.PrimitiveNearFieldData().Dimension(1) != T_DATA_DIM)
            {
               eprint(className()+" Constructor: T.PrimitiveNearFieldData().Dimension(1) != T_DATA_DIM ");
            }
            
            if constexpr ( diff_flag )
            {
                if( T.ThreadPrimitiveDNearFieldData().Dimension(2) != T_DATA_DIM )
                {
                    eprint(className()+" Constructor: T.ThreadPrimitiveDNearFieldData().Dimension(2) != T_DATA_DIM ");
                }
            }
        }
        
        // Use this constructor to make a copy of a configured kernel within the local thread.
        CLASS( const CLASS & other )
        :   BASE        ( other                                                              )
        ,   S_data      ( other.S_data                                                       )
        ,   S_D_data    ( other.S.ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   S_diag      ( other.S.NF_Accumulator().data(               omp_get_thread_num()) )
        ,   T_data      ( other.T_data                                                       )
        ,   T_D_data    ( other.T.ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   T_diag      ( other.T.NF_Accumulator().data(               omp_get_thread_num()) )
        {}
        
        virtual ~CLASS() = default;

    public:
        
        virtual void LoadS( const Int i_global_ ) override
        {
            i_global = i_global_;
            const Real * const restrict X = &S_data[S_DATA_DIM * i_global];
            
            a = X[0];
            
#ifdef NearField_S_Copy
            copy_buffer( &X[1+S_COORD_DIM], &P[0], PROJ_DIM );
#else
            P = &X[1+S_COORD_DIM];
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
            
            loadS();
        }
        
        virtual void LoadT( const Int j_global_ ) override
        {
            j_global = j_global_;
            const Real * const restrict Y = &T_data[T_DATA_DIM * j_global];
            
            b = Y[0];
            
#ifdef NearField_T_Copy
            copy_buffer( &Y[1+T_COORD_DIM], &Q[0], PROJ_DIM );
#else
            Q = &Y[1+S_COORD_DIM];
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
            
            loadT();
        }
        
        virtual void Prefetch( const Int j_next ) const override
        {
            prefetch_range<T_DATA_DIM,0,0>( &T_data[T_DATA_DIM * j_next] );
            
            if constexpr ( diff_flag )
            {
                prefetch_range<T_DATA_DIM,1,0>( &T_D_data[T_DATA_DIM * j_next] );
            }
        }
        
        virtual force_inline Real Compute( const Int k_global_ ) override
        {
            k_global = k_global_;
            return symmetry_factor * compute();
        }
            
        virtual void WriteS() override
        {
            if constexpr (diff_flag )
            {
                Real * restrict const to = &S_D_data[S_DATA_DIM * i_global];
                
                for( Int k = 0; k < S_DATA_DIM; ++k )
                {
                    to[k] += symmetry_factor * DX[k];
                }
            }
            
            writeS();
        }
        
        virtual void WriteT() override
        {
            if constexpr (diff_flag )
            {
                Real * restrict const to = &T_D_data[T_DATA_DIM * j_global];
                
                for( Int k = 0; k < T_DATA_DIM; ++k )
                {
                    to[k] += symmetry_factor * DY[k];
                }
            }
            
            writeT();
        }
        
    protected:

        
    public:
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        
    private:
        
        std::string className() const
        {
            return TO_STD_STRING(CLASS)+"<"
            + ToString(S_DOM_DIM) + ","
            + ToString(T_DOM_DIM) + ","
            + S.ClassName() + ","
            + ToString(is_symmetric) + ","
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(metric_flag) + ","
            ">";
        }

    };

} // namespace Repulsor

#undef BASE
#undef CLASS
