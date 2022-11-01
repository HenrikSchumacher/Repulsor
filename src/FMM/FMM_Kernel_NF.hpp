#pragma once

#define CLASS FMM_Kernel_NF
#define BASE  FMM_Kernel<BlockClusterTree_T_,energy_flag_,diff_flag_,metric_flag_>

namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename BlockClusterTree_T_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_
    >
    class CLASS : public BASE
    {
    public:
        
        using BlockClusterTree_T = BlockClusterTree_T_;
        
        using ClusterTree_T      = typename BlockClusterTree_T::ClusterTree_T;
        using Values_T           = typename BlockClusterTree_T::Values_T;
        using ValueContainer_T   = typename BlockClusterTree_T::ValueContainer_T;
        
        using Real               = typename BlockClusterTree_T::Real;
        using SReal              = typename BlockClusterTree_T::SReal;
        using ExtReal            = typename BlockClusterTree_T::ExtReal;
        using Int                = typename BlockClusterTree_T::Int;
        using LInt               = typename BlockClusterTree_T::LInt;
        
        using Configurator_T     = FMM_Configurator<BlockClusterTree_T>;
        
        using BASE::AMB_DIM;
        using BASE::PROJ_DIM;
        using BASE::symmetry_factor;
        using BASE::is_symmetric;
        using BASE::energy_flag;
        using BASE::diff_flag;
        using BASE::metric_flag;
        using BASE::GetS;
        using BASE::GetT;
        
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

        using BASE::bct;
        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
        using BASE::metric_data;
        
        static const constexpr Real S_scale = static_cast<Real>(1)/static_cast<Real>(S_DOM_DIM+1);
        static const constexpr Real T_scale = static_cast<Real>(1)/static_cast<Real>(T_DOM_DIM+1);
        
        using BASE::loadS;
        using BASE::loadT;
        using BASE::compute;
        using BASE::writeBlock;
        using BASE::writeT;
        using BASE::writeS;
        
    public:
        
        CLASS() = default;
        
        // To be used for configuration of kernel.
        CLASS( Configurator_T & conf )
        :   BASE        ( conf                                                              )
        ,   S_data      ( GetS().PrimitiveNearFieldData().data()                            )
        ,   S_D_data    ( GetS().ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   S_diag      ( GetS().NF_Accumulator().data(               omp_get_thread_num()) )
        ,   T_data      ( GetT().PrimitiveNearFieldData().data()                            )
        ,   T_D_data    ( GetT().ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   T_diag      ( GetT().NF_Accumulator().data(               omp_get_thread_num()) )
        {
            if( GetS().PrimitiveNearFieldData().Dimension(1) != S_DATA_DIM )
            {
                eprint(className()+" Constructor: GetS().PrimitiveNearFieldData().Dimension(1) != S_DATA_DIM");
            }
            
            if constexpr ( diff_flag )
            {
                if( GetS().ThreadPrimitiveDNearFieldData().Dimension(2) != S_DATA_DIM )
                {
                    eprint(className()+" Constructor: GetS().ThreadPrimitiveDNearFieldData().Dimension(2) != S_DATA_DIM");
                }
            }
            
            if( GetT().PrimitiveNearFieldData().Dimension(1) != T_DATA_DIM)
            {
               eprint(className()+" Constructor: GetT().PrimitiveNearFieldData().Dimension(1) != T_DATA_DIM ");
            }
            
            if constexpr ( diff_flag )
            {
                if( GetT().ThreadPrimitiveDNearFieldData().Dimension(2) != T_DATA_DIM )
                {
                    eprint(className()+" Constructor: GetT().ThreadPrimitiveDNearFieldData().Dimension(2) != T_DATA_DIM ");
                }
            }
        }
        
        // Use this constructor to make a copy of a configured kernel within the local thread.
        CLASS( const CLASS & other )
        :   BASE        ( other                                                                   )
        ,   S_data      ( other.S_data                                                            )
        ,   S_D_data    ( other.GetS().ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   S_diag      ( other.GetS().NF_Accumulator().data(               omp_get_thread_num()) )
        ,   T_data      ( other.T_data                                                            )
        ,   T_D_data    ( other.GetT().ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   T_diag      ( other.GetT().NF_Accumulator().data(               omp_get_thread_num()) )
        {}
        
        virtual ~CLASS() = default;

    public:
        
        virtual void LoadS( const Int i_global ) override
        {
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
            
            loadS( i_global );
        }
        
        virtual void LoadT( const Int j_global ) override
        {
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
            
            loadT( j_global );
        }
        
        virtual void Prefetch( const Int j_next ) const override
        {
            prefetch_range<T_DATA_DIM,0,0>( &T_data[T_DATA_DIM * j_next] );
            
            if constexpr ( diff_flag )
            {
                prefetch_range<T_DATA_DIM,1,0>( &T_D_data[T_DATA_DIM * j_next] );
            }
        }
        
        virtual force_inline Real Compute( const Int k_global ) override
        {
            const Real result = symmetry_factor * compute( k_global );
            
            writeBlock( k_global );
            
            return result;
        }
            
        virtual void WriteS( const Int i_global ) override
        {
            if constexpr (diff_flag )
            {
                Real * restrict const to = &S_D_data[S_DATA_DIM * i_global];
                
                for( Int k = 0; k < S_DATA_DIM; ++k )
                {
                    to[k] += symmetry_factor * DX[k];
                }
            }
            
            writeS( i_global );
        }
        
        virtual void WriteT( const Int j_global ) override
        {
            if constexpr (diff_flag )
            {
                Real * restrict const to = &T_D_data[T_DATA_DIM * j_global];
                
                for( Int k = 0; k < T_DATA_DIM; ++k )
                {
                    to[k] += symmetry_factor * DY[k];
                }
            }
            
            writeT( j_global );
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
            + bct.ClassName() + ","
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(metric_flag) + ","
            ">";
        }

    };

} // namespace Repulsor

#undef BASE
#undef CLASS
