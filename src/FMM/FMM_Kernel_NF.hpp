#pragma once

namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename ClusterTree_T_,
        bool is_symmetric_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_
    >
    class FMM_Kernel_NF : public FMM_Kernel<ClusterTree_T_,is_symmetric_,energy_flag_,diff_flag_,metric_flag_>
    {
    private:
        
        using Base_T = FMM_Kernel<ClusterTree_T_,is_symmetric_,energy_flag_,diff_flag_,metric_flag_>;
        
    public:
        
        using ClusterTree_T      = ClusterTree_T_;
        
        using Real               = typename Base_T::Real;
        using SReal              = typename Base_T::SReal;
        using ExtReal            = typename Base_T::ExtReal;
        using Int                = typename Base_T::Int;
        using LInt               = typename Base_T::LInt;
        
        using Configurator_T     = typename Base_T::Configurator_T;
        using Values_T           = typename Base_T::Values_T;
        using ValueContainer_T   = typename Base_T::ValueContainer_T;
        
        using Base_T::AMB_DIM;
        using Base_T::PROJ_DIM;
        using Base_T::symmetry_factor;
        using Base_T::is_symmetric;
        using Base_T::energy_flag;
        using Base_T::diff_flag;
        using Base_T::metric_flag;
        using Base_T::GetS;
        using Base_T::GetT;
        
        static constexpr Int S_DOM_DIM = S_DOM_DIM_;
        static constexpr Int T_DOM_DIM = T_DOM_DIM_;
        
        static constexpr Int S_COORD_DIM = (S_DOM_DIM+1)*AMB_DIM;
        static constexpr Int T_COORD_DIM = (T_DOM_DIM+1)*AMB_DIM;
        static constexpr Int S_DATA_DIM  = 1 + S_COORD_DIM + PROJ_DIM;
        static constexpr Int T_DATA_DIM  = 1 + T_COORD_DIM + PROJ_DIM;
        
    protected:
        
        Real * restrict const metric_data = nullptr;
        
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
        
        using Base_T::tri_i;
        using Base_T::tri_j;
        using Base_T::lin_k;
        using Base_T::thread;
    
        static const constexpr Real S_scale = static_cast<Real>(1)/static_cast<Real>(S_DOM_DIM+1);
        static const constexpr Real T_scale = static_cast<Real>(1)/static_cast<Real>(T_DOM_DIM+1);
        
    public:
        
        FMM_Kernel_NF() = default;
        
        // To be used for configuration of kernel.
        explicit FMM_Kernel_NF( Configurator_T & conf, const Int thread_ )
        :   Base_T      ( conf, thread_                                       )
        ,   S_data      ( GetS().PrimitiveNearFieldData().data()              )
        ,   S_D_data    ( GetS().ThreadPrimitiveDNearFieldData().data(thread) )
        ,   S_diag      ( GetS().NF_Accumulator().data(               thread) )
        ,   T_data      ( GetT().PrimitiveNearFieldData().data()              )
        ,   T_D_data    ( GetT().ThreadPrimitiveDNearFieldData().data(thread) )
        ,   T_diag      ( GetT().NF_Accumulator().data(               thread) )
        {
            if( GetS().PrimitiveNearFieldData().Dimension(1) != S_DATA_DIM )
            {
                eprint(ClassName()+" Constructor: GetS().PrimitiveNearFieldData().Dimension(1) != S_DATA_DIM");
            }
            
            if constexpr ( diff_flag )
            {
                if( GetS().ThreadPrimitiveDNearFieldData().Dimension(2) != S_DATA_DIM )
                {
                    eprint(ClassName()+" Constructor: GetS().ThreadPrimitiveDNearFieldData().Dimension(2) != S_DATA_DIM");
                }
            }
               
            if( GetT().PrimitiveNearFieldData().Dimension(1) != T_DATA_DIM)
            {
               eprint(ClassName()+" Constructor: GetT().PrimitiveNearFieldData().Dimension(1) != T_DATA_DIM ");
            }
            
            if constexpr ( diff_flag )
            {
                if( GetT().ThreadPrimitiveDNearFieldData().Dimension(2) != T_DATA_DIM )
                {
                    eprint(ClassName()+" Constructor: GetT().ThreadPrimitiveDNearFieldData().Dimension(2) != T_DATA_DIM ");
                }
            }
        }
        
        FMM_Kernel_NF( FMM_Kernel_NF & other, const Int thread_ )
        :   Base_T      ( other, thread_                                            )
        ,   metric_data ( other.OffDiag().data()                                    )
        ,   S_data      ( other.S_data                                              )
        ,   S_D_data    ( other.GetS().ThreadPrimitiveDNearFieldData().data(thread) )
        ,   S_diag      ( other.GetS().NF_Accumulator().data(               thread) )
        ,   T_data      ( other.T_data                                              )
        ,   T_D_data    ( other.GetT().ThreadPrimitiveDNearFieldData().data(thread) )
        ,   T_diag      ( other.GetT().NF_Accumulator().data(               thread) )
        {}
        
        ~FMM_Kernel_NF() = default;

    public:
        
        force_inline void Prefetch( const Int j_next ) const
        {
            prefetch_buffer<T_DATA_DIM,0,0>( &T_data[T_DATA_DIM * j_next] );
            
            if constexpr ( diff_flag )
            {
                prefetch_buffer<T_DATA_DIM,1,0>( &T_D_data[T_DATA_DIM * j_next] );
            }
        }
        
    protected:
        
        force_inline void loadS( const Int i_global )
        {
            const Real * const restrict X = &S_data[S_DATA_DIM * i_global];
            
            a = X[0];
            
#ifdef NearField_S_Copy
            copy_buffer<PROJ_DIM>( &X[1+S_COORD_DIM], &P[0] );
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
                zerofy_buffer<S_DATA_DIM>( &DX[0] );
            }
        }
        
        force_inline void loadT( const Int j_global )
        {
            const Real * const restrict Y = &T_data[T_DATA_DIM * j_global];
            
            b = Y[0];
            
#ifdef NearField_T_Copy
            copy_buffer<PROJ_DIM>( &Y[1+T_COORD_DIM], &Q[0] );
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
                zerofy_buffer<T_DATA_DIM>( &DY[0] );
            }
        }

        force_inline void writeS( const Int i_global )
        {
            if constexpr (diff_flag )
            {
                Real * restrict const to = &S_D_data[S_DATA_DIM * i_global];
                
                for( Int k = 0; k < S_DATA_DIM; ++k )
                {
                    to[k] += symmetry_factor * DX[k];
                }
            }
        }
        
        force_inline void writeT( const Int j_global )
        {
            if constexpr (diff_flag )
            {
                Real * restrict const to = &T_D_data[T_DATA_DIM * j_global];
                
                for( Int k = 0; k < T_DATA_DIM; ++k )
                {
                    to[k] += symmetry_factor * DY[k];
                }
            }
        }
        
        
    public:
        
        Values_T & OffDiag()
        {
            return this->metric_values["NF"];
        }
        
        const Values_T & OffDiag() const
        {
            return this->metric_values["NF"];
        }
        
        Values_T & Diag()
        {
            return this->metric_values["NF_diag"];
        }
        
        const Values_T & Diag() const
        {
            return this->metric_values["NF_diag"];
        }

        
        void Reduce( const FMM_Kernel_NF & ker )
        {}
        
        std::string ClassName() const
        {
            return "FMM_Kernel_NF<"
            + ToString(S_DOM_DIM) + ","
            + ToString(T_DOM_DIM) + ","
            + GetS().ClassName() + ","
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(metric_flag) + ","
            ">";
        }

    };

} // namespace Repulsor
