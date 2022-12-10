#pragma once

namespace Repulsor
{
    template<
        typename ClusterTree_T_,
        bool is_symmetric_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_
    >
    class FMM_Kernel_FF : public FMM_Kernel<ClusterTree_T_,is_symmetric_,energy_flag_,diff_flag_,metric_flag_>
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
        
        static constexpr Int S_DATA_DIM  = 1 + AMB_DIM + PROJ_DIM;
        static constexpr Int T_DATA_DIM  = 1 + AMB_DIM + PROJ_DIM;
        
    protected:
        
        
        Real * restrict const metric_data = nullptr;
        
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
        
        const  Real * restrict const S_data   = nullptr;
               Real * restrict const S_D_data = nullptr;
               Real * restrict const S_diag   = nullptr;

        const  Real * restrict const T_data   = nullptr;
               Real * restrict const T_D_data = nullptr;
               Real * restrict const T_diag   = nullptr;
        
        using Base_T::tri_i;
        using Base_T::tri_j;
        using Base_T::lin_k;
        
    public:
        
        FMM_Kernel_FF() = default;
        
        explicit FMM_Kernel_FF( Configurator_T & conf )
        :   Base_T        ( conf                                                           )
        ,   S_data      ( GetS().ClusterFarFieldData().data()                            )
        ,   S_D_data    ( GetS().ThreadClusterDFarFieldData().data(omp_get_thread_num()) )
        ,   S_diag      ( GetS().FF_Accumulator().data(            omp_get_thread_num()) )
        ,   T_data      ( GetT().ClusterFarFieldData().data()                            )
        ,   T_D_data    ( GetT().ThreadClusterDFarFieldData().data(omp_get_thread_num()) )
        ,   T_diag      ( GetT().FF_Accumulator().data(            omp_get_thread_num()) )
        {
            if( GetS().ClusterFarFieldData().Dimension(1) != S_DATA_DIM )
            {
                eprint(ClassName()+" Constructor: GetS().ClusterFarFieldData().Dimension(1) != S_DATA_DIM");
            }
            
            if constexpr ( diff_flag )
            {
                if( GetS().ThreadClusterDFarFieldData().Dimension(2) != S_DATA_DIM )
                {
                    eprint(ClassName()+" Constructor: GetS().ThreadClusterDFarFieldData().Dimension(2) != S_DATA_DIM");
                }
            }
               
            if( GetT().ClusterFarFieldData().Dimension(1) != T_DATA_DIM)
            {
               eprint(ClassName()+" Constructor: GetT().ClusterFarFieldData().Dimension(1) != T_DATA_DIM ");
            }
            
            if constexpr ( diff_flag )
            {
                if( GetT().ThreadClusterDFarFieldData().Dimension(2) != T_DATA_DIM )
                {
                    eprint(ClassName()+" Constructor: GetT().ThreadClusterDFarFieldData().Dimension(2) != T_DATA_DIM ");
                }
            }
        }
        
        FMM_Kernel_FF( FMM_Kernel_FF & other )
        :   Base_T      ( other                                                                )
        ,   metric_data ( other.OffDiag().data()                                               )
        ,   S_data      ( other.S_data                                                         )
        ,   S_D_data    ( other.GetS().ThreadClusterDFarFieldData().data(omp_get_thread_num()) )
        ,   S_diag      ( other.GetS().FF_Accumulator().data(            omp_get_thread_num()) )
        ,   T_data      ( other.T_data                                                         )
        ,   T_D_data    ( other.GetT().ThreadClusterDFarFieldData().data(omp_get_thread_num()) )
        ,   T_diag      ( other.GetT().FF_Accumulator().data(            omp_get_thread_num()) )
        {}
        
        ~FMM_Kernel_FF() = default;
        
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
#ifdef FarField_S_Copy
            copy_buffer<AMB_DIM >( &X[1],         &x[0] );
            copy_buffer<PROJ_DIM>( &X[1+AMB_DIM], &P[0] );
#else
            x = &X[1];
            P = &X[1 + AMB_DIM];
#endif
            
            if constexpr ( diff_flag )
            {
                zerofy_buffer<S_DATA_DIM>( &DX[0] );
            }
        }
        
        force_inline void loadT( const Int j_global )
        {
            const Real * const restrict Y = &T_data[T_DATA_DIM * j_global];
            
            b = Y[0];
#ifdef FarField_T_Copy
            copy_buffer<AMB_DIM >( &Y[1],         &y[0] );
            copy_buffer<PROJ_DIM>( &Y[1+AMB_DIM], &Q[0] );
#else
            y = &Y[1];
            Q = &Y[1 + AMB_DIM];
#endif
            if constexpr ( diff_flag )
            {
                zerofy_buffer<T_DATA_DIM>( &DY[0] );
            }
        }
            
        force_inline void writeS( const Int i_global )
        {
            if constexpr ( diff_flag )
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
            return this->metric_values["FF"];
        }
        
        const Values_T & OffDiag() const
        {
            return this->metric_values["FF"];
        }
        
        Values_T & Diag()
        {
            return this->metric_values["FF_diag"];
        }
        
        const Values_T & Diag() const
        {
            return this->metric_values["FF_diag"];
        }

        
        std::string ClassName() const
        {
            return "FMM_Kernel_FF<"
            + GetS().ClassName() + ","
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(metric_flag) +
            ">";
        }

    };

} // namespace Repulsor

