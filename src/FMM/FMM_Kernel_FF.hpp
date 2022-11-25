#pragma once

#define CLASS FMM_Kernel_FF
#define BASE  FMM_Kernel<ClusterTree_T_,is_symmetric_,energy_flag_,diff_flag_,metric_flag_>

namespace Repulsor
{
    template<
        typename ClusterTree_T_,
        bool is_symmetric_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_
    >
    class CLASS : public BASE
    {
    public:
        
        using ClusterTree_T      = ClusterTree_T_;
        
        using Real               = typename BASE::Real;
        using SReal              = typename BASE::SReal;
        using ExtReal            = typename BASE::ExtReal;
        using Int                = typename BASE::Int;
        using LInt               = typename BASE::LInt;
        
        using Configurator_T     = typename BASE::Configurator_T;
        using Values_T           = typename BASE::Values_T;
        using ValueContainer_T   = typename BASE::ValueContainer_T;
        
        using BASE::AMB_DIM;
        using BASE::PROJ_DIM;
        using BASE::symmetry_factor;
        using BASE::is_symmetric;
        using BASE::energy_flag;
        using BASE::diff_flag;
        using BASE::metric_flag;
        using BASE::GetS;
        using BASE::GetT;
        
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
        
        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
    public:
        
        CLASS() = default;
        
        explicit CLASS( Configurator_T & conf )
        :   BASE        ( conf                                                           )
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
        
        CLASS( const CLASS & other )
        :   BASE        ( other                                                                )
        ,   metric_data ( other.OffDiag().data()                                               )
        ,   S_data      ( other.S_data                                                         )
        ,   S_D_data    ( other.GetS().ThreadClusterDFarFieldData().data(omp_get_thread_num()) )
        ,   S_diag      ( other.GetS().FF_Accumulator().data(            omp_get_thread_num()) )
        ,   T_data      ( other.T_data                                                         )
        ,   T_D_data    ( other.GetT().ThreadClusterDFarFieldData().data(omp_get_thread_num()) )
        ,   T_diag      ( other.GetT().FF_Accumulator().data(            omp_get_thread_num()) )
        {}
        
        ~CLASS() = default;
        
    public:
        
        force_inline void Prefetch( const Int j_next ) const
        {
            prefetch_range<T_DATA_DIM,0,0>( &T_data[T_DATA_DIM * j_next] );
            
            if constexpr ( diff_flag )
            {
                prefetch_range<T_DATA_DIM,1,0>( &T_D_data[T_DATA_DIM * j_next] );
            }
        }

    protected:
        
        force_inline void loadS( const Int i_global )
        {
            const Real * const restrict X = &S_data[S_DATA_DIM * i_global];

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
        
        force_inline void loadT( const Int j_global )
        {
            const Real * const restrict Y = &T_data[T_DATA_DIM * j_global];
            
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

        Values_T & OffDiag() const
        {
            return this->metric_values["FF"];
        }
        
        Values_T & Diag() const
        {
            return this->metric_values["FF_diag"];
        }

        
        std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"
            + GetS().ClassName() + ","
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(metric_flag) +
            ">";
        }

    };

} // namespace Repulsor

#undef BASE
#undef CLASS

