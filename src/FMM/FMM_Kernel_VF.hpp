#pragma once


namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename ClusterTree_T_,
        bool is_symmetric_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_
    >
    class FMM_Kernel_VF : public FMM_Kernel<ClusterTree_T_,is_symmetric_,energy_flag_,diff_flag_,metric_flag_>
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
        
        using S_Tree_T = SimplexHierarchy<S_DOM_DIM,AMB_DIM,GJK_Real,Int,SReal>;
        using T_Tree_T = SimplexHierarchy<T_DOM_DIM,AMB_DIM,GJK_Real,Int,SReal>;
        
        using GJK_T    = GJK_Algorithm<AMB_DIM,GJK_Real,Int>;
        
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
        const SReal * restrict const S_ser     = nullptr;
        
        const  Real * restrict const T_data    = nullptr;
               Real * restrict const T_D_data  = nullptr;
               Real * restrict const T_diag    = nullptr;
        const SReal * restrict const T_ser     = nullptr;

#ifdef NearField_S_Copy
        mutable Real x_buffer [S_DATA_DIM] = {};
#else
        mutable Real const * restrict x_buffer  = nullptr;
#endif
        
#ifdef NearField_T_Copy
        mutable Real y_buffer [T_DATA_DIM] = {};
#else
        mutable Real const * restrict y_buffer  = nullptr;
#endif
        
        using Base_T::tri_i;
        using Base_T::tri_j;
        using Base_T::lin_k;

        
        mutable S_Tree_T S_Tree;
        mutable T_Tree_T T_Tree;
        
        const SReal * restrict const lambda = S_Tree.Center();
        const SReal * restrict const mu     = T_Tree.Center();
        
        mutable GJK_T gjk;
        
        const Real theta  = 10;
        const Real theta2 = 100;
        const Int  max_refinement = 20;
        
        mutable Int block_count = 0;
        mutable Int max_level_reached = 0;
        mutable Int bottom_count = 0;
        mutable Int inadmissable_count = 0;
        mutable Int primitive_count = 0;
        mutable Real total_sum = static_cast<Real>(0);

//#ifdef REPULSOR__PRINT_REPORTS_FOR_ADAPTIVE_KERNELS
//        mutable std::ofstream logfile;
//        mutable std::ofstream simplex_file;
//        mutable std::ofstream center_file;
//        mutable std::ofstream emb_simplex_file;
//        mutable std::ofstream emb_center_file;
//#endif
        
    public:
        
        FMM_Kernel_VF() = default;
        
        FMM_Kernel_VF( Configurator_T & conf, const Real theta_, const Int max_refinement_  )
        :   Base_T         ( conf                                                              )
        ,   S_data         ( GetS().PrimitiveNearFieldData().data()                            )
        ,   S_D_data       ( GetS().ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   S_diag         ( GetS().VF_Accumulator().data(               omp_get_thread_num()) )
        ,   S_ser          ( GetS().PrimitiveSerialized().data()                               )
        ,   T_data         ( GetT().PrimitiveNearFieldData().data()                            )
        ,   T_D_data       ( GetT().ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   T_diag         ( GetT().VF_Accumulator().data(               omp_get_thread_num()) )
        ,   T_ser          ( GetT().PrimitiveSerialized().data()                               )
        ,   theta          ( theta_                                                            )
        ,   theta2         ( theta_ * theta_                                                   )
        ,   max_refinement ( max_refinement_                                                   )
        {
            if( GetS().PrimitiveSerialized().Dimension(1) != S_Tree.SimplexPrototype().Size() )
            {
                eprint(ClassName()+" Constructor: GetS().PrimitiveSerialized().Dimension(1) != S_Tree.SimplexPrototype().Size()");
            }
            
            if( GetT().PrimitiveSerialized().Dimension(1) != T_Tree.SimplexPrototype().Size() )
            {
                eprint(ClassName()+" Constructor: GetT().PrimitiveSerialized().Dimension(1) != T_Tree.SimplexPrototype().Size()");
            }
        }
        
        FMM_Kernel_VF( FMM_Kernel_VF & other )
        :   Base_T        ( other                                                                   )
        ,   metric_data ( other.OffDiag().data()                                                  )
        ,   S_data      ( other.S_data                                                            )
        ,   S_D_data    ( other.GetS().ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   S_diag      ( other.GetS().VF_Accumulator().data(               omp_get_thread_num()) )
        ,   S_ser       ( other.S_ser                                                             )
        ,   T_data      ( other.T_data                                                            )
        ,   T_D_data    ( other.GetT().ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   T_diag      ( other.GetT().VF_Accumulator().data(               omp_get_thread_num()) )
        ,   T_ser       ( other.T_ser                                                             )
        ,   theta       ( other.theta                                                             )
        ,   theta2      ( other.theta2                                                            )
        ,   max_refinement ( other.max_refinement                                               )
        {}
        
        ~FMM_Kernel_VF()
        {
//#ifdef REPULSOR__PRINT_REPORTS_FOR_ADAPTIVE_KERNELS
//            logfile
//            << "\n"
//            << "Report for class                    = " << ClassName() << "\n"
//            << "Thread ID                           = " << omp_get_thread_num() << "\n"
//            << "Number of primitive pairs processed = " << primitive_count << "\n"
//            << "Total energy accumulated            = " << total_sum << "\n"
//            << "Number of quadrature points         = " << block_count << "\n"
//            << "Subdivision level reached           = " << max_level_reached << "\n"
//            ;
//            if( bottom_count > 0 )
//            {
//                logfile << "WARNING: Maximal subdivision level = "+ToString(max_level)+" reached "+ToString(bottom_count)+" times. Expect non-sufficent repulsive behavior. Consider refining the mesh.";
//            }
//            logfile << std::endl;
//#endif
        };

    public:
        
        force_inline void Prefetch( const Int j ) const
        {
            prefetch_buffer<T_DATA_DIM,0,0>( &T_data[T_DATA_DIM * j] );
            
            prefetch_buffer<T_Tree_T::SIZE,0,0>( &T_ser[T_Tree_T::SIZE * j] );
            
            if constexpr ( diff_flag )
            {
                prefetch_buffer<T_DATA_DIM,1,0>( &T_D_data[T_DATA_DIM * j] );
            }
        }
        
    protected:

        force_inline void loadS( const Int i_global )
        {
            const Real * const X = &S_data[S_DATA_DIM * i_global];
    
            S_Tree.RequireSimplex(S_ser, i_global);
            
            a = X[0];
        
#ifdef NearField_S_Copy
            copy_buffer<S_COORD_DIM>( &X[1],             &x_buffer[0] );
            copy_buffer<PROJ_DIM   >( &X[1+S_COORD_DIM], &P[0]        );
#else
            x_buffer = &X[1];
            P        = &X[1+T_COORD_DIM];
#endif
            if constexpr ( diff_flag )
            {
                zerofy_buffer<S_DATA_DIM>( &DX[0] );
            }
        }
        
        force_inline void loadT( const Int j_global )
        {
            const Real * const Y = &T_data[T_DATA_DIM * j_global];
    
            T_Tree.RequireSimplex(T_ser, j_global);
            
            b = Y[0];
        
#ifdef NearField_T_Copy
            copy_buffer<T_COORD_DIM>( &Y[1],             &y_buffer[0] );
            copy_buffer<PROJ_DIM   >( &Y[1+T_COORD_DIM], &Q[0]        );
#else
            y_buffer = &Y[1];
            Q        = &Y[1+T_COORD_DIM];
#endif
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
            return this->metric_values["VF"];
        }
        
        const Values_T & OffDiag() const
        {
            return this->metric_values["VF"];
        }
        
        Values_T & Diag()
        {
            return this->metric_values["VF_diag"];
        }
        
        const Values_T & Diag() const
        {
            return this->metric_values["VF_diag"];
        }
        
//        void CreateLogFile() const
//        {
//            DUMP(ClassName());
//            logprint("CreateLogFile");
//
//#ifdef REPULSOR__PRINT_REPORTS_FOR_ADAPTIVE_KERNELS
//
//            logprint("Really creating log file");
//            std:: string s = "./Repulsor__"+ClassName()+"_Report_"+ToString(omp_get_thread_num())+".txt";
//
//            DUMP(s);
//
//            logfile.open(s, std::ios_base::app);
//
//            logfile << "Log file for " << ClassName() << std::endl;
//
//            s = "./Repulsor__"+ClassName()+"_Simplices_"+ToString(omp_get_thread_num())+".txt";
//            simplex_file.open(s, std::ios_base::app);
//
//            s = "./Repulsor__"+ClassName()+"_Centers_"+ToString(omp_get_thread_num())+".txt";
//            center_file.open(s, std::ios_base::app);
//
//            s = "./Repulsor__"+ClassName()+"_EmbSimplices_"+ToString(omp_get_thread_num())+".txt";
//            emb_simplex_file.open(s, std::ios_base::app);
//
//            s = "./Repulsor__"+ClassName()+"_EmbCenters_"+ToString(omp_get_thread_num())+".txt";
//            emb_center_file.open(s, std::ios_base::app);
//
////            logprint("Writing to log file "+s+".");
//#endif
//        }

        
    public:
        
        Int MaxLevelReached() const
        {
            return max_level_reached;
        }
        
        void Reduce( FMM_Kernel_VF & ker )
        {
            max_level_reached = std::max( max_level_reached, ker.max_level_reached);
        }
        
        std::string ClassName() const
        {
            return "FMM_Kernel_VF<"
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

