#pragma once

#define CLASS FMM_Kernel_VF

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
        
        using S_Tree_T = SimplexHierarchy<S_DOM_DIM,AMB_DIM,GJK_Real,Int,SReal>;
        using T_Tree_T = SimplexHierarchy<T_DOM_DIM,AMB_DIM,GJK_Real,Int,SReal>;
        
        using GJK_T    = GJK_Algorithm<AMB_DIM,GJK_Real,Int>;
        
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

        using BASE::bct;
        
        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
        using BASE::metric_data;
        
        mutable S_Tree_T S_Tree;
        mutable T_Tree_T T_Tree;
        
        const SReal * restrict const lambda = S_Tree.Center();
        const SReal * restrict const mu     = T_Tree.Center();
        
        mutable GJK_T gjk;
        
        const Real theta  = 10;
        const Real theta2 = 100;
        const Int max_level = 30;
        
        mutable Int block_count = 0;
        mutable Int max_level_reached = 0;
        mutable Int bottom_count = 0;
        mutable Int inadmissable_count = 0;
        mutable Int primitive_count = 0;
        mutable Real total_sum = static_cast<Real>(0);

#ifdef REPULSOR__PRINT_REPORTS_FOR_ADAPTIVE_KERNELS
        mutable std::ofstream logfile;
        mutable std::ofstream simplex_file;
        mutable std::ofstream center_file;
        mutable std::ofstream emb_simplex_file;
        mutable std::ofstream emb_center_file;
#endif
        
        using BASE::loadS;
        using BASE::loadT;
        using BASE::compute;
        using BASE::writeBlock;
        using BASE::writeT;
        using BASE::writeS;
        
    public:
        
        CLASS() = default;
        
        CLASS( Configurator_T & conf, const Real theta_, const Int  max_level_ = 20 )
        :   BASE        ( conf                                                              )
        ,   S_data      ( GetS().PrimitiveNearFieldData().data()                            )
        ,   S_D_data    ( GetS().ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   S_diag      ( GetS().VF_Accumulator().data(               omp_get_thread_num()) )
        ,   S_ser       ( GetS().PrimitiveSerialized().data()                               )
        ,   T_data      ( GetT().PrimitiveNearFieldData().data()                            )
        ,   T_D_data    ( GetT().ThreadPrimitiveDNearFieldData().data(omp_get_thread_num()) )
        ,   T_diag      ( GetT().VF_Accumulator().data(               omp_get_thread_num()) )
        ,   T_ser       ( GetT().PrimitiveSerialized().data()                               )
        ,   theta       ( theta_                                                            )
        ,   theta2      ( theta_ * theta_                                                   )
        ,   max_level   ( max_level_                                                        )
        {
            if( GetS().PrimitiveSerialized().Dimension(1) != S_Tree.SimplexPrototype().Size() )
            {
                eprint(className()+" Constructor: GetS().PrimitiveSerialized().Dimension(1) != S_Tree.SimplexPrototype().Size()");
            }
            
            if( GetT().PrimitiveSerialized().Dimension(1) != T_Tree.SimplexPrototype().Size() )
            {
                eprint(className()+" Constructor: GetT().PrimitiveSerialized().Dimension(1) != T_Tree.SimplexPrototype().Size()");
            }
        }
        
        CLASS( const CLASS & other )
        :   BASE        ( other                                                                   )
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
        ,   max_level   ( other.max_level                                                         )
        {}
        
        virtual ~CLASS() override
        {
#ifdef REPULSOR__PRINT_REPORTS_FOR_ADAPTIVE_KERNELS
            logfile
            << "\n"
            << "Report for class                    = " << ClassName() << "\n"
            << "Thread ID                           = " << omp_get_thread_num() << "\n"
            << "Number of primitive pairs processed = " << primitive_count << "\n"
            << "Total energy accumulated            = " << total_sum << "\n"
            << "Number of quadrature points         = " << block_count << "\n"
            << "Subdivision level reached           = " << max_level_reached << "\n"
            ;
            if( bottom_count > 0 )
            {
                logfile << "WARNING: Maximal subdivision level = "+ToString(max_level)+" reached "+ToString(bottom_count)+" times. Expect non-sufficent repulsive behavior. Consider refining the mesh.";
            }
            logfile << std::endl;
#endif
        };

    public:
        
        virtual void LoadS( const Int i_global ) override
        {
            const Real * const X  = &S_data[S_DATA_DIM * i_global];
    
            S_Tree.RequireSimplex(S_ser, i_global);
            
            a = X[0];
        
#ifdef NearField_S_Copy
            copy_buffer( &X[1],             &x_buffer[0], S_COORD_DIM );
            copy_buffer( &X[1+S_COORD_DIM], &P[0]       , PROJ_DIM   );
#else
            x_buffer = &X[1];
            P        = &X[1+T_COORD_DIM];
#endif
            if constexpr ( diff_flag )
            {
                zerofy_buffer( &DX[0], S_DATA_DIM );
            }
            
            loadS( i_global );
        }
        
        virtual void LoadT( const Int j_global ) override
        {
            const Real * const Y  = &T_data[T_DATA_DIM * j_global];
    
            T_Tree.RequireSimplex(T_ser, j_global);
            
            b = Y[0];
        
#ifdef NearField_T_Copy
            copy_buffer( &Y[1],             &y_buffer[0], T_COORD_DIM );
            copy_buffer( &Y[1+T_COORD_DIM], &Q[0]       , PROJ_DIM   );
#else
            y_buffer = &Y[1];
            Q        = &Y[1+T_COORD_DIM];
#endif
            if constexpr ( diff_flag )
            {
                zerofy_buffer( &DY[0], T_DATA_DIM );
            }
            
            loadT( j_global );
        }
        
        virtual void Prefetch( const Int j ) const override
        {
            prefetch_range<T_DATA_DIM,0,0>( &T_data[T_DATA_DIM * j] );
            
            prefetch_range<T_Tree_T::SIZE,0,0>( &T_ser[T_Tree_T::SIZE * j] );
            
            if constexpr ( diff_flag )
            {
                prefetch_range<T_DATA_DIM,1,0>( &T_D_data[T_DATA_DIM * j] );
            }
        }
        
        virtual force_inline Real Compute( const Int k_global ) override
        {
            Real sum = static_cast<Real>(0);
            
            ++primitive_count;
            
            bool from_above = true;
            bool shall_continue = true;
            
            S_Tree.ToChild(0);
            T_Tree.ToChild(0);
        
            while( shall_continue )
            {
                if( from_above )
                {
                    if( S_Tree.Level() >= max_level )
                    {
                        // If at lowest level and inadmissable then we just compute the energy and move up.
                        max_level_reached = max_level;
                        block_count++;
                        bottom_count++;
                        sum += compute( k_global );
                        S_Tree.ToParent();
                        T_Tree.ToParent();
                        from_above = false;
                    }
                    else
                    {
                        // If not at lowest level, then we have to check for admissability.
                        auto & P = S_Tree.SimplexPrototype();
                        auto & Q = T_Tree.SimplexPrototype();
                        
                        const bool admissable = gjk.MultipoleAcceptanceCriterion(P, Q, theta2);
                        
                        if( admissable )
                        {
                            // We compute energy, go to parent, and prepare the next child of the parent.
                            max_level_reached = std::max( max_level_reached, S_Tree.Level() );
                            block_count++;
                            sum += compute( k_global );
                            S_Tree.ToParent();
                            T_Tree.ToParent();
                            from_above = false;
                        }
                        else
                        {
                            // If inadmissable, we go a level deeper.
   
                            S_Tree.ToChild(0);
                            T_Tree.ToChild(0);
                            from_above = true;
                        }
                    }
                }
                else
                {
                    // If we come from below, we have to find the next pair of simplices to visit.
                    
                    Int S_k = S_Tree.FormerChildID();
                    Int T_k = T_Tree.FormerChildID();
                    
//                    print("Coming from "+ToString(S_k)+"-th child of S and "+ToString(T_k)+"-th child of T.");
                    
                    if( T_k < T_Tree.ChildCount()-1 )
                    {
                        S_Tree.ToChild(S_k);
                        T_Tree.ToChild(T_k+1);
                        from_above = true;
                    }
                    else
                    {
                        if( S_k < S_Tree.ChildCount()-1 )
                        {
                            S_Tree.ToChild(S_k+1);
                            T_Tree.ToChild(0);
                            from_above = true;
                        }
                        else
                        {
                            // No further unvisited children. Either move up or break.
                            if( S_Tree.Level() == 0 )
                            {
                                shall_continue = false;
                                from_above = true;
                            }
                            else
                            {
                                S_Tree.ToParent();
                                T_Tree.ToParent();
                                from_above = false;
                            }
                        }
                    }
                    
                } // if( from_above )
                
            } // while( shall_continue )

            total_sum += symmetry_factor * sum;
            
            writeBlock( k_global );
            
            return symmetry_factor * sum;
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
        
        virtual LInt NonzeroCount() const override = 0;
        
        virtual void CreateLogFile() const
        {
            DUMP(ClassName());
            logprint("CreateLogFile");
            
#ifdef REPULSOR__PRINT_REPORTS_FOR_ADAPTIVE_KERNELS
            
            logprint("Really creating log file");
            std:: string s = "./Repulsor__"+ClassName()+"_Report_"+ToString(omp_get_thread_num())+".txt";
            
            DUMP(s);
            
            logfile.open(s, std::ios_base::app);
            
            logfile << "Log file for " << ClassName() << std::endl;
            
            s = "./Repulsor__"+ClassName()+"_Simplices_"+ToString(omp_get_thread_num())+".txt";
            simplex_file.open(s, std::ios_base::app);
            
            s = "./Repulsor__"+ClassName()+"_Centers_"+ToString(omp_get_thread_num())+".txt";
            center_file.open(s, std::ios_base::app);
            
            s = "./Repulsor__"+ClassName()+"_EmbSimplices_"+ToString(omp_get_thread_num())+".txt";
            emb_simplex_file.open(s, std::ios_base::app);
            
            s = "./Repulsor__"+ClassName()+"_EmbCenters_"+ToString(omp_get_thread_num())+".txt";
            emb_center_file.open(s, std::ios_base::app);
            
//            logprint("Writing to log file "+s+".");
#endif
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

