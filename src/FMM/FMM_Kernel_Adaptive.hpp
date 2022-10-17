#pragma once

#define CLASS FMM_Kernel_Adaptive
#define BASE  FMM_Kernel_NearField<S_DOM_DIM_,T_DOM_DIM_,ClusterTree_T_,is_symmetric_,energy_flag,diff_flag,hess_flag,metric_flag>

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
        using BASE::symmetry_factor;
        using BASE::is_symmetric;
        
        using BASE::S_DOM_DIM;
        using BASE::T_DOM_DIM;
        using BASE::S_COORD_DIM;
        using BASE::T_COORD_DIM;
        using BASE::S_DATA_DIM;
        using BASE::T_DATA_DIM;
        
        using S_Tree_T = SimplexHierarchy<S_DOM_DIM,AMB_DIM,GJK_Real,Int,SReal>;
        using T_Tree_T = SimplexHierarchy<T_DOM_DIM,AMB_DIM,GJK_Real,Int,SReal>;
        
        using GJK_T    = GJK_Algorithm<AMB_DIM,GJK_Real,Int>;
        
    protected:
        
        using BASE::S;
        using BASE::S_data;
        using BASE::S_D_data;
        const SReal * restrict const S_ser = nullptr;

        using BASE::T;
        using BASE::T_data;
        using BASE::T_D_data;
        const SReal * restrict const T_ser = nullptr;
        
        using BASE::DX;
        using BASE::DY;
        
        using BASE::a;
        using BASE::x;
        using BASE::P;

#ifdef NearField_S_Copy
        alignas( ALIGNMENT ) mutable Real x_buffer [S_DATA_DIM] = {};
#else
        mutable Real const * restrict x_buffer  = nullptr;
#endif
        
        using BASE::b;
        using BASE::y;
        using BASE::Q;
        
#ifdef NearField_T_Copy
        alignas( ALIGNMENT ) mutable Real y_buffer [T_DATA_DIM] = {};
#else
        mutable Real const * restrict y_buffer  = nullptr;
#endif
        
        using BASE::S_ID;
        using BASE::T_ID;
        
        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
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
        
    public:
        
        CLASS() = default;
        
        CLASS(
            const ClusterTree_T & S_,
            const ClusterTree_T & T_,
            const Real theta_,
            const Int  max_level_ = 20
        )
        :   BASE        ( S_, T_                         )
        ,   S_ser       ( S.PrimitiveSerialized().data() )
        ,   T_ser       ( T.PrimitiveSerialized().data() )
        ,   theta       ( theta_                         )
        ,   theta2      ( theta_ * theta_                )
        ,   max_level   ( max_level_                     )
        {
            
            if( S.PrimitiveSerialized().Dimension(1) != S_Tree.SimplexPrototype().Size() )
            {
                eprint(className()+" Constructor: S.PrimitiveSerialized().Dimension(1) != S_Tree.SimplexPrototype().Size()");
            }
            
            if( T.PrimitiveSerialized().Dimension(1) != T_Tree.SimplexPrototype().Size() )
            {
                eprint(className()+" Constructor: T.PrimitiveSerialized().Dimension(1) != T_Tree.SimplexPrototype().Size()");
            }
        }
        
        CLASS( const CLASS & other )
        :   BASE( other )
        ,   S_ser       ( other.S_ser     )
        ,   T_ser       ( other.T_ser     )
        ,   theta       ( other.theta     )
        ,   theta2      ( other.theta2    )
        ,   max_level   ( other.max_level )
        {}
        
        virtual ~CLASS() override
        {
#ifdef REPULSOR__PRINT_REPORTS_FOR_ADAPTIVE_KERNELS
            logfile
            << "\n"
            << "Report for class                    = " << this->ClassName() << "\n"
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
        
        virtual force_inline void LoadS( const Int i ) override
        {
            S_ID = i;
            const Real * const X  = &S_data[S_DATA_DIM * S_ID];
    
            S_Tree.RequireSimplex(S_ser, S_ID);
            
            a = X[0];
        
#ifdef NearField_S_Copy
            copy_buffer( &X[1],             &x_buffer[0], S_COORD_DIM );
            copy_buffer( &X[1+S_COORD_DIM], &Q[0]       , PROJ_DIM   );
#else
            x_buffer = &X[1];
            P        = &X[1+T_COORD_DIM];
#endif
            if constexpr ( diff_flag )
            {
                zerofy_buffer( &DX[0], S_DATA_DIM );
            }
        }
        
        virtual force_inline void LoadT( const Int j ) override
        {
            T_ID = j;
            const Real * const Y  = &T_data[T_DATA_DIM * T_ID];
    
            T_Tree.RequireSimplex(T_ser, T_ID);
            
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
        }

//        void PrefetchS( const Int i ) const
//        {}
        
        virtual force_inline void PrefetchT( const Int j ) const override
        {
            prefetch_range<T_DATA_DIM,0,0>( &T_data[T_DATA_DIM * j] );
            
            prefetch_range<T_Tree_T::SIZE,0,0>( &T_ser[T_Tree_T::SIZE * j] );
            
            if constexpr ( diff_flag )
            {
                prefetch_range<T_DATA_DIM,1,0>( &T_D_data[T_DATA_DIM * j] );
            }
        }
        
        virtual Real compute() override = 0;
        
        virtual Real Compute() override
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
                        sum += compute();
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
                            sum += compute();
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
            
            return symmetry_factor * sum;
        }
        
    public:
        
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
            + ToString(hess_flag) + ","
            + ToString(metric_flag) + ","
            ">";
        }

    };

} // namespace Repulsor

#undef BASE
#undef CLASS

