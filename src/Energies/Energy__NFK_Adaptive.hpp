#pragma once

#define CLASS Energy__NFK_Adaptive
#define BASE  Energy__NFK<DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Int,SReal>

#define ADMISSABILITY_CONDITION const bool admissable = gjk.MultipoleAcceptanceCriterion(P, Q, theta2);

namespace Repulsor
{
    template<int DOM_DIM1, int DOM_DIM2, int AMB_DIM, typename Real, typename Int, typename SReal>
    class CLASS : public BASE
    {
    public:
        
    protected:

        using BASE::S;
        using BASE::S_serialized;
        using BASE::S_near;
        using BASE::S_D_near;
        
        using BASE::T;
        using BASE::T_serialized;
        using BASE::T_near;
        using BASE::T_D_near;
        
        using BASE::DX;
        using BASE::DY;
        
        using BASE::a;
        using BASE::x;
        using BASE::p;
        using BASE::x_buffer;
        
        using BASE::b;
        using BASE::y;
        using BASE::q;
        using BASE::y_buffer;
        
        using BASE::S_ID;
        using BASE::T_ID;
        
        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
        using BASE::gjk;
        
    public:
        
        using BASE::NEAR_DIMS;
        using BASE::NEAR_DIMT;
        using BASE::COORD_DIMS;
        using BASE::COORD_DIMT;
        using BASE::PROJ_DIM;
        
    public:
        
        explicit CLASS( const AdaptivitySettings & settings_ = AdaptivitySettings() )
        : BASE()
        , settings(settings_)
        , theta(settings_.theta)
        , theta2(settings_.theta * settings_.theta)
        , intersection_theta2( settings_.intersection_theta * settings_.intersection_theta)
        {}

        // Copy constructor
        CLASS( const CLASS & other )
        : BASE(other)
        , settings      ( other.settings      )
        , theta         ( other.settings.theta )
        , theta2        ( other.settings.theta * other.settings.theta )
        , intersection_theta2( other.settings.intersection_theta * other.settings.intersection_theta)
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
                logfile << "WARNING: Maximal subdivision level = "+ToString(settings.max_level)+" reached "+ToString(bottom_count)+" times. Expect non-sufficent repulsive behavior. Consider refining the mesh.";
            }
            logfile << std::endl;
#endif
        }
        
        __ADD_CLONE_CODE_FOR_ABSTRACT_CLASS__(CLASS)
        
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
        
    protected :
        
        const AdaptivitySettings settings;
        const Real theta  = 10;
        const Real theta2 = 100;
        const Real intersection_theta2 = 10000000000.;
        
        //        const Real theta  =  10.;
        
        mutable Int block_count = 0;
        mutable Int max_level_reached = 0;
        mutable Int bottom_count = 0;
        mutable Int inadmissable_count = 0;
        mutable Int primitive_count = 0;
        mutable Real total_sum = static_cast<Real>(0);

        mutable std::ofstream logfile;
        mutable std::ofstream simplex_file;
        mutable std::ofstream center_file;
        mutable std::ofstream emb_simplex_file;
        mutable std::ofstream emb_center_file;
        
    public:
        

        virtual void LoadS( const Int i ) override
        {
            S_ID = i;
            const Real * const X  = S_near + NEAR_DIMS * S_ID;
    
            S.RequireSimplex(S_serialized, S_ID);
            
            a = X[0];
        
#ifdef NearField_S_Copy
            copy_buffer( &X[1],            &x_buffer[0], COORD_DIMS );
            copy_buffer( &X[1+COORD_DIMS], &p[0]       , PROJ_DIM   );
#else
            x_buffer = &X[1];
            p        = &X[1+AMB_DIM];
#endif
        }
        
        virtual void LoadT( const Int j ) override
        {
            T_ID = j;
            const Real * const Y  = T_near + NEAR_DIMT * T_ID;
    
            T.RequireSimplex(T_serialized, T_ID);
            
            b = Y[0];
        
#ifdef NearField_T_Copy
            copy_buffer( &Y[1],            &y_buffer[0], COORD_DIMT );
            copy_buffer( &Y[1+COORD_DIMT], &q[0]       , PROJ_DIM   );
#else
            y_buffer = &Y[1];
            q        = &Y[1+COORD_DIMT];
#endif
        }

        virtual Real Energy() override
        {
            Real sum = static_cast<Real>(0);
            
            ++primitive_count;
            
            bool from_above = true;
            bool shall_continue = true;
        
            S.ToChild(0);
            T.ToChild(0);
        
            while( shall_continue )
            {
                if( from_above )
                {
                    if( S.Level() >= settings.max_level )
                    {
                        // If at lowest level and inadmissable then we just compute the energy and move up.
                        max_level_reached = settings.max_level;
                        block_count++;
                        bottom_count++;
                        sum += energy();
                        S.ToParent();
                        T.ToParent();
                        from_above = false;
                    }
                    else
                    {
//                         If not at lowest level, then we have to check for admissability.
                        auto & P = S.SimplexPrototype();
                        auto & Q = T.SimplexPrototype();
                        
                        ADMISSABILITY_CONDITION;
                        
                        if( admissable )
                        {
                            // We compute energy, go to parent, and prepare the next child of the parent.
                            max_level_reached = std::max( max_level_reached, S.Level() );
                            block_count++;
                            sum += energy();
                            S.ToParent();
                            T.ToParent();
                            from_above = false;
                        }
                        else
                        {
                            // If inadmissable, we go a level deeper.
                            S.ToChild(0);
                            T.ToChild(0);
                            from_above = true;
                        }
                    }
                }
                else
                {
                    // If we come from below, we have to find the next pair of simplices to visit.
                    
                    Int S_k = S.FormerChildID();
                    Int T_k = T.FormerChildID();
                    
//                    print("Coming from "+ToString(S_k)+"-th child of S and "+ToString(T_k)+"-th child of T.");
                    
                    if( T_k < T.ChildCount()-1 )
                    {
                        S.ToChild(S_k);
                        T.ToChild(T_k+1);
                        from_above = true;
                    }
                    else
                    {
                        if( S_k < S.ChildCount()-1 )
                        {
                            S.ToChild(S_k+1);
                            T.ToChild(0);
                            from_above = true;
                        }
                        else
                        {
                            // No further unvisited children. Either move up or break.
                            if( S.Level() == 0 )
                            {
//                                inadmissable_count++;
                                shall_continue = false;
                                from_above = true;
                            }
                            else
                            {
                                S.ToParent();
                                T.ToParent();
                                from_above = false;
                            }
                        }
                    }
                    
                } // if( from_above )
                
            } // while( shall_continue )
            
            total_sum += sum;
            
            return sum;
        }
    
        virtual Real DEnergy() override
        {
            Real sum = static_cast<Real>(0);
            
            ++primitive_count;
            
            bool from_above = true;
            bool shall_continue = true;
            
            S.ToChild(0);
            T.ToChild(0);
        
            while( shall_continue )
            {
                if( from_above )
                {
                    if( S.Level() >= settings.max_level )
                    {
                        // If at lowest level and inadmissable then we just compute the energy and move up.
                        max_level_reached = settings.max_level;
                        block_count++;
                        bottom_count++;
                        sum += denergy();
                        S.ToParent();
                        T.ToParent();
                        from_above = false;
                    }
                    else
                    {
                        // If not at lowest level, then we have to check for admissability.
                        auto & P = S.SimplexPrototype();
                        auto & Q = T.SimplexPrototype();
                        
                        ADMISSABILITY_CONDITION
                        
                        if( admissable )
                        {
                            // We compute energy, go to parent, and prepare the next child of the parent.
                            max_level_reached = std::max( max_level_reached, S.Level() );
                            block_count++;
                            sum += denergy();
                            S.ToParent();
                            T.ToParent();
                            from_above = false;
                        }
                        else
                        {
                            // If inadmissable, we go a level deeper.
   
                            S.ToChild(0);
                            T.ToChild(0);
                            from_above = true;
                        }
                    }
                }
                else
                {
                    // If we come from below, we have to find the next pair of simplices to visit.
                    
                    Int S_k = S.FormerChildID();
                    Int T_k = T.FormerChildID();
                    
//                    print("Coming from "+ToString(S_k)+"-th child of S and "+ToString(T_k)+"-th child of T.");
                    
                    if( T_k < T.ChildCount()-1 )
                    {
                        S.ToChild(S_k);
                        T.ToChild(T_k+1);
                        from_above = true;
                    }
                    else
                    {
                        if( S_k < S.ChildCount()-1 )
                        {
                            S.ToChild(S_k+1);
                            T.ToChild(0);
                            from_above = true;
                        }
                        else
                        {
                            // No further unvisited children. Either move up or break.
                            if( S.Level() == 0 )
                            {
                                shall_continue = false;
                                from_above = true;
                            }
                            else
                            {
                                S.ToParent();
                                T.ToParent();
                                from_above = false;
                            }
                        }
                    }
                    
                } // if( from_above )
                
            } // while( shall_continue )

            total_sum += sum;
            
            return sum;
        }
        
        virtual Real  energy() const = 0;
        
        virtual Real denergy() const = 0;
        
        virtual bool IsRepulsive() const = 0;
        
    public:
        
        virtual std::string Stats() const override
        {
            return ClassName();
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM1)+","+ToString(DOM_DIM2)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
  
    };
    
} // namespace Repulsor

#undef CLASS
#undef BASE

#undef ADMISSABILITY_CONDITION
