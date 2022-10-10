#pragma once

#define CLASS Metric__NFK_Adaptive
#define BASE  Metric__NFK<DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Int,SReal>


namespace Repulsor
{
    template<int DOM_DIM1, int DOM_DIM2, int AMB_DIM, typename Real, typename Int, typename SReal>
    class CLASS : public BASE
    {
    public:
        
        using BASE::NEAR_DIMS;
        using BASE::NEAR_DIMT;
        using BASE::COORD_DIMS;
        using BASE::COORD_DIMT;
        using BASE::PROJ_DIM;
        
    protected:
        
        using BASE::S;
        using BASE::S_serialized;
        using BASE::S_near;
//        using BASE::S_D_near;
        
        using BASE::T;
        using BASE::T_serialized;
        using BASE::T_near;
//        using BASE::T_D_near;
        
//        using BASE::DX;
//        using BASE::DY;
        
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
        
        using BASE::lambda;
        using BASE::mu;

        
    public:
        
        explicit CLASS( const AdaptivitySettings & settings_ = AdaptivitySettings() )
        : BASE()
        , settings(settings_)
        , theta2(settings_.theta * settings_.theta)
        , intersection_theta2( settings_.intersection_theta * settings_.intersection_theta)
        {}

        // Copy constructor
        CLASS( const CLASS & other ) : BASE(other)
        , settings      ( other.settings      )
        , theta2        ( other.settings.theta * other.settings.theta )
        , intersection_theta2( other.settings.intersection_theta * other.settings.intersection_theta)
        {}
        
        
        virtual ~CLASS() override = default;
        
//        {
//#ifdef REPULSOR__PRINT_REPORTS_FOR_ADAPTIVE_KERNELS
//            if( block_count > 0 )
//            {
//                const Int thread = omp_get_thread_num();
//
//                std::string filename = "./Repulsor__"+this->ClassName()+"_Report_"+ToString( thread )+".txt";
//                std::ofstream s (filename);
//                s
//                    << "Report for class                  = " << this->ClassName() << "\n"
//                    << "Thread ID                         = " << thread << "\n"
////                    << "Colliding  simplex pairs found    = " << collision_count << "\n"
////                    << "Neighboring simplex pairs found   = " << neighbor_count << "\n"
////                    << "Admissable simplex pairs found    = " << admissable_count << "\n"
//                    << "Total number of quadrature points = " << block_count << "\n"
//                    << "Subdivision level reached         = " << max_reached_level << "\n"
//                ;
//                if( bottom_count > 0 )
//                {
//                    s << "WARNING: Maximal subdivision level = "+ToString(settings.max_level)+" reached "+ToString(bottom_count)+" times. Expect non-sufficent repulsive behavior. Consider refining the mesh.";
//                }
//                s << std::endl;
//            }
//#endif
//        }
        
        __ADD_CLONE_CODE_FOR_ABSTRACT_CLASS__(CLASS)
        
    protected:
        
        const AdaptivitySettings settings;
        const Real theta2 = 100.;
        const Real intersection_theta2 = 10000000000.;
        
        //        const Real theta  =  10.;
        
        mutable Int block_count = 0;
        mutable Int max_reached_level = 0;
        mutable Int bottom_count = 0;
        mutable Int primitive_count = 0;
        mutable Int collision_count = 0;
        

    public:
        
        virtual void LoadS( const Int i ) override
        {
            S_ID = i;
            const Real * const X = &S_near[NEAR_DIMS * S_ID];
    
            S.RequireSimplex(S_serialized, S_ID);
            
            a = X[0];
        
#ifdef NearField_S_Copy
            copy_buffer( &X[1],            &x_buffer[0], COORD_DIMS );
            copy_buffer( &X[1+COORD_DIMS], &p[0]       , PROJ_DIM   );
#else
            x_buffer = &X[1];
            p        = &X[1+COORD_DIMS];
#endif

//            for( Int k = 0; k < DOM_DIM+1; ++k )
//            {
//                lambda[k] = static_cast<Real>(S.Center()[k]);
//            }
        }
        
        virtual void LoadT( const Int j ) override
        {
            T_ID = j;
            const Real * const Y = &T_near[NEAR_DIMT * T_ID];
    
            T.RequireSimplex(T_serialized, T_ID);
            
            b = Y[0];
        
#ifdef NearField_T_Copy
            copy_buffer( &Y[1],            &y_buffer[0], NEAR_DIMT );
            copy_buffer( &Y[1+COORD_DIMT], &q[0]       , PROJ_DIM  );
#else
            y_buffer = &Y[1];
            q        = &Y[1+COORD_DIMT];
#endif
            
//            for( Int k = 0; k < DOM_DIM+1; ++k )
//            {
//                mu[k] = static_cast<Real>(T.Center()[k]);
//            }
        }
        
        virtual void computeBlock() = 0;

        virtual void StartEntry() = 0;

        virtual void FinishEntry( const Int pos ) = 0;
        
        virtual void ComputeBlock( const Int pos ) override
        {
            StartEntry();
            
            ++primitive_count;
            
            bool from_above = true;
            bool shall_continue = true;
            
            // This assumes that the primitive pair has to be subdivided (otherwise the calling function would not call this kernel).
            
            S.ToChild(0);
            T.ToChild(0);
        
            while( shall_continue )
            {
                if( from_above )
                {
                    if( S.Level() >= settings.max_level )
                    {
                        // If at lowest level and inadmissable then we just compute the energy and move up.
                        max_reached_level = settings.max_level;
                        block_count++;
                        bottom_count++;
                        computeBlock();
                        S.ToParent();
                        T.ToParent();
                        from_above = false;
                    }
                    else
                    {
                        // If not at lowest level, then we have to check for admissability.
                        
                        const bool admissable = gjk.MultipoleAcceptanceCriterion(
                            S.SimplexPrototype(),
                            T.SimplexPrototype(),
                            theta2
                        );
                        
                        if( admissable )
                        {
                            // We compute, go to parent, and prepare the next child of the parent.
                            max_reached_level = std::max( max_reached_level, S.Level() );
                            block_count++;
                            computeBlock();
                            S.ToParent();
                            T.ToParent();
                            from_above = false;
                        }
                        else
                        {
                            // If inadmissabile, we go a level deeper.

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
                
            } // while( true )
            
            FinishEntry(pos);
            
        } // ComputeBlock
        
        
        virtual void TransposeBlock( const Int from, const Int to ) override = 0;
        
        virtual void LoadOutputBuffer( const Real * restrict const Y ) override = 0;
        
        virtual void LoadInputBuffer ( const Real * restrict const X ) override = 0;

        virtual void ClearVector() override = 0;
        
        virtual void ApplyBlock( const Real alpha, const Int pos, const Int j ) override = 0;
        
        virtual void WriteVector( const Int i ) const override = 0;
        
        virtual Int BlockSize() const override = 0;
        
        virtual Int ValueSize() const override = 0;
        
        
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
