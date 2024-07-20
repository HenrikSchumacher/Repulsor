#pragma once


namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename ClusterTree_T_,
        bool symmetricQ_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_, bool density_flag_
    >
    class FMM_Kernel_VF : public FMM_Kernel<
            ClusterTree_T_,symmetricQ_,
            energy_flag_,diff_flag_,metric_flag_,density_flag_
    >
    {
    private:
        
        using Base_T = FMM_Kernel<
            ClusterTree_T_,symmetricQ_,
            energy_flag_,diff_flag_,metric_flag_,density_flag_
        >;

    public:
        
        using ClusterTree_T      = ClusterTree_T_;
        using Real               = typename Base_T::Real;
        using SReal              = typename Base_T::SReal;
        using ExtReal            = typename Base_T::ExtReal;
        using Int                = typename Base_T::Int;
        using LInt               = typename Base_T::LInt;
        
        using Configurator_T     = typename Base_T::Configurator_T;
        using ValueContainer_T   = typename Base_T::ValueContainer_T;
        using Values_T           = typename ValueContainer_T::Values_T;
        
        
        using Base_T::AMB_DIM;
        using Base_T::PROJ_DIM;
        using Base_T::symmetry_factor;
        using Base_T::symmetricQ;
        using Base_T::energy_flag;
        using Base_T::diff_flag;
        using Base_T::metric_flag;
        using Base_T::density_flag;
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
        
//        using GJK_T    = GJK_Algorithm<AMB_DIM,GJK_Real,Int>;
        using GJK_T    = GJK<AMB_DIM,GJK_Real,Int>;
        
    protected:
        
        Real * restrict const metric_data = nullptr;
        
        mutable Real sum = Scalar::Zero<Real>;
        
        mutable Real a   = Scalar::Zero<Real>;
        mutable Real b   = Scalar::Zero<Real>;
        
        mutable Tiny::Vector<AMB_DIM,Real,Int> x;
#ifdef NearField_S_Copy
        mutable Tiny::Vector<PROJ_DIM,Real,Int> P;
#else
        mutable Real const * restrict P = nullptr;
#endif
        
        mutable Tiny::Vector<AMB_DIM,Real,Int> y;
#ifdef NearField_T_Copy
        mutable Tiny::Vector<PROJ_DIM,Real,Int> Q;
#else
        mutable Real const * restrict Q = nullptr;
#endif
        
        mutable Tiny::Vector<S_DATA_DIM,Real,Int> DX;
        mutable Tiny::Vector<T_DATA_DIM,Real,Int> DY;
        
        cptr< Real> S_data    = nullptr;
        mptr< Real> S_D_data  = nullptr;
        mptr< Real> S_diag    = nullptr;
        cptr<SReal> S_ser     = nullptr;
        
        cptr< Real> T_data    = nullptr;
        mptr< Real> T_D_data  = nullptr;
        mptr< Real> T_diag    = nullptr;
        cptr<SReal> T_ser     = nullptr;

#ifdef NearField_S_Copy
        mutable Tiny::Vector<S_DATA_DIM,Real,Int> x_buffer;
#else
        mutable Real const * restrict x_buffer  = nullptr;
#endif
        
#ifdef NearField_T_Copy
        mutable Tiny::Vector<T_DATA_DIM,Real,Int> y_buffer;
#else
        mutable Real const * restrict y_buffer  = nullptr;
#endif
        using Base_T::thread;
        
        mutable S_Tree_T S_Tree;
        mutable T_Tree_T T_Tree;
        
        cptr<SReal> lambda = S_Tree.Center();
        cptr<SReal> mu     = T_Tree.Center();
        
        mutable GJK_T gjk;
        
        const Real theta  = 10;
        const Real theta2 = 100;
        const typename S_Tree_T::Child_T max_refinement = 20;
        
        mutable Int block_count = 0;
        mutable typename S_Tree_T::Child_T max_level_reached = 0;
        mutable Int bottom_count = 0;
        mutable Int inadmissable_count = 0;
        mutable Int primitive_count = 0;
        mutable Real total_sum = Scalar::Zero<Real>;
        
        mutable Size_T evaluations = 0;
        
//#ifdef REPULSOR__PRINT_REPORTS_FOR_ADAPTIVE_KERNELS
//        mutable std::ofstream logfile;
//        mutable std::ofstream simplex_file;
//        mutable std::ofstream center_file;
//        mutable std::ofstream emb_simplex_file;
//        mutable std::ofstream emb_center_file;
//#endif
        
    public:
        
//        FMM_Kernel_VF() = default;
        FMM_Kernel_VF() = delete;
        
        FMM_Kernel_VF( mref<Configurator_T> conf, const Int thread_, const Real theta_, const Int max_refinement_  )
        :   Base_T         ( conf, thread_                                       )
        ,   S_data         ( GetS().PrimitiveNearFieldData().data()              )
        ,   S_D_data       ( GetS().ThreadPrimitiveDNearFieldData().data(thread) )
        ,   S_diag         ( GetS().VF_Accumulator().data(               thread) )
        ,   S_ser          ( GetS().PrimitiveSerialized().data()                 )
        ,   T_data         ( GetT().PrimitiveNearFieldData().data()              )
        ,   T_D_data       ( GetT().ThreadPrimitiveDNearFieldData().data(thread) )
        ,   T_diag         ( GetT().VF_Accumulator().data(               thread) )
        ,   T_ser          ( GetT().PrimitiveSerialized().data()                 )
        ,   theta          ( theta_                                              )
        ,   theta2         ( theta_ * theta_                                     )
        ,   max_refinement ( Min(
                                  int_cast<typename S_Tree_T::Child_T>(max_refinement_),
                                  Min(S_Tree_T::MaxLevel(),S_Tree_T::MaxLevel())
                              )                                                                )
        {
            debug_print(std::string( "Initializing " + this->ClassName() + " from Configurator_T on thread " + ToString(thread)) );
            
            if( GetS().PrimitiveSerialized().Dimension(1) != S_Tree.SimplexPrototype().Size() )
            {
                eprint(ClassName()+" Constructor: GetS().PrimitiveSerialized().Dimension(1) != S_Tree.SimplexPrototype().Size()");
            }
            
            if( GetT().PrimitiveSerialized().Dimension(1) != T_Tree.SimplexPrototype().Size() )
            {
                eprint(ClassName()+" Constructor: GetT().PrimitiveSerialized().Dimension(1) != T_Tree.SimplexPrototype().Size()");
            }
            
            if( max_refinement_ > max_refinement )
            {
                wprint(ClassName()+" Constructor: Reduced max_refinement from "+ToString(max_refinement_)+" to "+ToString(max_refinement)+" to prevent integer overflow.");
            }
        }
        
        FMM_Kernel_VF( mref<FMM_Kernel_VF> other, const Int thread_ )
        :   Base_T      ( other, thread_                                            )
        ,   metric_data ( other.OffDiag().data()                                    )
        ,   S_data      ( other.S_data                                              )
        ,   S_D_data    ( other.GetS().ThreadPrimitiveDNearFieldData().data(thread) )
        ,   S_diag      ( other.GetS().VF_Accumulator().data(               thread) )
        ,   S_ser       ( other.S_ser                                               )
        ,   T_data      ( other.T_data                                              )
        ,   T_D_data    ( other.GetT().ThreadPrimitiveDNearFieldData().data(thread) )
        ,   T_diag      ( other.GetT().VF_Accumulator().data(               thread) )
        ,   T_ser       ( other.T_ser                                               )
        ,   theta       ( other.theta                                               )
        ,   theta2      ( other.theta2                                              )
        ,   max_refinement      ( other.max_refinement                              )
        ,   max_level_reached   ( other.max_level_reached                           )
        ,   evaluations         ( other.evaluations                                 )
        {
            debug_print(std::string( "Initializing "+ClassName()+" from "+ClassName()+" on thread " + ToString(thread)) );
        }
        
        ~FMM_Kernel_VF()
        {
//#ifdef REPULSOR__PRINT_REPORTS_FOR_ADAPTIVE_KERNELS
//            logfile
//            << "\n"
//            << "Report for class                    = " << ClassName() << "\n"
//            << "Thread ID                           = " << thread << "\n"
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
            
            if constexpr (diff_flag || density_flag)
            {
                prefetch_buffer<T_DATA_DIM,1,0>( &T_D_data[T_DATA_DIM * j] );
            }
        }
        
    protected:

        force_inline void loadS( const Int i_global )
        {
            cptr<Real> X = &S_data[S_DATA_DIM * i_global];
    
            S_Tree.RequireSimplex(S_ser, i_global);
            
            a = X[0];
        
#ifdef NearField_S_Copy
            x_buffer.Read( &X[1] );
            P.Read( &X[1+S_COORD_DIM] );
#else
            x_buffer = &X[1];
            P        = &X[1+T_COORD_DIM];
#endif
            if constexpr ( diff_flag || density_flag )
            {
                DX.SetZero();
            }
        }
        
        force_inline void loadT( const Int j_global )
        {
            cptr<Real> Y = &T_data[T_DATA_DIM * j_global];
    
            T_Tree.RequireSimplex(T_ser, j_global);
            
            b = Y[0];
        
#ifdef NearField_T_Copy
            y_buffer.Read( &Y[1] );
            Q.Read( &Y[1+T_COORD_DIM] );
#else
            y_buffer = &Y[1];
            Q        = &Y[1+T_COORD_DIM];
#endif
            if constexpr ( diff_flag || density_flag )
            {
                DY.SetZero();
            }
        }
        
        
        force_inline void writeS( const Int i_global )
        {
            if constexpr ( diff_flag || density_flag )
            {
                combine_buffers<Scalar::Flag::Generic,Scalar::Flag::Plus,S_DATA_DIM>(
                    symmetry_factor, DX.data(), Scalar::One<Real>, &S_D_data[S_DATA_DIM * i_global]
                );
            }
        }
        
        force_inline void writeT( const Int j_global )
        {
            if constexpr ( diff_flag || density_flag )
            {
                combine_buffers<Scalar::Flag::Generic,Scalar::Flag::Plus,T_DATA_DIM>(
                    symmetry_factor, DY.data(), Scalar::One<Real>, &T_D_data[T_DATA_DIM * j_global]
                );
            }
        }
        
    public:
        
        mref<Values_T> OffDiag()
        {
            return this->metric_values.VF;
        }
        
        cref<Values_T> OffDiag() const
        {
            return this->metric_values.VF;
        }
        
        mref<Values_T> Diag()
        {
            return this->metric_values.VF_diag;
        }
        
        cref<Values_T> Diag() const
        {
            return this->metric_values.VF_diag;
        }
        
//        void CreateLogFile() const
//        {
//            DUMP(ClassName());
//            logprint("CreateLogFile");
//
//#ifdef REPULSOR__PRINT_REPORTS_FOR_ADAPTIVE_KERNELS
//
//            logprint("Really creating log file");
//            std:: string s = "./Repulsor__"+ClassName()+"_Report_"+ToString(thread)+".txt";
//
//            DUMP(s);
//
//            logfile.open(s, std::ios_base::app);
//
//            logfile << "Log file for " << ClassName() << std::endl;
//
//            s = "./Repulsor__"+ClassName()+"_Simplices_"+ToString(thread)+".txt";
//            simplex_file.open(s, std::ios_base::app);
//
//            s = "./Repulsor__"+ClassName()+"_Centers_"+ToString(thread)+".txt";
//            center_file.open(s, std::ios_base::app);
//
//            s = "./Repulsor__"+ClassName()+"_EmbSimplices_"+ToString(thread)+".txt";
//            emb_simplex_file.open(s, std::ios_base::app);
//
//            s = "./Repulsor__"+ClassName()+"_EmbCenters_"+ToString(thread)+".txt";
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
        
        Size_T EvaluationCount() const
        {
            return evaluations;
        }
        
        void Reduce( cref<FMM_Kernel_VF> ker )
        {
            max_level_reached = Max( max_level_reached, ker.MaxLevelReached() );
            
            evaluations += ker.EvaluationCount();
        }
        
        std::string ClassName() const
        {
            return std::string("FMM_Kernel_VF<")
            + ToString(S_DOM_DIM) + ","
            + ToString(T_DOM_DIM) + ","
//            + GetS().ClassName() + ","
            + "...,"
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(metric_flag) +  ","
            + ToString(density_flag) +
            ">";
        }

    };

} // namespace Repulsor

