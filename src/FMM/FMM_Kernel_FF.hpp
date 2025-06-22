#pragma once

namespace Repulsor
{
    template<
        typename ClusterTree_T_,
        bool symmetricQ_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_, bool density_flag_
    >
    class FMM_Kernel_FF : public FMM_Kernel<ClusterTree_T_,symmetricQ_,energy_flag_,diff_flag_,metric_flag_,density_flag_>
    {
    private:
        
        using Base_T = FMM_Kernel<ClusterTree_T_,symmetricQ_,energy_flag_,diff_flag_,metric_flag_,density_flag_>;
      
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
        
        static constexpr Int S_DATA_DIM  = 1 + AMB_DIM + PROJ_DIM;
        static constexpr Int T_DATA_DIM  = 1 + AMB_DIM + PROJ_DIM;
        
    protected:
        
        
        mptr<Real> metric_data = nullptr;
        
        mutable Real sum = Scalar::Zero<Real>;
        
        mutable Real a   = Scalar::Zero<Real>;
        mutable Real b   = Scalar::Zero<Real>;
        
#ifdef FarField_S_Copy
        mutable Tiny::Vector<AMB_DIM, Real,Int> x;
        mutable Tiny::Vector<PROJ_DIM,Real,Int> P;
#else
        mutable Real const * restrict x = nullptr;
        mutable Real const * restrict P = nullptr;
#endif
        
#ifdef FarField_T_Copy
        mutable Tiny::Vector<AMB_DIM, Real,Int> y;
        mutable Tiny::Vector<PROJ_DIM,Real,Int> Q;
#else
        mutable Real const * restrict y = nullptr;
        mutable Real const * restrict Q = nullptr;
#endif
        
        mutable Tiny::Vector<S_DATA_DIM,Real,Int> DX;
        mutable Tiny::Vector<T_DATA_DIM,Real,Int> DY;
        
        cptr<Real> S_data   = nullptr;
        mptr<Real> S_D_data = nullptr;
        mptr<Real> S_diag   = nullptr;

        cptr<Real> T_data   = nullptr;
        mptr<Real> T_D_data = nullptr;
        mptr<Real> T_diag   = nullptr;
        
        using Base_T::thread;
        
    public:
        
        explicit FMM_Kernel_FF( mref<Configurator_T> conf, Int thread_ )
        :   Base_T      ( conf, thread_                                    )
        ,   S_data      ( GetS().ClusterFarFieldData().data()              )
        ,   S_D_data    ( GetS().ThreadClusterDFarFieldData().data(thread) )
        ,   S_diag      ( GetS().FF_Accumulator().data(            thread) )
        ,   T_data      ( GetT().ClusterFarFieldData().data()              )
        ,   T_D_data    ( GetT().ThreadClusterDFarFieldData().data(thread) )
        ,   T_diag      ( GetT().FF_Accumulator().data(            thread) )
        {
            TOOLS_DEBUG_PRINT(std::string( "Initializing " + this->ClassName()+" from Configurator_T on thread " + ToString(thread)) );
            
            if( GetS().ClusterFarFieldData().Dim(1) != S_DATA_DIM )
            {
                eprint(ClassName()+" Constructor: GetS().ClusterFarFieldData().Dim(1) != S_DATA_DIM");
            }
            
            if constexpr ( diff_flag || density_flag )
            {
                if( GetS().ThreadClusterDFarFieldData().Dim(2) != S_DATA_DIM )
                {
                    eprint(ClassName()+" Constructor: GetS().ThreadClusterDFarFieldData().Dim(2) != S_DATA_DIM");
                }
            }
               
            if( GetT().ClusterFarFieldData().Dim(1) != T_DATA_DIM)
            {
               eprint(ClassName()+" Constructor: GetT().ClusterFarFieldData().Dim(1) != T_DATA_DIM ");
            }
            
            if constexpr ( diff_flag || density_flag )
            {
                if( GetT().ThreadClusterDFarFieldData().Dim(2) != T_DATA_DIM )
                {
                    eprint(ClassName()+" Constructor: GetT().ThreadClusterDFarFieldData().Dim(2) != T_DATA_DIM ");
                }
            }
        }
        
        FMM_Kernel_FF( mref<FMM_Kernel_FF> other, Int thread_ )
        :   Base_T      ( other, thread_                                         )
        ,   metric_data ( other.OffDiag().data()                                 )
        ,   S_data      ( other.S_data                                           )
        ,   S_D_data    ( other.GetS().ThreadClusterDFarFieldData().data(thread) )
        ,   S_diag      ( other.GetS().FF_Accumulator().data(            thread) )
        ,   T_data      ( other.T_data                                           )
        ,   T_D_data    ( other.GetT().ThreadClusterDFarFieldData().data(thread) )
        ,   T_diag      ( other.GetT().FF_Accumulator().data(            thread) )
        {
            TOOLS_DEBUG_PRINT(std::string( "Initializing "+ClassName()+" from "+ClassName()+" on thread " + ToString(thread)) );
        }
        
        // Default constructor
        FMM_Kernel_FF() = delete;
        // Destructor
        virtual ~FMM_Kernel_FF() override = default;
        // Copy constructor
        FMM_Kernel_FF( const FMM_Kernel_FF & other ) = default;
        // Copy assignment operator
        FMM_Kernel_FF & operator=( const FMM_Kernel_FF & other ) = default;
        // Move constructor
        FMM_Kernel_FF( FMM_Kernel_FF && other ) = default;
        // Move assignment operator
        FMM_Kernel_FF & operator=( FMM_Kernel_FF && other ) = default;
        
    public:
        
        TOOLS_FORCE_INLINE void Prefetch( const Int j_next ) const
        {
            prefetch_buffer<T_DATA_DIM,0,0>( &T_data[T_DATA_DIM * j_next] );
            
            if constexpr ( diff_flag || density_flag )
            {
                prefetch_buffer<T_DATA_DIM,1,0>( &T_D_data[T_DATA_DIM * j_next] );
            }
        }

    protected:
        
        TOOLS_FORCE_INLINE void loadS( const Int i_global )
        {
            cptr<Real> X = &S_data[S_DATA_DIM * i_global];

            a = X[0];
#ifdef FarField_S_Copy
            x.Read( &X[1] );
            P.Read( &X[1+AMB_DIM] );
#else
            x = &X[1];
            P = &X[1 + AMB_DIM];
#endif
            
            if constexpr ( diff_flag || density_flag )
            {
                DX.SetZero();
            }
        }
        
        TOOLS_FORCE_INLINE void loadT( const Int j_global )
        {
            cptr<Real> Y = &T_data[T_DATA_DIM * j_global];
            
            b = Y[0];
#ifdef FarField_T_Copy
            y.Read( &Y[1] );
            Q.Read( &Y[1+AMB_DIM] );
#else
            y = &Y[1];
            Q = &Y[1 + AMB_DIM];
#endif
            if constexpr ( diff_flag || density_flag )
            {
                DY.SetZero();
            }
        }
            
        TOOLS_FORCE_INLINE void writeS( const Int i_global )
        {
            if constexpr ( diff_flag || density_flag )
            {
                combine_buffers<Scalar::Flag::Generic,Scalar::Flag::Plus,S_DATA_DIM>(
                    symmetry_factor, DX.data(), Scalar::One<Real>, &S_D_data[S_DATA_DIM * i_global]
                );
            }
        }
        
        TOOLS_FORCE_INLINE void writeT( const Int j_global )
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
            return this->metric_values.FF;
        }
        
        cref<Values_T> OffDiag() const
        {
            return this->metric_values.FF;
        }
        
        mref<Values_T> Diag()
        {
            return this->metric_values.FF_diag;
        }
        
        cref<Values_T> Diag() const
        {
            return this->metric_values.FF_diag;
        }

        void Reduce( cref<FMM_Kernel_FF> ker )
        {
            (void)ker;
        }
        
        std::string ClassName() const
        {
            return std::string("FMM_Kernel_FF<")
//            + GetS().ClassName() + ","
            + "...,"
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(metric_flag) + ","
            + ToString(density_flag) +
            ">";
        }

    };

} // namespace Repulsor

