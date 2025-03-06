#pragma once

namespace Repulsor
{
    template<
        int S_DOM_DIM_, int T_DOM_DIM_,
        typename ClusterTree_T_,
        bool symmetricQ_,
        bool energy_flag_, bool diff_flag_, bool metric_flag_, bool density_flag_
    >
    class FMM_Kernel_NF : public FMM_Kernel<ClusterTree_T_,symmetricQ_,energy_flag_,diff_flag_,metric_flag_, density_flag_>
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
        
        static constexpr Int S_DOM_DIM = S_DOM_DIM_;
        static constexpr Int T_DOM_DIM = T_DOM_DIM_;
        
        static constexpr Int S_COORD_DIM = (S_DOM_DIM+1)*AMB_DIM;
        static constexpr Int T_COORD_DIM = (T_DOM_DIM+1)*AMB_DIM;
        static constexpr Int S_DATA_DIM  = 1 + S_COORD_DIM + PROJ_DIM;
        static constexpr Int T_DATA_DIM  = 1 + T_COORD_DIM + PROJ_DIM;
        
    protected:
        
        mptr<Real> metric_data = nullptr;
        
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
        
        cptr<Real> S_data    = nullptr;
        mptr<Real> S_D_data  = nullptr;
        mptr<Real> S_diag    = nullptr;
        
        cptr<Real> T_data    = nullptr;
        mptr<Real> T_D_data  = nullptr;
        mptr<Real> T_diag    = nullptr;
        
        using Base_T::thread;
    
        static const constexpr Real S_scale = Inv<Real>(S_DOM_DIM+1);
        static const constexpr Real T_scale = Inv<Real>(T_DOM_DIM+1);
        
    public:
        
        FMM_Kernel_NF() = default;
        
        // To be used for configuration of kernel.
        explicit FMM_Kernel_NF( mref<Configurator_T> conf, const Int thread_ )
        :   Base_T      ( conf, thread_                                       )
        ,   S_data      ( GetS().PrimitiveNearFieldData().data()              )
        ,   S_D_data    ( GetS().ThreadPrimitiveDNearFieldData().data(thread) )
        ,   S_diag      ( GetS().NF_Accumulator().data(               thread) )
        ,   T_data      ( GetT().PrimitiveNearFieldData().data()              )
        ,   T_D_data    ( GetT().ThreadPrimitiveDNearFieldData().data(thread) )
        ,   T_diag      ( GetT().NF_Accumulator().data(               thread) )
        {
            TOOLS_DEBUG_PRINT(std::string( "Initializing "+ this->ClassName() + " from Configurator_T on thread " + ToString(thread)) );
            
            if( GetS().PrimitiveNearFieldData().Dim(1) != S_DATA_DIM )
            {
                eprint(ClassName()+" Constructor: GetS().PrimitiveNearFieldData().Dim(1) != S_DATA_DIM");
            }
            
            if constexpr ( diff_flag || density_flag )
            {
                if( GetS().ThreadPrimitiveDNearFieldData().Dim(2) != S_DATA_DIM )
                {
                    eprint(ClassName()+" Constructor: GetS().ThreadPrimitiveDNearFieldData().Dim(2) != S_DATA_DIM");
                }
            }
               
            if( GetT().PrimitiveNearFieldData().Dim(1) != T_DATA_DIM)
            {
               eprint(ClassName()+" Constructor: GetT().PrimitiveNearFieldData().Dim(1) != T_DATA_DIM ");
            }
            
            if constexpr ( diff_flag || density_flag )
            {
                if( GetT().ThreadPrimitiveDNearFieldData().Dim(2) != T_DATA_DIM )
                {
                    eprint(ClassName()+" Constructor: GetT().ThreadPrimitiveDNearFieldData().Dim(2) != T_DATA_DIM ");
                }
            }
        }
        
        FMM_Kernel_NF( mref<FMM_Kernel_NF> other, const Int thread_ )
        :   Base_T      ( other, thread_                                            )
        ,   metric_data ( other.OffDiag().data()                                    )
        ,   S_data      ( other.S_data                                              )
        ,   S_D_data    ( other.GetS().ThreadPrimitiveDNearFieldData().data(thread) )
        ,   S_diag      ( other.GetS().NF_Accumulator().data(               thread) )
        ,   T_data      ( other.T_data                                              )
        ,   T_D_data    ( other.GetT().ThreadPrimitiveDNearFieldData().data(thread) )
        ,   T_diag      ( other.GetT().NF_Accumulator().data(               thread) )
        {
            TOOLS_DEBUG_PRINT(std::string( "Initializing " + this->ClassName() + " from "+ClassName()+" on thread " + ToString(thread)) );
        }
        
        ~FMM_Kernel_NF() = default;

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
            TOOLS_MAKE_FP_FAST()
            
            cptr<Real> X = &S_data[S_DATA_DIM * i_global];
            
            a = X[0];
            
#ifdef NearField_S_Copy
            P.Read( &X[1+S_COORD_DIM] );
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
            
            if constexpr ( diff_flag || density_flag )
            {
                DX.SetZero();
            }
        }
        
        TOOLS_FORCE_INLINE void loadT( const Int j_global )
        {
            TOOLS_MAKE_FP_FAST()
            
            cptr<Real> Y = &T_data[T_DATA_DIM * j_global];
            
            b = Y[0];
            
#ifdef NearField_T_Copy
            Q.Read( &Y[1+T_COORD_DIM] );
#else
            Q = &Y[1+T_COORD_DIM];
#endif
            for( Int k = 0; k < AMB_DIM; ++k )
            {
                y[k] = T_scale * Y[1 + AMB_DIM * 0 + k];
                
                for( Int kk = 1; kk < T_DOM_DIM+1; ++kk )
                {
                    y[k] += T_scale * Y[1 + AMB_DIM * kk + k];
                }
            }
            
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
            return this->metric_values.NF;
        }
        
        cref<Values_T> OffDiag() const
        {
            return this->metric_values.NF;
        }
        
        mref<Values_T> Diag()
        {
            return this->metric_values.NF_diag;
        }
        
        cref<Values_T> Diag() const
        {
            return this->metric_values.NF_diag;
        }

        
        void Reduce( cref<FMM_Kernel_NF> ker )
        {
            (void)ker;
        }
        
        std::string ClassName() const
        {
            return std::string("FMM_Kernel_NF<")
            + ToString(S_DOM_DIM) + ","
            + ToString(T_DOM_DIM) + ","
//            + GetS().ClassName() + ","
            + + "...,"
            + ToString(energy_flag) + ","
            + ToString(diff_flag) + ","
            + ToString(metric_flag) + ","
            + ToString(density_flag) +
            ">";
        }

    };

} // namespace Repulsor
