#pragma once

#define CLASS FMM_Kernel

namespace Repulsor
{
    template<
        typename ClusterTree_T_,
        bool is_symmetric_,
        bool energy_flag_, bool diff_flag_, bool hess_flag_, bool metric_flag_
    >
    class alignas( OBJECT_ALIGNMENT ) CLASS
    {
    public:
        
        using ClusterTree_T = ClusterTree_T_;
        using Root_T        = CLASS;
        
        using Real    = typename ClusterTree_T::Real;
        using Int     = typename ClusterTree_T::Int;
        using SReal   = typename ClusterTree_T::SReal;
        using ExtReal = typename ClusterTree_T::ExtReal;
        
        using ValueContainer_T = Tensor2<Real,Int>;
        
        using Configurator_T = FMM_Configurator<ClusterTree_T>;
        
        static constexpr bool is_symmetric = is_symmetric_;
        
        static constexpr Int AMB_DIM   = ClusterTree_T::AMB_DIM;
        static constexpr Int PROJ_DIM  = (AMB_DIM*(AMB_DIM+1))/2;
        
        static constexpr bool energy_flag = energy_flag_;
        static constexpr bool diff_flag   = diff_flag_;
        static constexpr bool hess_flag   = hess_flag_;
        static constexpr bool metric_flag = metric_flag_;
        
    protected:
        
        static constexpr Real zero = static_cast<Real>(0);
        static constexpr Real one  = static_cast<Real>(1);
        static constexpr Real two  = static_cast<Real>(2);
        
        static constexpr Real symmetry_factor = one / (one + static_cast<Real>(is_symmetric) );
                                                        
        Int tri_i [PROJ_DIM] = {};
        Int tri_j [PROJ_DIM] = {};
        Int lin_k [AMB_DIM][AMB_DIM] = {};
        
        void Init()
        {
            Int k = 0;
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                lin_k[i][i] = k;
                tri_i[k]    = i;
                tri_j[k]    = i;
                ++k;
                for( Int j = i+1; j < AMB_DIM; ++j )
                {
                    tri_i[k]    = i;
                    tri_j[k]    = j;
                    lin_k[i][j] = lin_k[j][i] = k;
                    ++k;
                }
            }
        }
        
        void CheckInit()
        {
            print(ClassName()+"::CheckInit");
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                for( Int j = 0; j < AMB_DIM; ++j )
                {
                    print("{ "+ToString(i)+"," +ToString(j)+" } -> " + ToString(lin_k[i][j]) );
                }
            }

            for( Int k = 0; k < PROJ_DIM; ++k )
            {
                print( ToString(k) +" -> { "+ToString(tri_i[k])+"," +ToString(tri_j[k])+" }");
            }
        }
        
        
    public:

        CLASS() = delete;
        
        // To be used for configuration of kernel. CANNOT BE USED FOR COMPUTE MODE!
        CLASS( Configurator_T & conf )
        :   S             ( conf.GetS()                 )
        ,   T             ( conf.GetT()                 )
        ,   metric_values ( conf.MetricValues()         )   // In configure mode, kernels needs access to refs.
        ,   prec_values   ( conf.PreconditionerValues() )   // for allocation.
        {
            Init();
        }
        
        // Copy constructor. Must be used for compute mode!
        CLASS( const CLASS & other )
        :   S                  ( other.S                  )
        ,   T                  ( other.T                  )
        ,   metric_values      ( other.metric_values      )
        ,   metric_values_data ( other.metric_values_data ) // In compute the pointers are also needed!
        ,   prec_values        ( other.prec_values        )
        ,   prec_values_data   ( other.prec_values_data   ) // In compute the pointers are also needed!
        {
            Init();
        }
        
        virtual ~CLASS() = default;

    public:
        
        void Allocate( const Int nnz )
        {
            if constexpr ( metric_flag )
            {
                if(
                   metric_values.Dimensions(0) != nnz
                   ||
                   metric_values.Dimensions(1) != MetricNonzeroCount()
                   )
                {
                    metric_values = Tensor2<Real,Int>( nnz, MetricNonzeroCount() );
                }
                
                if(
                   prec_values.Dimensions(0) != nnz
                   ||
                   prec_values.Dimensions(1) != PreconditionerNonzeroCount()
                   )
                {
                    prec_values = Tensor2<Real,Int>( nnz, PreconditionerNonzeroCount() );
                }
            }
        }
        
        bool PointersInititializedQ() const
        {
            if( metric_values_data == nullptr )
            {
                eprint(ClassName()+"::PointersInititializedQ: Pointers for metric_values not initialized. Make sure to initialize the kernel via the copy constructor.");
                
                return false;
            }
            
            if( metric_values_data == nullptr )
            {
                eprint(ClassName()+"::PointersInititializedQ: Pointers for prec_values not initialized. Make sure to initialize the kernel via the copy constructor.");
                
                return false;
            }
            
            return true;
        }
        
        virtual Int MetricNonzeroCount() const = 0;

        virtual Int PreconditionerNonzeroCount() const = 0;
        
        virtual void LoadS( const Int i ) = 0;
        
        virtual void LoadT( const Int j ) = 0;

//        void PrefetchS( const Int i ) const
//        {}
        
        virtual void PrefetchT( const Int j ) const = 0;
        
        virtual Real compute( const Int block_ID ) = 0;
                                                       
        virtual Real Compute( const Int block_ID ) = 0;
        
        virtual void WriteS() = 0;
        
        virtual void WriteT() = 0;
        
    protected:
        
        const ClusterTree_T & S;
        const ClusterTree_T & T;
        
        ValueContainer_T & metric_values;
        ValueContainer_T & prec_values;
        
        Real * restrict const metric_values_data = nullptr;
        Real * restrict const   prec_values_data = nullptr;
        
        Int S_ID = -1;
        Int T_ID = -1;
        
    public:
        
        virtual std::string ClassName() const
        {
            return className();
        }
        
    private:
        
        std::string className() const
        {
            return TO_STD_STRING(CLASS)+"<"+S.ClassName()+">";
        }

    };

} // namespace Repulsor

#undef CLASS
