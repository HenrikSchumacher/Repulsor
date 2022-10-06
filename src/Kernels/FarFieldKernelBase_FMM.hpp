#pragma once

#define CLASS FarFieldKernelBase_FMM

// This sets up the basic I/O routines for all FMM-like far field kernels, both for energies and metrics.
namespace Repulsor
{
    template<int AMB_DIM, int DEGREE, typename Real, typename Int>
    class alignas( OBJECT_ALIGNMENT ) CLASS
    {
        ASSERT_FLOAT(Real);
        ASSERT_INT  (Int );
        
    public:
        
        static constexpr Int PROJ_DIM = (AMB_DIM*(AMB_DIM+1))/2;
        static constexpr Int FAR_DIM  = 1 + AMB_DIM + PROJ_DIM;
        
        CLASS()
        {
            Init();
        }
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   S_far       ( other.S_far   )
        ,   S_D_far     ( other.S_D_far )
        ,   T_far       ( other.T_far   )
        ,   T_D_far     ( other.T_D_far )
        {
            copy_buffer( &other.tri_i[0],    &tri_i[0],    PROJ_DIM        );
            copy_buffer( &other.tri_j[0],    &tri_j[0],    PROJ_DIM        );
            copy_buffer( &other.lin_k[0][0], &lin_k[0][0], AMB_DIM*AMB_DIM );
        }

        virtual ~CLASS() = default;
        
        __ADD_CLONE_CODE_FOR_BASE_CLASS__(CLASS)
        
    public:
        
        void Init()
        {
            Int k = 0;
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                lin_k[i][i] = k;
                tri_i[k] = i;
                tri_j[k] = i;
                ++k;
                for( Int j = i+1; j < AMB_DIM; ++j )
                {
                    tri_i[k] = i;
                    tri_j[k] = j;
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
        
    protected :
        
        const Real * restrict S_far   = nullptr;
              Real * restrict S_D_far = nullptr;
        
        const Real * restrict T_far   = nullptr;
              Real * restrict T_D_far = nullptr;
        
#ifdef FarField_FMM_S_Copy
        mutable Real a = 0.;
        mutable Real x [AMB_DIM]  = {};
        mutable Real p [PROJ_DIM] = {};
#else
        mutable Real a = 0.;
        mutable Real const * restrict x = nullptr;
        mutable Real const * restrict p = nullptr;
#endif
        
#ifdef FarField_FMM_T_Copy
        mutable Real b = 0.;
        mutable Real y [AMB_DIM]  = {};
        mutable Real q [PROJ_DIM] = {};
#else
        mutable Real b = 0.;
        mutable Real const * restrict y = nullptr;
        mutable Real const * restrict q = nullptr;
#endif
        
        mutable Real DX [1 + AMB_DIM + PROJ_DIM] = {};
        mutable Real DY [1 + AMB_DIM + PROJ_DIM] = {};

        Int tri_i [PROJ_DIM]         = {};
        Int tri_j [PROJ_DIM]         = {};
        Int lin_k [AMB_DIM][AMB_DIM] = {};
        
        mutable Int S_ID = -1;
        mutable Int T_ID = -1;
        
    public :
        
        bool LoadFarField( const Tensor2<Real,Int> & S_far_, const Tensor2<Real,Int> & T_far_ )
        {
            bool check = true;
            if( (S_far_.Dimension(1) != FAR_DIM) || (T_far_.Dimension(1) != FAR_DIM) )
            {
                eprint(ClassName()+"::LoadFarField : (S_far_.Dimension(1) != FarDim()) || (T_far_.Dimension(1) != FarDim()).");
                check = false;
            }
            else
            {
                S_far = S_far_.data();
                T_far = T_far_.data();
            }
            return check;
        }

        
        bool LoadDFarField( Tensor2<Real,Int> & S_D_far_, Tensor2<Real,Int> & T_D_far_ )
        {
            bool check = true;
            if( (S_D_far_.Dimension(1) != FAR_DIM) || (T_D_far_.Dimension(1) != FAR_DIM) )
            {
                eprint(ClassName()+"::LoadFarField : (S_D_far_.Dimension(1) != FarDim()) || (T_D_far_.Dimension(1) != FarDim()).");
                check = false;
            }
            else
            {
                S_D_far = S_D_far_.data();
                T_D_far = T_D_far_.data();
            }
            return check;
        }
        
        constexpr Int AmbDim() const
        {
            return AMB_DIM;
        }
        
        constexpr Int CoordDim() const
        {
            return AMB_DIM;
        }
        
        constexpr Int ProjectorDim() const
        {
            return PROJ_DIM;
        }
        
        constexpr Int FarDim() const
        {
            return FAR_DIM;
        }
        
        virtual void LoadS( const Int i )
        {
            S_ID = i;
            const Real * const restrict X = &S_far[FAR_DIM * S_ID];

            a = X[0];
#ifdef FarField_FMM_S_Copy
            copy_buffer( &X[1],         &x[0], AMB_DIM  );
            copy_buffer( &X[1+AMB_DIM], &p[0], PROJ_DIM );
#else
            x = &X[1];
            p = &X[1 + AMB_DIM];
#endif
        }
        
        
        virtual void CleanseDBufferS()
        {
            zerofy_buffer(&DX[0],FAR_DIM);
        }
        
        virtual void WriteDBufferS() const
        {
            Real * restrict const to = &S_D_far[FAR_DIM * S_ID];
            
            for( Int k = 0; k < FAR_DIM; ++k )
            {
                to[k] += DX[k];
            }
        }
        
        virtual void LoadT( const Int j )
        {
            T_ID = j;
            const Real * const restrict Y = &T_far[FAR_DIM * T_ID];
            
            b = Y[0];
#ifdef FarField_FMM_T_Copy
            copy_buffer( &Y[1],         &y[0], AMB_DIM  );
            copy_buffer( &Y[1+AMB_DIM], &q[0], PROJ_DIM );
#else
            y = &Y[1];
            q = &Y[1 + AMB_DIM];
#endif
        }

        virtual void CleanseDBufferT()
        {
            zerofy_buffer( &DY[0], FAR_DIM );
        }
        
        virtual void WriteDBufferT() const
        {
            Real * restrict const to = &T_D_far[FAR_DIM * T_ID];
            
            for( Int k = 0; k < FAR_DIM; ++k )
            {
                to[k] += DY[k];
            }
        }
        
    public:
        
        virtual std::string Stats() const
        {
            return this->ClassName();
        }
        
        virtual std::string ClassName() const
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
  
    };
    
} // namespace Repulsor

#undef CLASS

