#pragma once

#define CLASS SquaredDistanceMatrix__NFK_Adaptive
#define BASE  Metric__NFK_Adaptive<DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Int,SReal>

namespace Repulsor
{
    template<int DOM_DIM1, int DOM_DIM2, int AMB_DIM, typename Real, typename Int, typename SReal>
    class CLASS : public BASE
    {
    protected:
        
        using BASE::S_near;
        using BASE::S_D_near;
        using BASE::S_serialized;
        
        using BASE::T_near;
        using BASE::T_D_near;
        using BASE::T_serialized;
        
        using BASE::a;
        using BASE::x_buffer;
        using BASE::p;
        
        using BASE::b;
        using BASE::y_buffer;
        using BASE::q;

        using BASE::settings;
        using BASE::S;
        using BASE::T;
        
        using BASE::S_ID;
        using BASE::T_ID;
        
        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
        using BASE::lambda;
        using BASE::mu;

    public:
        
        using BASE::NEAR_DIMS;
        using BASE::NEAR_DIMT;
        using BASE::COORD_DIMS;
        using BASE::COORD_DIMT;
        using BASE::PROJ_DIM;
        
    public:
        
        CLASS() : BASE() {}

        // Copy constructor
        CLASS( const CLASS & other )
        :   BASE( other )
        ,   values( other.values )
        {}

        // Move constructor
        CLASS( CLASS && other ) noexcept 
        :   BASE( other )
        ,   values( other.values )
        {
            other.values = nullptr;
        }
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)
        
    protected:
                        
        Real * restrict values = nullptr;
        
        Real value = static_cast<Real>(0);
        
    public:
        
        virtual void AllocateValueBuffers(
            std::map<KernelType, Tensor1<Real,Int>> & near_values, const Int nnz ) override
        {
            near_values[KernelType::SquaredDistance] = Tensor1<Real,Int> (nnz);
            
            values = near_values[KernelType::SquaredDistance].data();
        }
        
        virtual void LoadValueBuffers(
            std::map<KernelType, Tensor1<Real,Int>> & near_values ) override
        {
            values = near_values[KernelType::SquaredDistance].data();
        }
        
        virtual void LoadS( const Int i ) override
        {
            S_ID = i;
            const Real * const X  = &S_near[NEAR_DIMS * S_ID];
    
            S.RequireSimplex(S_serialized, S_ID);
            
            a = X[0];
        
#ifdef NearField_S_Copy
            copy_buffer( &X[1], &x_buffer[0], COORD_DIMS );
#else
            x_buffer = &X[1];
#endif
        }
        
        virtual void LoadT( const Int j ) override
        {
            T_ID = j;
            const Real * const Y  = T_near + NEAR_DIMT * T_ID;
    
            T.RequireSimplex(T_serialized, T_ID);
            
            b = Y[0];
        
#ifdef NearField_T_Copy
            copy_buffer( &Y[1], &y_buffer[0], COORD_DIMT );
#else
            y_buffer = &Y[1];
#endif
        }
        
        
        virtual void StartRow() override {}
        
        virtual void StartEntry() override
        {
            value = static_cast<Real>(0);
        }
        
        virtual void FinishEntry( const Int pos ) override
        {
            values[pos] = value;
        }

        
        virtual void FinishRow( const Int diag_pos ) override {}

        virtual void metric() override
        {
            Real x [AMB_DIM] = {};
            Real y [AMB_DIM] = {};
            Real v [AMB_DIM] = {};
            
            Real r2 = static_cast<Real>(0);

            for( Int i = 0; i < AMB_DIM; ++i )
            {
                x[i] = lambda[0] * x_buffer[AMB_DIM*0+i];
                y[i] = mu    [0] * y_buffer[AMB_DIM*0+i];

                for( Int ii = 1; ii < DOM_DIM1+1; ++ii )
                {
                    x[i] += lambda[ii] * x_buffer[AMB_DIM*ii+i];
                }
                
                for( Int ii = 1; ii < DOM_DIM2+1; ++ii )
                {
                    y[i] += mu    [ii] * y_buffer[AMB_DIM*ii+i];
                }
            }
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                v[i] = y[i] - x[i];
                r2  += v[i] * v[i];
            }

            const Real w = S.Weight() * T.Weight();
            
            value += w * r2;
        }
        
    public:
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM1)+","+ToString(DOM_DIM2)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+">";
        }
  
    };
    
} // namespace Repulsor

#undef BASE
#undef CLASS

