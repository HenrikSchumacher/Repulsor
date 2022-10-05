#pragma once

#define CLASS SquaredDistanceMatrix__FFK_FMM
#define BASE  Metric__FFK_FMM<AMB_DIM,DEGREE,Real,Int>

namespace Repulsor
{
    template<int AMB_DIM, int DEGREE, typename Real, typename  Int>
    class CLASS : public BASE
    {
    protected:
        
        using BASE::S_far;
        using BASE::T_far;
        
        using BASE::a;
        using BASE::x;
        using BASE::p;
        
        using BASE::b;
        using BASE::y;
        using BASE::q;
        
        using BASE::S_ID;
        using BASE::T_ID;
        
        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
    public:

        using BASE::LoadS;
        using BASE::LoadT;
        using BASE::FarDim;
        using BASE::CoordDim;
        using BASE::ProjectorDim;
        
    public:
        
        CLASS() : BASE() {}

        // Copy constructor
        CLASS( const CLASS & other )
        :   BASE( other )
        ,   values ( other.values )
        {}

        // Move constructor
        CLASS( CLASS && other ) noexcept 
        :   BASE   ( other )
        ,   values ( other.values )
        {
            other.values = nullptr;
        }
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)
        
    protected:
        
        Real * restrict values = nullptr;
        
    public:

        virtual void AllocateValueBuffers(
            std::map<KernelType, Tensor1<Real,Int>> & far_values, const Int nnz ) override
        {
            far_values[KernelType::SquaredDistance] = Tensor1<Real,Int> (nnz);
            
            values = far_values[KernelType::SquaredDistance].data();
        }
        
        virtual void LoadValueBuffers(
            std::map<KernelType, Tensor1<Real,Int>> & far_values ) override
        {
            values = far_values[KernelType::SquaredDistance].data();
        }
        
        
        virtual void Metric( const Int pos ) override
        {
            Real v [AMB_DIM] = {};
            
            Real r2 = static_cast<Real>(0);
            
            for( Int i = 0; i < AMB_DIM; ++i )
            {
                v[i] = y[i] - x[i];
                r2  += v[i] * v[i];
            }
            
            values[pos] = r2;
        }
        
    public:
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(AMB_DIM)+","+ToString(DEGREE)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
  
    };
    
} // namespace Repulsor

#undef CLASS
#undef BASE
