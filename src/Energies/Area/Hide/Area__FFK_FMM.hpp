#pragma once

#define CLASS Area__FFK_FMM
#define BASE  Energy__FFK_FMM<AMB_DIM,DEGREE,Real,Int>

namespace Repulsor
{
    template<int AMB_DIM, int DEGREE, typename Real, typename Int>
    class CLASS : public BASE
    {
    protected:
        
        using BASE::S_far;
        using BASE::S_D_far;
        
        using BASE::T_far;
        using BASE::T_D_far;
        
        using BASE::DX;
        using BASE::DY;
        
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
    
        CLASS() : BASE() {}
        
        // Copy constructor
        CLASS( const Area__FFK_FMM & other ) : BASE( other ) {}
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)
        
    protected:

    public:
        
        virtual Real Energy() override
        {
            return a * static_cast<Real>(2) * b;
        }
        
        virtual Real DEnergy() override
        {
            DX[0] += static_cast<Real>(2) * b;
            DY[0] += a * static_cast<Real>(2);
            
            return a * static_cast<Real>(2) * b;
        }

    public:
        
        virtual std::string Stats() const override
        {
            return ClassName();
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(AMB_DIM)+","+ToString(DEGREE)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
  
    };
    
} // namespace Repulsor

#undef BASE
#undef CLASS
