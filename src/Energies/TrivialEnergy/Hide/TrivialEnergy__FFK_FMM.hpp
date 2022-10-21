#pragma once

#define CLASS TrivialEnergy__FFK_FMM
#define BASE Energy__FFK_FMM<AMB_DIM,DEGREE,Real,Int>

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
//        using BASE::x;
        using BASE::p;
        
        using BASE::b;
//        using BASE::y;
        using BASE::q;
        
    public:
        
        CLASS() : BASE() {}
        
        // Copy constructor
        CLASS( const CLASS & other ) : BASE( other ) {}
        
        virtual ~CLASS() override = default;

        __ADD_CLONE_CODE__(CLASS)
        
    public :
        
        virtual Real Energy() override
        {
            return static_cast<Real>(0);
        }
        
        virtual Real DEnergy() override
        {
            return static_cast<Real>(0);
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
