#pragma once

#define CLASS Area__NFK_Adaptive
#define BASE  Energy__NFK_Adaptive<DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Int,SReal>

namespace Repulsor
{
    template<int DOM_DIM1, int DOM_DIM2, int AMB_DIM, typename Real, typename Int, typename SReal>
    class CLASS : public BASE
    {
    protected:
        
        using BASE::S_near;
        using BASE::S_D_near;
        
        using BASE::T_near;
        using BASE::T_D_near;
        
        using BASE::DX;
        using BASE::DY;
        
        using BASE::a;
        using BASE::x;
        using BASE::x_buffer;
        using BASE::p;
        
        using BASE::b;
        using BASE::y;
        using BASE::y_buffer;
        using BASE::q;
        
        using BASE::theta2;
        using BASE::S;
        using BASE::T;
        
        using BASE::S_ID;
        using BASE::T_ID;
        
        using BASE::tri_i;
        using BASE::tri_j;
        using BASE::lin_k;
        
    public:
        
        explicit CLASS( const AdaptivitySettings & settings_ = AdaptivitySettings() )
        :
            BASE(settings_)
        {}

        // Copy constructor
        CLASS( const CLASS & other ) : BASE(other) {}
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)
        
    protected:
    
    public:
        
        virtual Real energy() const override
        {
            const Real w = S.Weight() * T.Weight();
         
            return a * static_cast<Real>(2) * b * w;
        }
        
        virtual Real denergy() const override
        {
            const Real w = S.Weight() * T.Weight();
         
            DX[0] += static_cast<Real>(2) * b * w;
            DY[0] += a * static_cast<Real>(2) * w;
            
            return a * static_cast<Real>(2) * b * w;
        }
        
    public:
        
        virtual bool IsRepulsive() const override
        {
            return false;
        }
        
        virtual std::string Stats() const override
        {
            return ClassName()+": theta = "+ToString(std::sqrt(theta2));
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM1)+","+ToString(DOM_DIM2)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+">";
        }
  
    };

} // namespace Repulsor

#undef CLASS
#undef BASE
