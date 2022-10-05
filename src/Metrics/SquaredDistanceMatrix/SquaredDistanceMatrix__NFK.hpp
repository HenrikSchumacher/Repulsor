#pragma once

#define CLASS SquaredDistanceMatrix__NFK
#define BASE  Metric__NFK<DOM_DIM1,DOM_DIM2,AMB_DIM,Real,Int,SReal>


namespace Repulsion
{
    
    template<int DOM_DIM1, int DOM_DIM2, int AMB_DIM, typename Real, typename Int, typename SReal>
    class CLASS : public BASE
    {
    protected:
        
        using BASE::S_near;
        using BASE::T_near;
        
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
        using BASE::NearDim;
        using BASE::CoordDim;
        using BASE::ProjectorDim;
        
    public:
        
        CLASS() : BASE() {}

        // Copy constructor
        CLASS( const CLASS & other ) : BASE ( other ) {}
        
        virtual ~CLASS() override = default;
        
        REPULSION__ADD_CLONE_CODE(CLASS)
        
    protected:
        
        Real * restrict values = nullptr;
                
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
        
        virtual void StartRow() override {}
        
        virtual void FinishRow( const Int diag_pos ) override
        {
            values[diag_pos] = static_cast<Real>(0);
        }
        
        virtual void Metric_Symmetric( const Int pos ) override
        {
            Real v [AMB_DIM] = {};
            
            Real r2 = static_cast<Real>(0);

            for( Int i = 0; i < AMB_DIM; ++i )
            {
                v[i] = y[i] - x[i];
                r2  += v[i] * v[i];
            }

            values[pos] = r2;
            
//            print("{i,j} = { "+ToString(S_ID)+","+ToString(T_ID)+" }");
//
//            print("x = { "+ToString(x[0])+","+ToString(x[1])+","+ToString(x[2])+" }");
//            print("y = { "+ToString(y[0])+","+ToString(y[1])+","+ToString(y[2])+" }");
//
//            print("r = "+ToString(r2));
        }
        
        
        virtual void Metric_Asymmetric( const Int pos ) override
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
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM1)+","+ToString(DOM_DIM1)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
  
    };
    
} // namespace Repulsion

#undef CLASS
#undef BASE
