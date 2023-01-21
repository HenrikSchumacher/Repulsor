#pragma once

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal>
    class EnergyDimRestricted : public EnergyBase<Real,Int,SReal,ExtReal>
    {
    private:
        
        using Base_T = EnergyBase<Real,Int,SReal,ExtReal>;
        
    public:
        
        using MeshBase_T        = typename Base_T::MeshBase_T;

        using Values_T          = typename Base_T::Values_T;
        using ValueContainer_T  = typename Base_T::ValueContainer_T;
        using TangentVector_T   = typename Base_T::TangentVector_T;
        using CotangentVector_T = typename Base_T::CotangentVector_T;
        
        using Mesh_T            = SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,SReal,ExtReal>;
        
        EnergyDimRestricted() = default;

        virtual ~EnergyDimRestricted() override = default;
        
        // Do a down cast and delegate implementation further to descendant class.
        ExtReal value( const MeshBase_T & M ) const override
        {
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                        
            if( Q != nullptr )
            {
                return value(*Q);
            }
            else
            {
                eprint(ClassName()+"::value: Input could not be downcast to compatible type. Doing nothing.");
                return static_cast<ExtReal>(0);
            }
        }
        
        // Actual implementation to be specified by descendants.
        virtual ExtReal value( const Mesh_T & M ) const = 0;

        
        // Do a down cast and delegate implementation further to descendant class.
        void differential( const MeshBase_T & M ) const override
        {
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                        
            if( Q != nullptr )
            {
                differential(*Q);
            }
            else
            {
                eprint(ClassName()+"::differential: Input could not be downcast to compatible type. Doing nothing.");
            }
        }
        
        // Actual implementation to be specified by descendants.
        virtual void differential( const Mesh_T & M ) const = 0;

    public:
        
        static std::string className()
        {
            return "EnergyDimRestricted<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        
    }; // class EnergyDimRestricted

}// namespace Repulsor
