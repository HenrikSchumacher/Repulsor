#pragma once

namespace Repulsor
{
    template<typename Mesh_T> class EnergyDimRestricted {};

    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename LInt, typename SReal, typename ExtReal >
    class EnergyDimRestricted<SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>>
        : public EnergyBase<SimplicialMeshBase<Real,Int,LInt,SReal,ExtReal>>
    {
    public:
        
        using Mesh_T     = SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,LInt,SReal,ExtReal>;
        using MeshBase_T = typename Mesh_T::Base_T;
        
        using Base_T     = EnergyBase<MeshBase_T>;
        using Root_T     = EnergyBase<MeshBase_T>;

        using ValueContainer_T  = typename Base_T::ValueContainer_T;
        using TangentVector_T   = typename Base_T::TangentVector_T;
        using CotangentVector_T = typename Base_T::CotangentVector_T;
        
        EnergyDimRestricted() = default;

        virtual ~EnergyDimRestricted() override = default;
        
        using Base_T::Value;
        using Base_T::Differential;
        
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
        ExtReal differential( const MeshBase_T & M ) const override
        {
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                    
            if( Q != nullptr )
            {
                ExtReal en = differential(*Q);
                
                return en;
                
            }
            else
            {
                eprint(ClassName()+"::differential: Input could not be downcast to compatible type. Doing nothing.");
                return 0;
            }
        }
        
        // Actual implementation to be specified by descendants.
        virtual ExtReal differential( const Mesh_T & M ) const = 0;
        
        // Do a down cast and delegate implementation further to descendant class.
        ExtReal density( const MeshBase_T & M ) const override
        {
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                    
            if( Q != nullptr )
            {
                ExtReal en = density(*Q);
                
                return en;
                
            }
            else
            {
                eprint(ClassName()+"::density: Input could not be downcast to compatible type. Doing nothing.");
                return 0;
            }
        }
        
        // Actual implementation to be specified by descendants.
        virtual ExtReal density( const Mesh_T & M ) const = 0;

    public:
        
        static std::string className()
        {
            return "EnergyDimRestricted<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+","+TypeName<SReal>+","+TypeName<ExtReal>+">";
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        
    }; // class EnergyDimRestricted

}// namespace Repulsor
