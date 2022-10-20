#pragma once

#define CLASS EnergyDimRestricted
#define BASE  EnergyBase<Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:
        using MeshBase_T       = typename BASE::MeshBase_T;

        using ValueContainer_T = typename BASE::ValueContainer_T;
        using Differential_T   = typename BASE::Differential_T;
        using TangentVector_T  = typename BASE::TangentVector_T;
        
        using Mesh_T           = SimplicialMesh<DOM_DIM,AMB_DIM,Real,Int,SReal,ExtReal>;

        
        explicit CLASS( ExtReal weight_ = static_cast<ExtReal>(1) )
        :   BASE(weight_)
        {}

        virtual ~CLASS() override = default;
        
        // Actual implementation to be specified by descendants.
        ExtReal compute( const MeshBase_T & M ) const override
        {
            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
                        
            if( Q != nullptr )
            {
                return value(*Q);
            }
            else
            {
                eprint(ClassName()+"::compute: Input could not be downcast to compatible type. Doing nothing.");
                return static_cast<ExtReal>(0);
            }
        }
        
        // Actual implementation to be specified by descendants.
        virtual ExtReal compute( const Mesh_T & M ) const = 0;
        
        // Do a down casy and delegate implementation further to descendant class.
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

        
        // Do a down casy and delegate implementation further to descendant class.
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

//
//        void SimplexEnergies( const MeshBase_T & M, ExtReal * output, bool addTo = false ) const override
//        {
//            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
//
//            if( Q != nullptr )
//            {
//                SimplexEnergies(*Q, output, addTo);
//            }
//            else
//            {
//                eprint(ClassName()+"::SimplexEnergies: Could not downcast input to compatible type. Doing nothing.");
//            }
//        }
        
//        virtual void SimplexEnergies( const Mesh_T & M, ExtReal * output, bool addTo = false  ) const = 0;
        
//        void Density( const MeshBase_T & M, ExtReal * output, bool addTo = false ) const override
//        {
//            const Mesh_T * Q = dynamic_cast<const Mesh_T *>(&M);
//
//            if( Q != nullptr )
//            {
//                Density(*Q, output, addTo);
//            }
//            else
//            {
//                eprint(ClassName()+"::Density: Could not downcast input to compatible type. Doing nothing.");
//            }
//        }
//
//        virtual void Density( const Mesh_T & M, ExtReal * output, bool addTo = false  ) const = 0;
//
    public:
        
//        virtual std::string Stats() const override = 0;
        
        static std::string className()
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
        virtual std::string ClassName() const override
        {
            return className();
        }
        
    };

}// namespace Repulsor

#undef BASE
#undef CLASS
