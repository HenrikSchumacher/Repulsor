#pragma once

#define CLASS Energy_FMM
#define BASE  Energy_Restricted<DOM_DIM,AMB_DIM,Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, int DEGREE, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
    public:

        using Mesh_T     = typename BASE::Mesh_T;
        using N_Kernel_T = Energy__NFK     <DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>;
        using F_Kernel_T = Energy__FFK_FMM <AMB_DIM,DEGREE,Real,Int>;
        
        CLASS(
            const N_Kernel_T & N,
            const F_Kernel_T & F,
            ExtReal weight_ = static_cast<ExtReal>(1)
        )
        :   BASE(weight_)
        ,   N_ker(N.Clone())
        ,   F_ker(F.Clone())
        {}

        virtual ~CLASS() override = default;
        
        N_Kernel_T & NearFieldKernel()
        {
            return *N_ker;
        }
        
        const N_Kernel_T & NearFieldKernel() const
        {
            return *N_ker;
        }
        
        F_Kernel_T & FarFieldKernel()
        {
            return *F_ker;
        }
        
        const F_Kernel_T & FarFieldKernel() const
        {
            return *F_ker;
        }
        
        // Returns the current value of the energy.
        ExtReal Value( const Mesh_T & M ) const override
        {
            ptic(ClassName()+"::Value");
            
            DUMP(Stats());
            
            const Real near = NearField(M);
//            const Real near = 0;
            const Real  far =  FarField(M);
//            const Real  far = 0;
            
            DUMP(near);
            DUMP(far);
            
            ptoc(ClassName()+"::Value");
            return weight * static_cast<ExtReal>(near+far);
        }
        
        // Returns the current differential of the energy, stored in the given
        // V x AMB_DIM matrix, where each row holds the differential (a AMB_DIM-vector) with
        // respect to the corresponding vertex.
        // Returns the current value of the energy.
        ExtReal Differential( const Mesh_T & M, ExtReal * output, bool addTo = false ) const override
        {
            ptic(ClassName()+"::Differential");
            
            DUMP(Stats());
            
            M.GetClusterTree().CleanseDerivativeBuffers();
            
            const Real near = DNearField(M);
//            const Real near = 0;
            const Real  far =  DFarField(M);
//            const Real far = 0;
            DUMP(near);
            DUMP(far);
            
            M.Assemble_ClusterTree_Derivatives(output, weight, addTo);

            
            ptoc(ClassName()+"::Differential");
            return weight * (near+far);
        }
        
        virtual void SimplexEnergies( const Mesh_T & M, ExtReal * output, bool addTo = false  ) const override
        {
            ptic(ClassName()+"::SimplexEnergies");
            
            DUMP(Stats());
            
            M.GetClusterTree().CleanseDerivativeBuffers();
            
            NearField_SimplexEnergies(M);
             FarField_SimplexEnergies(M);
            
            M.Assemble_ClusterTree_SimplexEnergies(output, weight, addTo);
            
            ptoc(ClassName()+"::SimplexEnergies");
        }
        
        virtual void Density( const Mesh_T & M, ExtReal * output, bool addTo = false  ) const override
        {
            ptic(ClassName()+"::Density");
            
            DUMP(Stats());
            
            M.GetClusterTree().CleanseDerivativeBuffers();
            
            NearField_SimplexEnergies(M);
             FarField_SimplexEnergies(M);
            
            M.Assemble_ClusterTree_Density(output, weight, addTo);
            
            ptoc(ClassName()+"::Density");
        }
        
    public:
        
        bool use_blocking = false;
        
    protected:
        
        using BASE::weight;
        
        mutable std::unique_ptr<N_Kernel_T> N_ker;
        mutable std::unique_ptr<F_Kernel_T> F_ker;
        
       virtual Real NearField( const Mesh_T & M ) const
        {
            ptic(ClassName()+"::NearField");

            auto & bct  = M.GetBlockClusterTree();
            auto & near = bct.Near();
            
            if( near.NonzeroCount() <= 0 )
            {
                wprint(ClassName()+"::NearField: no near field found. Returning 0." );
                ptoc(ClassName()+"::NearField");
                return static_cast<Real>(0);
            }

            Real sum = static_cast<Real>(0);
            
            const auto & job_ptr = near.UpperTriangularJobPtr();
            
            // Make sure that near.Diag() is created before the parallel region.
            (void)near.Diag();
            
            const Int thread_count = job_ptr.Size()-1;
            
            #pragma omp parallel for num_threads(thread_count) reduction( + : sum )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                Int const * restrict const diag  = near.Diag().data();
                Int const * restrict const outer = near.Outer().data();
                Int const * restrict const inner = near.Inner().data();
                
                std::unique_ptr<N_Kernel_T> N = N_ker->Clone();
                
                // This looks massive, but is only exchange of pointers.
                (void)N->LoadNearField(
                    bct.GetS().PrimitiveNearFieldData(),
                    bct.GetT().PrimitiveNearFieldData()
                );
//                (void)N->LoadDNearField(
//                    bct.GetS().ThreadPrimitiveDNearFieldData()[thread],
//                    bct.GetT().ThreadPrimitiveDNearFieldData()[thread]
//                );
                (void)N->LoadPrimitiveSerializedData(
                    bct.GetS().PrimitiveSerializedData(),
                    bct.GetT().PrimitiveSerializedData()
                );
                
                Real local_sum = static_cast<Real>(0);

                const Int i_begin = job_ptr[thread  ];
                const Int i_end   = job_ptr[thread+1];

                for( Int i = i_begin; i < i_end; ++i )
                {
                    N->LoadS(i);
                    
                    const Int k_begin = diag [i]+1;
                    const Int k_end   = outer[i+1];
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        const Int j = inner[k];
                        
                        N->LoadT(j);
                    
                        local_sum += N->Energy();
                        
                    } // for( Int k = k_begin; k < k_end; ++k )
                
                } // for( Int i = i_begin; i < i_end; ++i )
                
                sum += local_sum;
    
            } // #pragma omp parallel
            
            ptoc(ClassName()+"::NearField");
            return sum;
        }
        
        virtual Real DNearField( const Mesh_T & M ) const
        {
            ptic(ClassName()+"::DNearField");

            auto & bct  = M.GetBlockClusterTree();
            auto & near = bct.Near();
            
            if( near.NonzeroCount() <= 0 )
            {
                wprint(ClassName()+"::DNearField: no near field found. Returning 0." );
                ptoc(ClassName()+"::DNearField");
                return static_cast<Real>(0);
            }

            Real sum = static_cast<Real>(0);
            
            const auto & job_ptr = near.UpperTriangularJobPtr();
            
            // Make sure that near.Diag() is created before the parallel region.
            (void)near.Diag();
            
            const Int thread_count = job_ptr.Size()-1;
            
            #pragma omp parallel for num_threads(thread_count) reduction( + : sum )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                Int const * restrict const diag  = near.Diag().data();
                Int const * restrict const outer = near.Outer().data();
                Int const * restrict const inner = near.Inner().data();
                
                std::unique_ptr<N_Kernel_T> N = N_ker->Clone();
                
                // This looks massive, but is only exchange of pointers.
                (void)N->LoadNearField(
                    bct.GetS().PrimitiveNearFieldData(),
                    bct.GetT().PrimitiveNearFieldData()
                );
                (void)N->LoadDNearField(
                    bct.GetS().ThreadPrimitiveDNearFieldData()[thread],
                    bct.GetT().ThreadPrimitiveDNearFieldData()[thread]
                );
                (void)N->LoadPrimitiveSerializedData(
                    bct.GetS().PrimitiveSerializedData(),
                    bct.GetT().PrimitiveSerializedData()
                );
                
                Real local_sum = static_cast<Real>(0);

                const Int i_begin = job_ptr[thread  ];
                const Int i_end   = job_ptr[thread+1];

                for( Int i = i_begin; i < i_end; ++i )
                {
                    N->LoadS(i);
                    N->CleanseDBufferS();
                    
                    const Int k_begin = diag [i]+1;
                    const Int k_end   = outer[i+1];
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        const Int j = inner[k];
                        
                        N->LoadT(j);
                        N->CleanseDBufferT();
                        
                        local_sum += N->DEnergy();
                        
                        N->WriteDBufferT();
                        
                    } // for( Int k = k_begin; k < k_end; ++k )
                    
                    N->WriteDBufferS();
                
                } // for( Int i = i_begin; i < i_end; ++i )
                
                sum += local_sum;
                
            } // #pragma omp parallel
            
            ptoc(ClassName()+"::DNearField");
            return sum;
        }
        
        virtual void NearField_SimplexEnergies( const Mesh_T & M ) const
        {
            ptic(ClassName()+"::NearField_SimplexEnergies");
            
            auto & bct  = M.GetBlockClusterTree();
            auto & near = bct.Near();
            
            if( near.NonzeroCount() <= 0 )
            {
                wprint(ClassName()+"::NearField_SimplexEnergies: no near field found. Returning 0." );
                ptoc(ClassName()+"::NearField_SimplexEnergies");
                return;
            }
            
            const auto & job_ptr = near.UpperTriangularJobPtr();
            
            // Make sure that near.Diag() is created before the parallel region.
            (void)near.Diag();
            
            const Int thread_count = job_ptr.Size()-1;
            
            #pragma omp parallel for num_threads(thread_count)
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                Int const * restrict const diag  = near.Diag().data();
                Int const * restrict const outer = near.Outer().data();
                Int const * restrict const inner = near.Inner().data();
                
                std::unique_ptr<N_Kernel_T> N = N_ker->Clone();
                
                // This looks massive, but is only exchange of pointers.
                (void)N->LoadNearField(
                    bct.GetS().PrimitiveNearFieldData(),
                    bct.GetT().PrimitiveNearFieldData()
                );
                //                (void)N->LoadDNearField(
                //                    bct.GetS().ThreadPrimitiveDNearFieldData()[thread],
                //                    bct.GetT().ThreadPrimitiveDNearFieldData()[thread]
                //                );
                (void)N->LoadPrimitiveSerializedData(
                    bct.GetS().PrimitiveSerializedData(),
                    bct.GetT().PrimitiveSerializedData()
                );
                
                auto & X = bct.GetS().ThreadPrimitiveDNearFieldData()[thread];
                auto & Y = bct.GetT().ThreadPrimitiveDNearFieldData()[thread];
                
                const Int i_begin = job_ptr[thread  ];
                const Int i_end   = job_ptr[thread+1];
                
                for( Int i = i_begin; i < i_end; ++i )
                {
                    N->LoadS(i);
                    
                    const Int k_begin = diag [i]+1;
                    const Int k_end   = outer[i+1];
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        const Int j = inner[k];
                        
                        N->LoadT(j);
                        
                        const Real En = N->Energy();
                        
                        X(i,0) += En;
                        Y(j,0) += En;
                        
                    } // for( Int k = k_begin; k < k_end; ++k )
                    
                } // for( Int i = i_begin; i < i_end; ++i )
                
            } // #pragma omp parallel
            
            ptoc(ClassName()+"::NearField_SimplexEnergies");
        }
        
        
        
        
        
        Real FarField( const Mesh_T & M ) const
        {
            ptic(ClassName()+"::FarField");

            auto & bct = M.GetBlockClusterTree();
            auto & far = bct.Far();

            if( far.NonzeroCount() <= 0 )
            {
//                wprint(ClassName()+"::FarField: no admissible blocks found. Returning 0." );
                ptoc(ClassName()+"::FarField");
                return static_cast<Real>(0);
            }
          
            Real sum = static_cast<Real>(0);
            
            const auto & job_ptr = far.UpperTriangularJobPtr();
            
            // Make sure that far.Diag() is created before the parallel region.
            (void)far.Diag();
            
            const Int thread_count = job_ptr.Size()-1;
            
            #pragma omp parallel for num_threads(thread_count) reduction( + : sum )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                Int const * restrict const diag  = far.Diag().data();
                Int const * restrict const outer = far.Outer().data();
                Int const * restrict const inner = far.Inner().data();
                
                std::unique_ptr<F_Kernel_T> F = F_ker->Clone();
                
                (void)F->LoadFarField(
                    bct.GetS().ClusterFarFieldData(),
                    bct.GetT().ClusterFarFieldData()
                );
//                (void)F->LoadDFarField(
//                    bct.GetS().ThreadClusterDFarFieldData()[thread],
//                    bct.GetT().ThreadClusterDFarFieldData()[thread]
//                );
                
                Real local_sum = static_cast<Real>(0);

                const Int i_begin = job_ptr[thread  ];
                const Int i_end   = job_ptr[thread+1];
                
                for( Int i = i_begin; i < i_end; ++i )
                {
                    F->LoadS(i);
                    
                    const Int k_begin = diag[i];
                    const Int k_end   = outer[i+1];
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        const Int j = inner[k];

                        F->LoadT(j);
                    
                        local_sum += F->Energy();

                    } // for( Int k = k_begin; k < k_end; ++k )

                } // for( Int i = i_begin; i < i_end; ++i )

                sum += local_sum;
                
            } // #pragma omp parallel
            
            ptoc(ClassName()+"::FarField");
            return sum;
        }
        
        Real DFarField( const Mesh_T & M ) const
        {
            ptic(ClassName()+"::DFarField");

            auto & bct = M.GetBlockClusterTree();
            auto & far = bct.Far();

            if( far.NonzeroCount() <= 0 )
            {
//                wprint(ClassName()+"::DFarField: no admissible blocks found. Returning 0." );
                ptoc(ClassName()+"::DFarField");
                return static_cast<Real>(0);
            }
          
            Real sum = static_cast<Real>(0);
            
            const auto & job_ptr = far.UpperTriangularJobPtr();
            
            // Make sure that far.Diag() is created before the parallel region.
            (void)far.Diag();
            
            const Int thread_count = job_ptr.Size()-1;
            
            #pragma omp parallel for num_threads(thread_count) reduction( + : sum )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                Int const * restrict const diag  = far.Diag().data();
                Int const * restrict const outer = far.Outer().data();
                Int const * restrict const inner = far.Inner().data();
                
                std::unique_ptr<F_Kernel_T> F = F_ker->Clone();
                
                (void)F->LoadFarField(
                    bct.GetS().ClusterFarFieldData(),
                    bct.GetT().ClusterFarFieldData()
                );
                (void)F->LoadDFarField(
                    bct.GetS().ThreadClusterDFarFieldData()[thread],
                    bct.GetT().ThreadClusterDFarFieldData()[thread]
                );
                
                Real local_sum = static_cast<Real>(0);

                const Int i_begin = job_ptr[thread  ];
                const Int i_end   = job_ptr[thread+1];
                
                for( Int i = i_begin; i < i_end; ++i )
                {
                    F->LoadS(i);
                    F->CleanseDBufferS();

                    const Int k_begin = diag[i];
                    const Int k_end   = outer[i+1];
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        const Int j = inner[k];

                        F->LoadT(j);
                        F->CleanseDBufferT();
                        
                        local_sum += F->DEnergy();

                        F->WriteDBufferT();
                        
                    } // for( Int k = k_begin; k < k_end; ++k )

                    F->WriteDBufferS();
                    
                } // for( Int i = i_begin; i < i_end; ++i )

                sum += local_sum;
                
            } // #pragma omp parallel
            
            ptoc(ClassName()+"::DFarField");
            return sum;
        }
        
        void FarField_SimplexEnergies( const Mesh_T & M ) const
        {
            ptic(ClassName()+"::FarField_SimplexEnergies");

            auto & bct = M.GetBlockClusterTree();
            auto & far = bct.Far();

            if( far.NonzeroCount() <= 0 )
            {
//                wprint(ClassName()+"::FarField_SimplexEnergies: no admissible blocks found. Returning 0." );
                ptoc(ClassName()+"::FarField_SimplexEnergies");
                return;
            }

            const auto & job_ptr = far.UpperTriangularJobPtr();
            
            // Make sure that far.Diag() is created before the parallel region.
            (void)far.Diag();
            
            const Int thread_count = job_ptr.Size()-1;
            
            #pragma omp parallel for num_threads(thread_count)
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                Int const * restrict const diag  = far.Diag().data();
                Int const * restrict const outer = far.Outer().data();
                Int const * restrict const inner = far.Inner().data();
                
                std::unique_ptr<F_Kernel_T> F = F_ker->Clone();
                
                (void)F->LoadFarField(
                    bct.GetS().ClusterFarFieldData(),
                    bct.GetT().ClusterFarFieldData()
                );
//                (void)F->LoadDFarField(
//                    bct.GetS().ThreadClusterDFarFieldData()[thread],
//                    bct.GetT().ThreadClusterDFarFieldData()[thread]
//                );
                
                auto & X = bct.GetS().ThreadClusterDFarFieldData()[thread];
                auto & Y = bct.GetT().ThreadClusterDFarFieldData()[thread];
                
                const Int i_begin = job_ptr[thread  ];
                const Int i_end   = job_ptr[thread+1];
                
                for( Int i = i_begin; i < i_end; ++i )
                {
                    F->LoadS(i);
                    
                    const Int k_begin = diag[i];
                    const Int k_end   = outer[i+1];
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        const Int j = inner[k];

                        F->LoadT(j);
                    
                        const Real En = F->Energy();
                        
                        X(i,0) += En;
                        Y(j,0) += En;

                    } // for( Int k = k_begin; k < k_end; ++k )

                } // for( Int i = i_begin; i < i_end; ++i )
                
            } // #pragma omp parallel
            
            ptoc(ClassName()+"::FarField_SimplexEnergies");
        }
        
    public:
        
        virtual std::string Stats() const override
        {
            std::stringstream s;
            
            s << ClassName() << ": \n\n";
            
            if( N_ker != nullptr )
            {
                s << "          NearFieldKernel = " << N_ker->Stats() << "\n\n";
            }
            
            if( F_ker != nullptr )
            {
                s << "          FarFieldKernel  = " << F_ker->Stats() << "\n";
            }
            return s.str();
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+ToString(DEGREE)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    };
    
}// namespace Repulsor

#undef BASE
#undef CLASS
