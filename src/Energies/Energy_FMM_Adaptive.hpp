#pragma once

#define CLASS Energy_FMM_Adaptive
#define BASE  Energy_FMM<DOM_DIM,AMB_DIM,DEGREE,Real,Int,SReal,ExtReal>

namespace Repulsion
{
    template<int DOM_DIM, int AMB_DIM, int DEGREE, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
        
        using Mesh_T     = typename BASE::Mesh_T;
        using N_Kernel_T = typename BASE::N_Kernel_T;
        using F_Kernel_T = typename BASE::F_Kernel_T;
        using A_Kernel_T = Energy__NFK_Adaptive<DOM_DIM,DOM_DIM,AMB_DIM,Real,Int,SReal>;
        
    protected:
        
        using BASE::use_blocking;
        using BASE::N_ker;
        using BASE::F_ker;
        using BASE::weight;
        
        mutable std::unique_ptr<A_Kernel_T> A_ker;
        
    public:

        CLASS(
            const N_Kernel_T & N,
            const A_Kernel_T & A,
            const F_Kernel_T & F,
            ExtReal weight_ = static_cast<ExtReal>(1)
        )
        :   BASE(N,F,weight_)
        ,   A_ker( A.Clone() )
        {}

        virtual ~CLASS() override = default;

            
        virtual std::string Stats() const override
        {
            std::stringstream s;
            
            s << ClassName() <<": \n\n";
            if( N_ker != nullptr )
            {
                s <<  "          NearFieldKernel          = " << N_ker->Stats() << "\n\n";
            }
            if( A_ker != nullptr )
            {
                s <<  "          NearFieldKernel_Adaptive = " << A_ker->Stats() << "\n\n";
            }
            if( F_ker != nullptr )
            {
                s <<  "          FarFieldKernel           = " << F_ker->Stats() << "\n";
            }
            return s.str();
        }
        
        A_Kernel_T & NearFieldKernelAdaptive()
        {
            return *A_ker;
        }
        
        const A_Kernel_T & NearFieldKernelAdaptive() const
        {
            return *A_ker;
        }

        
    public:
        
    protected:
        
        Real NearField( const Mesh_T & M ) const override
        {
            Real sum = NearField_Adaptive(M) + NearField_Nonadaptive(M);
            return sum;
        }
        
        
        Real DNearField( const Mesh_T & M ) const override
        {
            Real sum = DNearField_Adaptive(M) + DNearField_Nonadaptive(M);
            return sum;
        }
        
        void NearField_SimplexEnergies( const Mesh_T & M ) const override
        {
            NearField_SimplexEnergies_Adaptive(M);

            NearField_SimplexEnergies_Nonadaptive(M);
        }
        
        
        Real NearField_Nonadaptive( const Mesh_T & M ) const
        {
            ptic(ClassName()+"::NearField_Nonadaptive");
            
            auto & bct  = M.GetBlockClusterTree();
            auto & near = bct.AdaptiveNoSubdivisionData();
            
            if( near.NonzeroCount() <= 0 )
            {
//                wprint(ClassName()+"::NearField_Nonadaptive: no near field found. Returning 0." );
                logprint("No near field found. Returning 0.");
                ptoc(ClassName()+"::NearField_Nonadaptive");
                return static_cast<Real>(0);
            }

            Real sum = static_cast<Real>(0);
            
            const auto & job_ptr = near.UpperTriangularJobPtr();
            
            // Make sure that near.Diag() is created before the parallel region.
            (void)near.Diag();
            
            const Int thread_count = job_ptr.Size()-1;
            
//            print( "near.Diag()= " + near.Diag().ToString() );
            
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
                    
                    const Int k_begin = diag[i]+1;
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

            ptoc(ClassName()+"::NearField_Nonadaptive");
            return sum;
        }
        
        
        Real DNearField_Nonadaptive( const Mesh_T & M ) const
        {
            ptic(ClassName()+"::DNearField_Nonadaptive");

            auto & bct  = M.GetBlockClusterTree();
            auto & near = bct.AdaptiveNoSubdivisionData();
            
            if( near.NonzeroCount() <= 0 )
            {
                logprint("No near field found. Returning 0.");
                ptoc(ClassName()+"::DNearField_Nonadaptive");
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
                    }
                    
                    N->WriteDBufferS();
                
                } // for( Int i = i_begin; i < i_end; ++i )
                
                sum += local_sum;
                
            } // #pragma omp parallel
            
            ptoc(ClassName()+"::DNearField_Nonadaptive");
            return sum;
        }
        
        
        void NearField_SimplexEnergies_Nonadaptive( const Mesh_T & M ) const
        {
            ptic(ClassName()+"::NearField_SimplexEnergies_Nonadaptive");

            auto & bct  = M.GetBlockClusterTree();
            auto & near = bct.AdaptiveNoSubdivisionData();
            
            if( near.NonzeroCount() <= 0 )
            {
//                wprint(ClassName()+"::NearField_Nonadaptive: no near field found. Returning 0." );
                logprint("No near field found. Returning 0.");
                ptoc(ClassName()+"::NearField_SimplexEnergies_Nonadaptive");
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
            
            ptoc(ClassName()+"::NearField_SimplexEnergies_Nonadaptive");
        }
        
        
        Real NearField_Adaptive( const Mesh_T & M ) const
        {
            ptic(ClassName()+"::NearField_Adaptive");

            auto & bct  = M.GetBlockClusterTree();
            auto & near = bct.AdaptiveSubdivisionData();

            
            if( near.NonzeroCount() <= 0 )
            {
                logprint("No primitive pairs to subdivide found. Returning 0.");
                ptoc(ClassName()+"::NearField_Adaptive");
                return static_cast<Real>(0);
            }
            
            if( A_ker->IsRepulsive() )
            {
                if( bct.PrimitiveIntersectionCount() > 0 )
                {
                    ptoc(ClassName()+"::NearField_Adaptive");
                    return std::numeric_limits<Real>::max();
                }
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
                
                std::unique_ptr<A_Kernel_T> A = A_ker->Clone();
                
#ifdef REPULSION__PRINT_REPORTS_FOR_ADAPTIVE_KERNELS
                A->CreateLogFile();
#endif
                
                // This looks massive, but is only exchange of pointers.
                (void)A->LoadNearField(
                    bct.GetS().PrimitiveNearFieldData(),
                    bct.GetT().PrimitiveNearFieldData()
                );
//                (void)N->LoadDNearField(
//                    bct.GetS().ThreadPrimitiveDNearFieldData()[thread],
//                    bct.GetT().ThreadPrimitiveDNearFieldData()[thread]
//                );
                (void)A->LoadPrimitiveSerializedData(
                    bct.GetS().PrimitiveSerializedData(),
                    bct.GetT().PrimitiveSerializedData()
                );
                
                Real local_sum = static_cast<Real>(0);

                const Int i_begin = job_ptr[thread  ];
                const Int i_end   = job_ptr[thread+1];

                for( Int i = i_begin; i < i_end; ++i )
                {
                    A->LoadS(i);
                    
                    const Int k_begin = diag [i];
                    const Int k_end   = outer[i+1];
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        const Int j = inner[k];
                        
                        A->LoadT(j);
                    
                        local_sum += A->Energy();
                        
                    } // for( Int k = k_begin; k < k_end; ++k )
                
                } // for( Int i = i_begin; i < i_end; ++i )
                
                sum += local_sum;
    
            } // #pragma omp parallel
            
            ptoc(ClassName()+"::NearField_Adaptive");
            return sum;
        }
        
        
        Real DNearField_Adaptive( const Mesh_T & M ) const
        {
            ptic(ClassName()+"::DNearField_Adaptive");

            auto & bct  = M.GetBlockClusterTree();
            auto & near = bct.AdaptiveSubdivisionData();
            
            if( near.NonzeroCount() <= 0 )
            {
                logprint("No primitive pairs to subdivide found. Returning 0.");
                ptoc(ClassName()+"::DNearField_Adaptive");
                return static_cast<Real>(0);
            }
            
            if( A_ker->IsRepulsive() )
            {
                if( bct.PrimitiveIntersectionCount() > 0 )
                {
                    ptoc(ClassName()+"::NearField_Adaptive");
                    return std::numeric_limits<Real>::max();
                }
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
                
                std::unique_ptr<A_Kernel_T> A = A_ker->Clone();
                
                // This looks massive, but is only exchange of pointers.
                (void)A->LoadNearField(
                    bct.GetS().PrimitiveNearFieldData(),
                    bct.GetT().PrimitiveNearFieldData()
                );
                (void)A->LoadDNearField(
                    bct.GetS().ThreadPrimitiveDNearFieldData()[thread],
                    bct.GetT().ThreadPrimitiveDNearFieldData()[thread]
                );
                (void)A->LoadPrimitiveSerializedData(
                    bct.GetS().PrimitiveSerializedData(),
                    bct.GetT().PrimitiveSerializedData()
                );
                
                Real local_sum = static_cast<Real>(0);

                const Int i_begin = job_ptr[thread  ];
                const Int i_end   = job_ptr[thread+1];

                for( Int i = i_begin; i < i_end; ++i )
                {
                    A->LoadS(i);
                    A->CleanseDBufferS();
                    
                    const Int k_begin = diag [i];
                    const Int k_end   = outer[i+1];
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        const Int j = inner[k];
                        
                        A->LoadT(j);
                        A->CleanseDBufferT();
                        
                        local_sum += A->DEnergy();
                        
                        A->WriteDBufferT();
                    }
                    
                    A->WriteDBufferS();
                
                } // for( Int i = i_begin; i < i_end; ++i )
                
                sum += local_sum;
                
            } // #pragma omp parallel
            
            ptoc(ClassName()+"::DNearField_Adaptive");
            return sum;
        }
        
        
        void NearField_SimplexEnergies_Adaptive( const Mesh_T & M ) const
        {
            ptic(ClassName()+"::NearField_SimplexEnergies_Adaptive");

            auto & bct  = M.GetBlockClusterTree();
            auto & near = bct.AdaptiveSubdivisionData();

            if( near.NonzeroCount() <= 0 )
            {
                logprint("No primitive pairs to subdivide found. Returning 0.");
                ptoc(ClassName()+"::NearField_SimplexEnergies_Adaptive");
                return;
            }
            
//            if( A_ker->IsRepulsive() )
//            {
//                if( bct.CollisionCount() > 0 )
//                {
//                    ptoc(ClassName()+"::NearField_SimplexEnergies_Adaptive");
//                    return;
//                }
//            }
            
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
                
                std::unique_ptr<A_Kernel_T> A = A_ker->Clone();
                
                // This looks massive, but is only exchange of pointers.
                (void)A->LoadNearField(
                    bct.GetS().PrimitiveNearFieldData(),
                    bct.GetT().PrimitiveNearFieldData()
                );
//                (void)N->LoadDNearField(
//                    bct.GetS().ThreadPrimitiveDNearFieldData()[thread],
//                    bct.GetT().ThreadPrimitiveDNearFieldData()[thread]
//                );
                (void)A->LoadPrimitiveSerializedData(
                    bct.GetS().PrimitiveSerializedData(),
                    bct.GetT().PrimitiveSerializedData()
                );
                
                auto & X = bct.GetS().ThreadPrimitiveDNearFieldData()[thread];
                auto & Y = bct.GetT().ThreadPrimitiveDNearFieldData()[thread];
                
                const Int i_begin = job_ptr[thread  ];
                const Int i_end   = job_ptr[thread+1];

                for( Int i = i_begin; i < i_end; ++i )
                {
                    A->LoadS(i);
                    
                    const Int k_begin = diag [i];
                    const Int k_end   = outer[i+1];
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        const Int j = inner[k];
                        
                        A->LoadT(j);
                    
                        const Real En = A->Energy();
                        
                        X(i,0) += En;
                        Y(j,0) += En;
                        
                    } // for( Int k = k_begin; k < k_end; ++k )
                
                } // for( Int i = i_begin; i < i_end; ++i )
    
            } // #pragma omp parallel
            
            ptoc(ClassName()+"::NearField_SimplexEnergies_Adaptive");
        }
        
    public:
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+ToString(DEGREE)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    };
    
}// namespace Repulsion

#undef BASE
#undef CLASS
