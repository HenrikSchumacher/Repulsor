#pragma once

#define CLASS TP_SingularMetric_FMM_Adaptive
#define BASE  Metric_FMM_Adaptive<DOM_DIM,AMB_DIM,DEGREE,Real,Int,SReal,ExtReal>

namespace Repulsor
{
    template<int DOM_DIM, int AMB_DIM, int DEGREE, typename Real, typename Int, typename SReal, typename ExtReal>
    class CLASS : public BASE
    {
        using BASE::bct;
        using BASE::near_values;
        using BASE::far_values;
        using BASE::metrics_initialized;
        using BASE::NearField;
        using BASE::FarField;
        
    public:
        
        static constexpr Int nnz_per_block = AMB_DIM+2;
        
        using BlockClusterTree_T = typename BASE::BlockClusterTree_T;

        CLASS() : BASE() {}
        
        CLASS(
            const BlockClusterTree_T & bct_,
            const Real alpha_,
            const Real beta_,
            const ExtReal weight = static_cast<ExtReal>(1),
            const AdaptivitySettings & settings_ = AdaptivitySettings()
        )
        :
        BASE(
            bct_,
            *std::make_unique<
                TP_SingularMetric__NFK          <DOM_DIM,AMB_DIM,Real,Int,SReal>
            > (alpha_,beta_),
            *std::make_unique<
                TP_SingularMetric__NFK_Adaptive <DOM_DIM,AMB_DIM,Real,Int,SReal>
            >(alpha_,beta_,settings_),
            *std::make_unique<
                TP_SingularMetric__FFK_FMM       <DOM_DIM,AMB_DIM,DEGREE,Real,Int>
            >(alpha_,beta_),
            weight
        )
        {
            this->RequireMetrics();
        }
        
        virtual ~CLASS() override = default;
        
    protected:
        
        void ComputeDiagonals() const override
        {
            ptic(ClassName()+"::ComputeDiagonals");

            eprint(ClassName()+"::ComputeDiagonals is yet to be implemented.");
            
//            Int cols = 1;
//
//            auto & T = bct->GetT();
//            auto & S = bct->GetS();
//
//            T.RequireBuffers(cols);
//            S.RequireBuffers(cols);
//
////            S.CleanseBuffers();
////            T.CleanseBuffers();
//
////            S.PrimitiveOutputBuffer().Fill( static_cast<Real>(0) );
//
//            {
//                      Real * restrict const in_buffer = T.PrimitiveInputBuffer().data();
//                const Real * restrict const near_data = T.PrimitiveNearFieldData().data();
//                const Int near_dim = T.NearDim();
//                const Int primitive_count = T.PrimitiveCount();
//
//                for( Int i = 0; i < primitive_count; ++i )
//                {
//                    in_buffer[i] = near_data[ near_dim * i ];
//                }
//            }
//
//            T.PrimitivesToClusters(false);
//
//            T.PercolateUp();
//
//            for( auto const & f : near_values)
//            {
//                ApplyNearFieldKernel( static_cast<Real>(1), f.first );
//
//                ApplyFarFieldKernel ( static_cast<Real>(1), f.first );
//
//                S.PercolateDown();
//
//                S.ClustersToPrimitives( true );
//
//                Real const * restrict const u = S.PrimitiveOutputBuffer().data();
//                Real const * restrict const a = S.PrimitiveNearFieldData().data();
//                Real       * restrict const values = near_values.find(f.first)->second.data();
//                Int  const * restrict const diag_ptr = bct->Near().Diag().data();
//
//                const Int m        = S.PrimitiveCount();
//                const Int near_dim = S.NearDim();
//
//                #pragma omp parallel for
//                for( Int i = 0; i < m; ++i )
//                {
//                    values[diag_ptr[i]] -= u[i]/a[ near_dim * i ];
//                }
//            }
            
            ptoc(ClassName()+"::ComputeDiagonals");
            
        } // ComputeDiagonals
        
           
    public:
        
        virtual void ApplyNearFieldKernel( const Real factor, const KernelType type ) const override
        {
            ptic(ClassName()+"::ApplyNearFieldKernel ("+KernelTypeName[type]+")" );
            
            eprint(ClassName()+"::ApplyNearFieldKernel is yet to be implemented.");
            
//            bct->Near().Multiply_DenseMatrix(
//                factor,
//                near_values.find(type)->second.data(),
//                bct->GetT().PrimitiveInputBuffer().data(),
//                static_cast<Real>(0),
//                bct->GetS().PrimitiveOutputBuffer().data(),
//                bct->GetT().BufferDimension()
//            );
            
            ptoc(ClassName()+"::ApplyNearFieldKernel ("+KernelTypeName[type]+")" );
        }
        
        virtual void ApplyFarFieldKernel( const Real factor, const KernelType type ) const override
        {
            ptic(ClassName()+"::ApplyFarFieldKernel ("+KernelTypeName[type]+")");
            
            wprint(ClassName()+"::ApplyFarFieldKernel is yet to be debugged.");
            
            // Caution: Dimensions are hard-coded!
            constexpr Int cols = (1 + AMB_DIM) * AMB_DIM;
            
            if( bct->GetT().BufferDimension() != cols )
            {
                eprint(ClassName()+"::ApplyFarFieldKernel: bct->GetT().BufferDimension() != "+ToString(cols));
                ptoc(ClassName()+"::ApplyFarFieldKernel ("+KernelTypeName[type]+")");
                return;
            }
            
            const Int m = bct->Far().RowCount();
//            const Int n = bct->Far().ColCount();
            const auto & job_ptr = bct->Far().JobPtr();
            
            Int  const * restrict const outer = bct->Far().Outer().data();
            Int  const * restrict const inner = bct->Far().Inner().data();
            Real const * restrict const A     = far_values.find(KernelType::MixedOrder)->second.data();
            
            Real       * restrict const U = bct->GetS().ClusterOutputBuffer().data();
            Real const * restrict const V = bct->GetT().ClusterInputBuffer().data();
            
            valprint("total nnz = ", far_values.find(KernelType::MixedOrder)->second.Size());
            
            if( factor == static_cast<Real>(0) )
            {
                zerofy_buffer( U, m * cols );
            }
            else
            {
                // The target buffer U may contain nan, so we have to _overwrite_ instead of multiply by 0 and add to it!
                
                #pragma omp parallel for num_threads( job_ptr.Size()-1 )
                for( Int thread = 0; thread < job_ptr.Size()-1; ++thread )
                {
                    const Int i_begin = job_ptr[thread  ];
                    const Int i_end   = job_ptr[thread+1];
                    
                    for( Int i = i_begin; i < i_end; ++i )
                    {
                        Real u_loc [AMB_DIM+1][AMB_DIM] = {{}};
                        
                        const Int idx_begin = outer[i  ];
                        const Int idx_end   = outer[i+1];
                        
                        for( Int idx = idx_begin; idx < idx_end; ++idx )
                        {
                            const Int  j = inner[idx];
                            
                            const Real K_ij = A[nnz_per_block * idx    ];
                            const Real K_ji = A[nnz_per_block * idx + 1];
                            
                            Real factor_ij [ AMB_DIM ];
                            Real factor_ji [ AMB_DIM ];
                            
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                const Real a = A[nnz_per_block * idx + 2 + k ];
                                factor_ij[k] = K_ij * a;
                                factor_ji[k] = K_ji * a;
                            }
                            
                            Real v_loc [AMB_DIM+1][AMB_DIM] = {{}};
                            
                            copy_buffer( &V[cols * j], &v_loc[0][0], cols );
                            
                            for( Int k = 0; k < AMB_DIM; ++k )
                            {
                                u_loc[0][k] -= (K_ij + K_ji) * v_loc[0][k];
                                
                                for( Int l = 0; l < AMB_DIM; ++l )
                                {
                                    u_loc[0] [k] += factor_ji[l] * v_loc[l+1][0];
                                    u_loc[k+1][l] += -factor_ij[l] * v_loc[0][k];
                                }
                            }
                        }
                        
                        for( Int k = 0; k < AMB_DIM+1; ++k )
                        {
                            for( Int l = 0; l < AMB_DIM; ++l )
                            {
                                u_loc[k][l] *= factor;
                            }
                        }
                        
                        copy_buffer( &u_loc[0][0], &U[cols * i], cols );
                    }
                }
            }
            
            ptoc(ClassName()+"::ApplyFarFieldKernel ("+KernelTypeName[type]+")");
        }
        
    public:
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+ToString(DEGREE)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+TypeName<SReal>::Get()+","+TypeName<ExtReal>::Get()+">";
        }
        
    };
    
}// namespace Repulsor

#undef BASE
#undef CLASS
