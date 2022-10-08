#pragma once

//##############################################################################################
//      TangentPointSingularMetric
//##############################################################################################

protected:

    void RequireTangentPointSingularMetric() const
    {
        if( !tpsm_initialized )
        {
            tpsm = std::make_unique<
                TP_SingularMetric_FMM_Adaptive<DOM_DIM,AMB_DIM,0,Real,Int,SReal,ExtReal>
            >( GetBlockClusterTree(), tp_alpha, tp_beta, tp_weight, adaptivity_settings );
                
            tpsm_initialized = true;
        }
    }
public:

    virtual void TangentPointSingularMetric_Multiply(
        const ExtReal alpha, const ExtReal * U,
        const ExtReal  beta,       ExtReal * V,
        Int cols
    ) const override
    {
        RequireTangentPointSingularMetric();
        
        tpsm->Multiply_DenseMatrix( static_cast<Real>(alpha), U, beta, V, cols, KernelType::MixedOrder );
    }


    virtual void TangentPointSingularMetric_Multiply(
        const ExtReal alpha, const Tensor1<ExtReal,Int> & U,
        const ExtReal  beta,       Tensor1<ExtReal,Int> & V
    ) const override
    {
        TangentPointSingularMetric_Multiply(alpha, U.data(), beta, V.data(), 1 );
    }

    virtual void TangentPointSingularMetric_Multiply(
        const ExtReal alpha, const Tensor2<ExtReal,Int> & U,
        const ExtReal  beta,       Tensor2<ExtReal,Int> & V
    ) const override
    {
        const Int cols = std::min( U.Dimension(1), V.Dimension(1) );
        TangentPointSingularMetric_Multiply(alpha, U.data(), beta, V.data(), cols );
    }

    virtual const Tensor1<Real,Int> & TangentPointSingularMetric_Values(
        const bool farQ
    ) const override
    {
        RequireTangentPointSingularMetric();
        
        if( farQ )
        {
            return tpsm->FarFieldValues().find(KernelType::MixedOrder)->second;
        }
        else
        {
            return tpsm->NearFieldValues().find(KernelType::MixedOrder)->second;
        }
    }


    virtual void TangentPointSingularMetric_ApplyKernel(
        const bool farQ
    ) const override
    {
        RequireTangentPointSingularMetric();

        if( farQ )
        {
            tpsm->ApplyFarFieldKernel ( static_cast<Real>(1), KernelType::MixedOrder );
        }
        else
        {
            tpsm->ApplyNearFieldKernel( static_cast<Real>(1), KernelType::MixedOrder );
        }
    }
