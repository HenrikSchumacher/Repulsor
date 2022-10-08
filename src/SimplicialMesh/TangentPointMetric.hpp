#pragma once

//##############################################################################################
//      TangentPointMetric
//##############################################################################################
       
protected:

    void RequireTangentPointMetric() const
    {
        if( !tpm_initialized )
        {
            tpm = std::make_unique<
                TP_Metric_FMM_Adaptive<DOM_DIM,AMB_DIM,0,Real,Int,SReal,ExtReal>
            >( GetBlockClusterTree(), tp_alpha, tp_beta, tp_weight, adaptivity_settings );
                
            tpm_initialized = true;
        }
    }

public:

    virtual void TangentPointMetric_Multiply(
        const ExtReal alpha, const ExtReal * U,
        const ExtReal  beta,       ExtReal * V,
        Int cols,
        KernelType kernel
    ) const override
    {
        RequireTangentPointMetric();
        
        tpm->Multiply_DenseMatrix( static_cast<Real>(alpha), U, beta, V, cols, kernel );
    }


    virtual void TangentPointMetric_Multiply(
        const ExtReal alpha, const Tensor1<ExtReal,Int> & U,
        const ExtReal  beta,       Tensor1<ExtReal,Int> & V,
        KernelType kernel
    ) const override
    {
        TangentPointMetric_Multiply(alpha, U.data(), beta, V.data(), 1, kernel);
    }

    virtual void TangentPointMetric_Multiply(
        const ExtReal alpha, const Tensor2<ExtReal,Int> & U,
        const ExtReal  beta,       Tensor2<ExtReal,Int> & V,
        KernelType kernel
    ) const override
    {
        const Int cols = std::min( U.Dimension(1), V.Dimension(1) );
        TangentPointMetric_Multiply(alpha, U.data(), beta, V.data(), cols, kernel);
    }

    virtual const Tensor1<Real,Int> & TangentPointMetric_Values(
        const bool farQ,
        const KernelType kernel
    ) const override
    {
        RequireTangentPointMetric();
        
        if( farQ )
        {
            return tpm->FarFieldValues().find(kernel)->second;
        }
        else
        {
            return tpm->NearFieldValues().find(kernel)->second;
        }
    }

    
    virtual void TangentPointMetric_Multiply(
        const ExtReal alpha, const ExtReal * U,
        const ExtReal  beta,       ExtReal * V,
        Int cols
    ) const override
    {
        RequireTangentPointMetric();
        
        tpm->Multiply_DenseMatrix(
            static_cast<Real>(alpha), U, beta, V, cols, KernelType::HighOrder );
        
        tpm->Multiply_DenseMatrix(
            static_cast<Real>(alpha), U, static_cast<Real>(1), V, cols, KernelType::LowOrder );
    }
        
    virtual void TangentPointMetric_Multiply(
        const ExtReal alpha, const Tensor1<ExtReal,Int> & U,
        const ExtReal  beta,       Tensor1<ExtReal,Int> & V
    ) const override
    {
        TangentPointMetric_Multiply(alpha, U.data(), beta, V.data(), 1);
    }

    virtual void TangentPointMetric_Multiply(
        const ExtReal alpha, const Tensor2<ExtReal,Int> & U,
        const ExtReal  beta,       Tensor2<ExtReal,Int> & V
    ) const override
    {
        const Int cols = std::min( U.Dimension(1), V.Dimension(1) );
        TangentPointMetric_Multiply(alpha, U.data(), beta, V.data(), cols);
    }

    virtual void TangentPointMetric_ApplyKernel(
        const bool farQ,
        const KernelType kernel
    ) const override
    {
        RequireTangentPointMetric();
        
        if( farQ )
        {
            tpm->ApplyFarFieldKernel ( static_cast<Real>(1), kernel );
        }
        else
        {
            tpm->ApplyNearFieldKernel( static_cast<Real>(1), kernel );
        }
    }
