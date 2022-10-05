#pragma once

//##############################################################################################
//      TangentPointEnergy
//##############################################################################################
           
protected:

    void RequireTangentPointEnergy() const
    {
        if( !tpe_initialized )
        {
            tpe = std::make_unique<
                TP_Energy_FMM_Adaptive<DOM_DIM,AMB_DIM,0,Real,Int,SReal,ExtReal>
            >( tp_alpha, tp_beta, tp_weight, adaptivity_settings );
            
            tpe_initialized = true;
        }
    }

public:
        
    virtual ExtReal TangentPointEnergy() const override
    {
        RequireTangentPointEnergy();
        
        return tpe->Value(*this);
    }

    virtual ExtReal TangentPointEnergy_Differential( ExtReal * output, bool addTo = false ) const override
    {
        RequireTangentPointEnergy();
        
        return tpe->Differential(*this, output, addTo);
    }

    virtual ExtReal TangentPointEnergy_Differential( Tensor1<ExtReal,Int> & output, bool addTo = false ) const override
    {
        return TangentPointEnergy_Differential( output.data(), addTo );
    }
    
    virtual ExtReal TangentPointEnergy_Differential( Tensor2<ExtReal,Int> & output, bool addTo = false ) const override
    {
        return TangentPointEnergy_Differential( output.data(), addTo );
    }
    
    virtual void TangentPointEnergy_Density( ExtReal * output, bool addTo = false ) const override
    {
        RequireTangentPointEnergy();
        
        tpe->Density(*this, output, addTo);
    }

    virtual void TangentPointEnergy_Density( Tensor1<ExtReal,Int> & output, bool addTo = false ) const override
    {
        return TangentPointEnergy_Density( output.data(), addTo );
    }
    
    virtual void TangentPointEnergy_SimplexEnergies( ExtReal * output, bool addTo = false ) const override
    {
        RequireTangentPointEnergy();
        
        tpe->SimplexEnergies(*this, output, addTo);
    }

    virtual void TangentPointEnergy_SimplexEnergies( Tensor1<ExtReal,Int> & output, bool addTo = false ) const override
    {
        return TangentPointEnergy_SimplexEnergies( output.data(), addTo );
    }
    
    
