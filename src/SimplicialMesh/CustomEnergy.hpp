#pragma once

//##############################################################################################
//      Custom energy
//##############################################################################################

protected:
    
    mutable ExtReal custom_e_weight = 1;
    
    mutable std::unique_ptr<EnergyBase<Real,Int,SReal,ExtReal>> custom_e;
    

public:
    
    void LoadCustomEnergy( std::unique_ptr<EnergyBase<Real,Int,SReal,ExtReal>> e ) const override
    {
        custom_e = std::move(e);
        if( custom_e != nullptr )
        {
            custom_e->SetWeight(custom_e_weight);
        }
    }
    
    virtual ExtReal GetCustomEnergyWeight() const override
    {
       return custom_e_weight;
    }

    virtual void SetCustomEnergyWeight( const ExtReal weight ) const override
    {
        custom_e_weight = weight;
        
        if( custom_e != nullptr )
        {
            custom_e->SetWeight(custom_e_weight);
        }
    }

    virtual ExtReal CustomEnergy() const override
    {
        if( custom_e != nullptr )
        {
            return custom_e->Value(*this);
        }
        else
        {
            return static_cast<ExtReal>(0);
        }
    };

    virtual ExtReal CustomEnergy_Differential( ExtReal * output, bool addTo = false ) const override
    {
        if( custom_e != nullptr )
        {
            return custom_e->Differential(*this, output, addTo);
        }
        else
        {
            if( !addTo )
            {
                for( Int i = 0; i < DofCount(); ++i )
                {
                    output[i] = static_cast<ExtReal>(0);
                }
            }
            return static_cast<ExtReal>(0);
        }
    }
    
    virtual ExtReal CustomEnergy_Differential( Tensor1<ExtReal,Int> & output, bool addTo = false ) const override
    {
        return CustomEnergy_Differential( output.data(), addTo );
    }
    
    virtual ExtReal CustomEnergy_Differential( Tensor2<ExtReal,Int> & output, bool addTo = false ) const override
    {
        return CustomEnergy_Differential( output.data(), addTo );
    }
    
    virtual void CustomEnergy_Density( ExtReal * output, bool addTo = false ) const override
    {
        custom_e->Density(*this, output, addTo );
    }
    
    virtual void CustomEnergy_Density( Tensor1<ExtReal,Int> & output, bool addTo = false ) const override
    {
        CustomEnergy_Density( output.data(), addTo );
    }
    
    virtual void CustomEnergy_SimplexEnergies( ExtReal * output, bool addTo = false ) const override
    {
        custom_e->SimplexEnergies(*this, output, addTo );
    }
    
    virtual void CustomEnergy_SimplexEnergies( Tensor1<ExtReal,Int> & output, bool addTo = false ) const override
    {
        CustomEnergy_SimplexEnergies( output.data(), addTo );
    }
