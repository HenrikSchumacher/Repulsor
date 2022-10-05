#pragma once

//##############################################################################################
//      TangentPointObstacleEnergy
//##############################################################################################

protected:


    void RequireTangentPointObstacleEnergy() const
    {
        if( !tpo_initialized )
        {
            if( obstacle_initialized )
            {
                Energy_T * r = nullptr;
                
                switch( GetObstacle().DomDim() )
                {
                    case 1:
                    {
                        r = new TP_ObstacleEnergy_FMM_Adaptive<DOM_DIM,1,AMB_DIM,0,Real,Int,SReal,ExtReal> ( tp_alpha, tp_beta, tp_weight, adaptivity_settings );
                        break;
                    }
                    case 2:
                    {
                        r = new TP_ObstacleEnergy_FMM_Adaptive<DOM_DIM,2,AMB_DIM,0,Real,Int,SReal,ExtReal> ( tp_alpha, tp_beta, tp_weight, adaptivity_settings );
                        break;
                    }
                    default:
                    {
                        wprint(className()+"RequireTangentPointObstacleEnergy : domain dimension "+ToString(GetObstacle().DomDim())+" invalid. Only obstacles of domain dimension 1 and 2 are implemented.");
                    }

                }
                
                if( r != nullptr )
                {
                    tpo = std::unique_ptr<Energy_T>(r);
                    
                    tpo_initialized = true;
                }
            }
        }
    }

public:
    virtual ExtReal TangentPointObstacleEnergy() const override
    {
        RequireTangentPointObstacleEnergy();
        
        return tpo->Value(*this);
    };

    virtual ExtReal TangentPointObstacleEnergy_Differential( ExtReal * output, bool addTo = false ) const override
    {
        RequireTangentPointObstacleEnergy();
        
        return tpo->Differential(*this, output, addTo);
    }

    virtual ExtReal TangentPointObstacleEnergy_Differential( Tensor1<ExtReal,Int> & output, bool addTo = false ) const override
    {
        return TangentPointObstacleEnergy_Differential( output.data(), addTo );
    }
    
    virtual ExtReal TangentPointObstacleEnergy_Differential( Tensor2<ExtReal,Int> & output, bool addTo = false ) const override
    {
        return TangentPointObstacleEnergy_Differential( output.data(), addTo );
    }

    virtual void TangentPointObstacleEnergy_Density( ExtReal * output, bool addTo = false ) const override
    {
        RequireTangentPointObstacleEnergy();
        
        tpo->Density(*this, output, addTo);
    }
    
    virtual void TangentPointObstacleEnergy_Density( Tensor1<ExtReal,Int> & output, bool addTo = false ) const override
    {
        return TangentPointObstacleEnergy_Density( output.data(), addTo );
    }
    
    virtual void TangentPointObstacleEnergy_SimplexEnergies( ExtReal * output, bool addTo = false ) const override
    {
        RequireTangentPointObstacleEnergy();
        
        tpo->SimplexEnergies(*this, output, addTo);
    }

    virtual void TangentPointObstacleEnergy_SimplexEnergies( Tensor1<ExtReal,Int> & output, bool addTo = false ) const override
    {
        return TangentPointObstacleEnergy_SimplexEnergies( output.data(), addTo );
    }
    
