#pragma once

//##############################################################################################
//      Trivial energy
//##############################################################################################
    
protected:
    
    mutable ExtReal trivial_e_weight = 1;
    
    mutable bool trivial_e_initialized = false;
    mutable std::unique_ptr<TrivialEnergy_FMM_Adaptive<DOM_DIM,AMB_DIM,0,Real,Int,SReal,ExtReal>> trivial_e;
    
    mutable bool trivial_oe_initialized = false;
    mutable std::unique_ptr<Energy_T> trivial_oe;
    
    void RequireTrivialEnergy() const
    {
        if( !trivial_e_initialized )
        {
            trivial_e = std::make_unique< TrivialEnergy_FMM_Adaptive<DOM_DIM,AMB_DIM,0,Real,Int,SReal,ExtReal>
            >(
                trivial_e_weight,
                adaptivity_settings
            );
            
            trivial_e_initialized = true;
        }
    }
    
//        void RequireTrivialObstacleEnergy() const
//        {
//            if( !trivial_oe_initialized )
//            {
//                trivial_oe = std::make_unique< TrivialObstacleEnergy_FMM_Adaptive<DOM_DIM,DOM_DIM,AMB_DIM,0,Real,Int,SReal,ExtReal>
//                >(
//                    trivial_e_weight,
//                    adaptivity_settings
//                );
//
//                trivial_oe_initialized = true;
//            }
//        }
    
    void RequireTrivialObstacleEnergy() const
    {
        if( !trivial_oe_initialized )
        {
            if( obstacle_initialized )
            {
                Energy_T * r = nullptr;
                
                switch( GetObstacle().DomDim() )
                {
                    case 0:
                    {
                        r = new TrivialObstacleEnergy_FMM_Adaptive<DOM_DIM,0,AMB_DIM,0,Real,Int,SReal,ExtReal>
                            (
                                trivial_e_weight,
                                adaptivity_settings
                            );
                        break;
                    }
                    case 1:
                    {
                        r = new TrivialObstacleEnergy_FMM_Adaptive<DOM_DIM,1,AMB_DIM,0,Real,Int,SReal,ExtReal>
                            (
                                trivial_e_weight,
                                adaptivity_settings
                            );
                        break;
                    }
                    case 2:
                    {
                        r = new TrivialObstacleEnergy_FMM_Adaptive<DOM_DIM,2,AMB_DIM,0,Real,Int,SReal,ExtReal>
                            (
                                trivial_e_weight,
                                adaptivity_settings
                            );
                        break;
                    }
                    default:
                    {
                        wprint(className()+"RequireTrivialObstacleEnergy : domain dimension "+ToString(GetObstacle().DomDim())+" invalid. Only obstacles of domain dimension 1 and 2 are implemented.");
                    }
    
                }
                
                if( r != nullptr )
                {
                    trivial_oe = std::unique_ptr<Energy_T>(r);
                    
                    trivial_oe_initialized = true;
                }
            }
        }
    }

public:
    
    virtual ExtReal GetTrivialEnergyWeight() const override
    {
       return trivial_e_weight;
    }

    virtual void SetTrivialEnergyWeight( const ExtReal weight ) const override
    {
        trivial_e_weight = weight;
        if( trivial_e_initialized )
        {
            trivial_e->SetWeight(tp_weight);
        }
    }

    virtual ExtReal TrivialEnergy() const override
    {
        RequireTrivialEnergy();
        
        return trivial_e->Value(*this);
    };

    virtual ExtReal TrivialEnergy_Differential( ExtReal * output, bool addTo = false ) const override
    {
        RequireTrivialEnergy();
        
        return trivial_e->Differential(*this, output, addTo);
    };
    
    virtual ExtReal TrivialEnergy_Differential( Tensor1<ExtReal,Int> & output, bool addTo = false ) const override
    {
        return TrivialEnergy_Differential(output.data(), addTo );
    }
    
    virtual  ExtReal TrivialEnergy_Differential( Tensor2<ExtReal,Int> & output, bool addTo = false ) const override
    {
        return TrivialEnergy_Differential(output.data(), addTo );
    }
    
    virtual ExtReal TrivialObstacleEnergy() const override
    {
        RequireTrivialObstacleEnergy();
        
        return trivial_oe->Value(*this);
    }

    virtual ExtReal TrivialObstacleEnergy_Differential( ExtReal * output, bool addTo = false ) const override
    {
        RequireTrivialObstacleEnergy();
        
        return trivial_oe->Differential(*this, output, addTo);
    }
    
    virtual ExtReal TrivialObstacleEnergy_Differential( Tensor1<ExtReal,Int> & output, bool addTo = false ) const override
    {
        return TrivialObstacleEnergy_Differential( output.data(), addTo);
    }
    
    virtual ExtReal TrivialObstacleEnergy_Differential( Tensor2<ExtReal,Int> & output, bool addTo = false ) const override
    {
        return TrivialObstacleEnergy_Differential( output.data(), addTo);
    }
