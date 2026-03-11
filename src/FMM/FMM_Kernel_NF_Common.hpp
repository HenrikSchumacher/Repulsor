public:

    void PrintReport( float time ) const
    {
        logprint(this->ClassName()+":\t time elapsed: \t" + ToStringFPGeneral(time)
            + "\n"
        );
    }

