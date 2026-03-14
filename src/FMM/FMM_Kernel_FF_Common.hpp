public:

    void PrintReport( float time ) const
    {
        logprint(this->ClassName()+":\t time elapsed: \t" + ToString(time)
            + "\n"
        );
    }
