protected:
    
    template<bool addtoQ = false>
    void DNearToHulls(
        Tensor3<Real,Int> & restrict buffer,
    ) const
    {
        ptic(ClassName()+"::DNearToHulls");
        
        if constexpr (!addtoQ)
        {
            buffer.SetZero();
        }
        
        ParallelDo(
            [&]( const Int thread )
            {
                
            },
            ThreadCount()
        );
    
        ptoc(ClassName()+"::DNearToHulls");
};

