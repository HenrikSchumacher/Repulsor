protected:
    
    void ComputeNearFarData(
        mref<Tensor2<Real,Int>> P_near,
        mref<Tensor2<Real,Int>> P_far
    ) const
    {
        TOOLS_PTIMER(timer,ClassName()+"::ComputeNearFarData");
        
        ParallelDo(
            [&]( const Int thread )
            {
                SimplexDataKernel_T ker ( V_coords, simplices, V_charges );
                
                const Int i_begin = JobPointer<Int>(SimplexCount(),ThreadCount(),thread  );
                const Int i_end   = JobPointer<Int>(SimplexCount(),ThreadCount(),thread+1);
                
                for( Int i = i_begin; i < i_end; ++i )
                {
                    ker.ReadSimplex(i);
                    ker.WriteNear( P_near.data(i) );
                    ker.WriteFar ( P_far.data(i)  );
                }
            },
            ThreadCount()
        );
};
