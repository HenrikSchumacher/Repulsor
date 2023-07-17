protected:
    
    void ComputeNearFarData(
        Tensor2<Real,Int> & restrict P_near,
        Tensor2<Real,Int> & restrict P_far
    ) const
    {
        ptic(ClassName()+"::ComputeNearFarData");
        
        ParallelDo(
            [&]( const Int thread )
            {
                PrimitiveDataKernel_T ker ( V_coords, simplices, V_charges );
                
                const Int i_begin = JobPointer<Int>(SimplexCount(),ThreadCount(),thread  );
                const Int i_end   = JobPointer<Int>(SimplexCount(),ThreadCount(),thread+1);
                
                for( Int i = i_begin; i < i_end; ++i )
                {
                    ker.ReadPrimitive(i);
                    
                    ker.WriteNear( P_near.data(i) );
                    ker.WriteFar ( P_far.data(i)  );
                }
            },
            ThreadCount()
        );
    
        ptoc(ClassName()+"::ComputeNearFarData");
};
