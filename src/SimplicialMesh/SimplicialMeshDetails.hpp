// Machine generated code. Don't edit this file!

#pragma once

namespace Repulsor
{
	template<int DOM_DIM, int AMB_DIM, typename Real, typename Int>
    struct SimplicialMeshDetails
    {
		inline std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<")+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+">";
        }

		
	void ComputeNearFarData(
		const Tensor2<Real,Int> & V_coords,
		const Tensor2<Int ,Int> & simplices,
			  Tensor2<Real,Int> & P_near,
			  Tensor2<Real,Int> & P_far
	) const
    {
        ptic(ClassName()+"::ComputeNearFarData");
        eprint(ClassName()+"::ComputeNearFarData not implemented. Doing nothing.");
        ptoc(ClassName()+"::ComputeNearFarData");
    }

	void ComputeNearFarDataOps( 
		const Tensor2<Real,Int> & V_coords,
        const Tensor2<Int ,Int> & simplices,
		      Tensor2<Real,Int> & P_coords,
		      Tensor3<Real,Int> & P_hull_coords,
		      Tensor2<Real,Int> & P_near,
		      Tensor2<Real,Int> & P_far,
		Sparse::MatrixCSR<Real,Int,Int> & DiffOp,
		Sparse::MatrixCSR<Real,Int,Int> & AvOp 
	) const
    {
        ptic(ClassName()+"::ComputeNearFarDataOps");
        eprint(ClassName()+"::ComputeNearFarDataOps not implemented. Doing nothing.");
        ptoc(ClassName()+"::ComputeNearFarDataOps");
    }

	void DNearToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_near, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DNearToHulls");
        eprint(ClassName()+"::DNearToHulls not implemented. Returning 0.");
		
		if(!addTo)
		{
			buffer.Fill(static_cast<Real>(0));
		}
        ptoc(ClassName()+"::DNearToHulls");
    }

	void DFarToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_far,
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DFarToHulls");
                eprint(ClassName()+"::DFarToHulls not implemented. Returning 0.");
		
		if(!addTo)
		{
			buffer.Fill(static_cast<Real>(0));
		}
        ptoc(ClassName()+"::DFarToHulls");
    }

	}; // SimplicialMeshDetails<DOM_DIM,AMB,Real,Int>

//----------------------------------------------------------------------------------------------

    
    template<typename Real, typename Int>
    struct SimplicialMeshDetails<0,1,Real,Int>
    {
	private:

		const Int thread_count = 1;

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(std::max(static_cast<Int>(1),thread_count_))
		{}
	
		inline Int ThreadCount() const
		{
			return thread_count;
		}
	
		inline std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<0,1,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void ComputeNearFarData( 
		const Tensor2<Real,Int> & V_coords,
		const Tensor2<Int ,Int> & simplices,
			  Tensor2<Real,Int> & P_near,
			  Tensor2<Real,Int> & P_far
	) const
    {
        ptic(ClassName()+"::ComputeNearFarData");
        
        //Int size       = 1;
        //Int amb_dim    = 1;
        //Int dom_dim    = 0;

        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			ptr<Real> V_coords__      = V_coords.data();	
			ptr<Int>  simplices__     = simplices.data();

			Real hull    [1][1];

			Int simplex  [1];
			
			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {
				ptr<Real> near = P_near.data(i);                    
				ptr<Real> far  = P_far.data(i);   
            
				simplex[0] = simplices__[1*i +0];

				near[1] = hull[0][0] = V_coords__[1*simplex[0]+0];

				far[1] = static_cast<Real>(1.) * ( hull[0][0] );

				near[0] = far[0] = static_cast<Real>(1);
				near[2] = far[2] = static_cast<Real>(1);

            } // for( Int i = i_begin; i < i_end; ++i )

		} // #pragma omp parallel num_threads( ThreadCount() )

        ptoc(ClassName()+"::ComputeNearFarData");
    }

    void ComputeNearFarDataOps( 
		const Tensor2<Real,Int> & V_coords,
        const Tensor2<Int ,Int> & simplices,
		      Tensor2<Real,Int> & P_coords,
		      Tensor3<Real,Int> & P_hull_coords,
		      Tensor2<Real,Int> & P_near,
		      Tensor2<Real,Int> & P_far,
		Sparse::MatrixCSR<Real,Int,Int> & DiffOp,
		Sparse::MatrixCSR<Real,Int,Int> & AvOp 
	) const
    {
        ptic(ClassName()+"::ComputeNearFarDataOps");
        
        //Int size       = 1;
        //Int amb_dim    = 1;
        //Int dom_dim    = 0;
        
        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			ptr<Int> AvOp_outer = AvOp.Outer().data();
			ptr<Int> AvOp_inner = AvOp.Inner().data();
			AvOp.Value().Fill(static_cast<Real>(1));

			mut<Int> DiffOp_outer = DiffOp.Outer().data();
			mut<Int> DiffOp_inner = DiffOp.Inner().data();
			DiffOp.Value().SetZero();

			ptr<Real> V_coords__      = V_coords.data();
			
			ptr<Int>  simplices__     = simplices.data();
			mut<Real> P_hull_coords__ = P_hull_coords.data();
		    mut<Real> P_coords__      = P_coords.data();

			Int simplex        [1];
			Int sorted_simplex [1];

			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {

				mut<Real> near = P_near.data(i);                    
				mut<Real> far  = P_far.data(i);

				simplex[0] = sorted_simplex[0] = simplices__[1*i +0];
                  
                // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
                std::sort( sorted_simplex, sorted_simplex + 1 );

				AvOp_outer[i+1] = (i+1) * 1;  
                      
				AvOp_inner[1*i+0] = sorted_simplex[0];

				DiffOp_outer[1*i+0] = (1 * i + 0) * 1;

				DiffOp_inner[(i*1+0)*1+0] = sorted_simplex[0];

				near[1] = P_hull_coords__[1*i+0] = V_coords__[1*simplex[0]+0];

				far[1] = P_coords__[1*i+0] = 1. * ( P_hull_coords__[1*i+0] );

                near[0] = far[0] = static_cast<Real>(1);
				near[2] = far[2] = static_cast<Real>(1);

            }
        }
        ptoc(ClassName()+"::ComputeNearFarDataOps");
    }

	void DNearToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_near, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DNearToHulls");

        if( P_D_near.Dimension(1) != 3 )
        {
            eprint("in DNearToHulls: P_D_near.Dimension(1) != 3. Aborting");
        }

		//ptr<Real> V_coords__  = V_coords.data();
		//ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_near__  = P_D_near.data();
		mut<Real> buffer__    = buffer.data();
        
        if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[1*i+0] += P_D_near__[3*i+1];
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[1*i+0] = P_D_near__[3*i+1];
			}
		}

        ptoc(ClassName()+"::DNearToHulls");
        
    }

	void DFarToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_far, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DFarToHulls");

        if( P_D_far.Dimension(1) != 3 )
        {
            eprint("in DFarToHulls: P_D_far.Dimension(1) != 3. Aborting");
        }

		//ptr<Real> V_coords__  = V_coords.data();
		//ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_far__   = P_D_far.data();
		mut<Int>  buffer__    = buffer.data();
        
        if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[1*i+0] += P_D_far__[3*i+1];
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[1*i+0] = P_D_far__[3*i+1];
			}
		}

        ptoc(ClassName()+"::DFarToHulls");
        
    }

	}; // SimplicialMeshDetails<0,1,Real,Int>

//----------------------------------------------------------------------------------------------

    template<typename Real, typename Int>
    struct SimplicialMeshDetails<0,2,Real,Int>
    {
	private:

		const Int thread_count = 1;

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(std::max(static_cast<Int>(1),thread_count_))
		{}
	
		inline Int ThreadCount() const
		{
			return thread_count;
		}
	
		inline std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<0,2,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void ComputeNearFarData( 
		const Tensor2<Real,Int> & V_coords,
		const Tensor2<Int ,Int> & simplices,
			  Tensor2<Real,Int> & P_near,
			  Tensor2<Real,Int> & P_far
	) const
    {
        ptic(ClassName()+"::ComputeNearFarData");
        
        //Int size       = 1;
        //Int amb_dim    = 2;
        //Int dom_dim    = 0;

        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			ptr<Real> V_coords__      = V_coords.data();	
			ptr<Int>  simplices__     = simplices.data();

			Real hull    [1][2];

			Int simplex  [1];
			
			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {
				ptr<Real> near = P_near.data(i);                    
				ptr<Real> far  = P_far.data(i);   
            
				simplex[0] = simplices__[1*i +0];

				near[1] = hull[0][0] = V_coords__[2*simplex[0]+0];
				near[2] = hull[0][1] = V_coords__[2*simplex[0]+1];

				far[1] = static_cast<Real>(1.) * ( hull[0][0] );
				far[2] = static_cast<Real>(1.) * ( hull[0][1] );

				near[0] = far[0] = static_cast<Real>(1);
				near[3] = far[3] = static_cast<Real>(1);
				near[4] = far[4] = static_cast<Real>(0);
				near[5] = far[5] = static_cast<Real>(1);

            } // for( Int i = i_begin; i < i_end; ++i )

		} // #pragma omp parallel num_threads( ThreadCount() )

        ptoc(ClassName()+"::ComputeNearFarData");
    }

    void ComputeNearFarDataOps( 
		const Tensor2<Real,Int> & V_coords,
        const Tensor2<Int ,Int> & simplices,
		      Tensor2<Real,Int> & P_coords,
		      Tensor3<Real,Int> & P_hull_coords,
		      Tensor2<Real,Int> & P_near,
		      Tensor2<Real,Int> & P_far,
		Sparse::MatrixCSR<Real,Int,Int> & DiffOp,
		Sparse::MatrixCSR<Real,Int,Int> & AvOp 
	) const
    {
        ptic(ClassName()+"::ComputeNearFarDataOps");
        
        //Int size       = 1;
        //Int amb_dim    = 2;
        //Int dom_dim    = 0;
        
        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			ptr<Int> AvOp_outer = AvOp.Outer().data();
			ptr<Int> AvOp_inner = AvOp.Inner().data();
			AvOp.Value().Fill(static_cast<Real>(1));

			mut<Int> DiffOp_outer = DiffOp.Outer().data();
			mut<Int> DiffOp_inner = DiffOp.Inner().data();
			DiffOp.Value().SetZero();

			ptr<Real> V_coords__      = V_coords.data();
			
			ptr<Int>  simplices__     = simplices.data();
			mut<Real> P_hull_coords__ = P_hull_coords.data();
		    mut<Real> P_coords__      = P_coords.data();

			Int simplex        [1];
			Int sorted_simplex [1];

			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {

				mut<Real> near = P_near.data(i);                    
				mut<Real> far  = P_far.data(i);

				simplex[0] = sorted_simplex[0] = simplices__[1*i +0];
                  
                // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
                std::sort( sorted_simplex, sorted_simplex + 1 );

				AvOp_outer[i+1] = (i+1) * 1;  
                      
				AvOp_inner[1*i+0] = sorted_simplex[0];

				DiffOp_outer[2*i+0] = (2 * i + 0) * 1;
				DiffOp_outer[2*i+1] = (2 * i + 1) * 1;

				DiffOp_inner[(i*2+0)*1+0] = sorted_simplex[0];
				DiffOp_inner[(i*2+1)*1+0] = sorted_simplex[0];

				near[1] = P_hull_coords__[2*i+0] = V_coords__[2*simplex[0]+0];
				near[2] = P_hull_coords__[2*i+1] = V_coords__[2*simplex[0]+1];

				far[1] = P_coords__[2*i+0] = 1. * ( P_hull_coords__[2*i+0] );
				far[2] = P_coords__[2*i+1] = 1. * ( P_hull_coords__[2*i+1] );

                near[0] = far[0] = static_cast<Real>(1);
				near[3] = far[3] = static_cast<Real>(1);
				near[4] = far[4] = static_cast<Real>(0);
				near[5] = far[5] = static_cast<Real>(1);

            }
        }
        ptoc(ClassName()+"::ComputeNearFarDataOps");
    }

	void DNearToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_near, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DNearToHulls");

        if( P_D_near.Dimension(1) != 6 )
        {
            eprint("in DNearToHulls: P_D_near.Dimension(1) != 6. Aborting");
        }

		//ptr<Real> V_coords__  = V_coords.data();
		//ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_near__  = P_D_near.data();
		mut<Real> buffer__    = buffer.data();
        
        if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[2*i+0] += P_D_near__[6*i+1];
				buffer__[2*i+1] += P_D_near__[6*i+2];
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[2*i+0] = P_D_near__[6*i+1];
				buffer__[2*i+1] = P_D_near__[6*i+2];
			}
		}

        ptoc(ClassName()+"::DNearToHulls");
        
    }

	void DFarToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_far, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DFarToHulls");

        if( P_D_far.Dimension(1) != 6 )
        {
            eprint("in DFarToHulls: P_D_far.Dimension(1) != 6. Aborting");
        }

		//ptr<Real> V_coords__  = V_coords.data();
		//ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_far__   = P_D_far.data();
		mut<Int>  buffer__    = buffer.data();
        
        if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[2*i+0] += P_D_far__[6*i+1];
				buffer__[2*i+1] += P_D_far__[6*i+2];
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[2*i+0] = P_D_far__[6*i+1];
				buffer__[2*i+1] = P_D_far__[6*i+2];
			}
		}

        ptoc(ClassName()+"::DFarToHulls");
        
    }

	}; // SimplicialMeshDetails<0,2,Real,Int>

//----------------------------------------------------------------------------------------------

    template<typename Real, typename Int>
    struct SimplicialMeshDetails<0,3,Real,Int>
    {
	private:

		const Int thread_count = 1;

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(std::max(static_cast<Int>(1),thread_count_))
		{}
	
		inline Int ThreadCount() const
		{
			return thread_count;
		}
	
		inline std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<0,3,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void ComputeNearFarData( 
		const Tensor2<Real,Int> & V_coords,
		const Tensor2<Int ,Int> & simplices,
			  Tensor2<Real,Int> & P_near,
			  Tensor2<Real,Int> & P_far
	) const
    {
        ptic(ClassName()+"::ComputeNearFarData");
        
        //Int size       = 1;
        //Int amb_dim    = 3;
        //Int dom_dim    = 0;

        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			ptr<Real> V_coords__      = V_coords.data();	
			ptr<Int>  simplices__     = simplices.data();

			Real hull    [1][3];

			Int simplex  [1];
			
			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {
				ptr<Real> near = P_near.data(i);                    
				ptr<Real> far  = P_far.data(i);   
            
				simplex[0] = simplices__[1*i +0];

				near[1] = hull[0][0] = V_coords__[3*simplex[0]+0];
				near[2] = hull[0][1] = V_coords__[3*simplex[0]+1];
				near[3] = hull[0][2] = V_coords__[3*simplex[0]+2];

				far[1] = static_cast<Real>(1.) * ( hull[0][0] );
				far[2] = static_cast<Real>(1.) * ( hull[0][1] );
				far[3] = static_cast<Real>(1.) * ( hull[0][2] );

				near[ 0] = far[ 0] = static_cast<Real>(1);
				near[ 4] = far[ 4] = static_cast<Real>(1);
				near[ 5] = far[ 5] = static_cast<Real>(0);
				near[ 6] = far[ 6] = static_cast<Real>(0);
				near[ 7] = far[ 7] = static_cast<Real>(1);
				near[ 8] = far[ 8] = static_cast<Real>(0);
				near[ 9] = far[ 9] = static_cast<Real>(1);

            } // for( Int i = i_begin; i < i_end; ++i )

		} // #pragma omp parallel num_threads( ThreadCount() )

        ptoc(ClassName()+"::ComputeNearFarData");
    }

    void ComputeNearFarDataOps( 
		const Tensor2<Real,Int> & V_coords,
        const Tensor2<Int ,Int> & simplices,
		      Tensor2<Real,Int> & P_coords,
		      Tensor3<Real,Int> & P_hull_coords,
		      Tensor2<Real,Int> & P_near,
		      Tensor2<Real,Int> & P_far,
		Sparse::MatrixCSR<Real,Int,Int> & DiffOp,
		Sparse::MatrixCSR<Real,Int,Int> & AvOp 
	) const
    {
        ptic(ClassName()+"::ComputeNearFarDataOps");
        
        //Int size       = 1;
        //Int amb_dim    = 3;
        //Int dom_dim    = 0;
        
        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			ptr<Int> AvOp_outer = AvOp.Outer().data();
			ptr<Int> AvOp_inner = AvOp.Inner().data();
			AvOp.Value().Fill(static_cast<Real>(1));

			mut<Int> DiffOp_outer = DiffOp.Outer().data();
			mut<Int> DiffOp_inner = DiffOp.Inner().data();
			DiffOp.Value().SetZero();

			ptr<Real> V_coords__      = V_coords.data();
			
			ptr<Int>  simplices__     = simplices.data();
			mut<Real> P_hull_coords__ = P_hull_coords.data();
		    mut<Real> P_coords__      = P_coords.data();

			Int simplex        [1];
			Int sorted_simplex [1];

			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {

				mut<Real> near = P_near.data(i);                    
				mut<Real> far  = P_far.data(i);

				simplex[0] = sorted_simplex[0] = simplices__[1*i +0];
                  
                // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
                std::sort( sorted_simplex, sorted_simplex + 1 );

				AvOp_outer[i+1] = (i+1) * 1;  
                      
				AvOp_inner[1*i+0] = sorted_simplex[0];

				DiffOp_outer[3*i+0] = (3 * i + 0) * 1;
				DiffOp_outer[3*i+1] = (3 * i + 1) * 1;
				DiffOp_outer[3*i+2] = (3 * i + 2) * 1;

				DiffOp_inner[(i*3+0)*1+0] = sorted_simplex[0];
				DiffOp_inner[(i*3+1)*1+0] = sorted_simplex[0];
				DiffOp_inner[(i*3+2)*1+0] = sorted_simplex[0];

				near[1] = P_hull_coords__[3*i+0] = V_coords__[3*simplex[0]+0];
				near[2] = P_hull_coords__[3*i+1] = V_coords__[3*simplex[0]+1];
				near[3] = P_hull_coords__[3*i+2] = V_coords__[3*simplex[0]+2];

				far[1] = P_coords__[3*i+0] = 1. * ( P_hull_coords__[3*i+0] );
				far[2] = P_coords__[3*i+1] = 1. * ( P_hull_coords__[3*i+1] );
				far[3] = P_coords__[3*i+2] = 1. * ( P_hull_coords__[3*i+2] );

                near[ 0] = far[ 0] = static_cast<Real>(1);
				near[ 4] = far[ 4] = static_cast<Real>(1);
				near[ 5] = far[ 5] = static_cast<Real>(0);
				near[ 6] = far[ 6] = static_cast<Real>(0);
				near[ 7] = far[ 7] = static_cast<Real>(1);
				near[ 8] = far[ 8] = static_cast<Real>(0);
				near[ 9] = far[ 9] = static_cast<Real>(1);

            }
        }
        ptoc(ClassName()+"::ComputeNearFarDataOps");
    }

	void DNearToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_near, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DNearToHulls");

        if( P_D_near.Dimension(1) != 10 )
        {
            eprint("in DNearToHulls: P_D_near.Dimension(1) != 10. Aborting");
        }

		//ptr<Real> V_coords__  = V_coords.data();
		//ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_near__  = P_D_near.data();
		mut<Real> buffer__    = buffer.data();
        
        if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[3*i+0] += P_D_near__[10*i+1];
				buffer__[3*i+1] += P_D_near__[10*i+2];
				buffer__[3*i+2] += P_D_near__[10*i+3];
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[3*i+0] = P_D_near__[10*i+1];
				buffer__[3*i+1] = P_D_near__[10*i+2];
				buffer__[3*i+2] = P_D_near__[10*i+3];
			}
		}

        ptoc(ClassName()+"::DNearToHulls");
        
    }

	void DFarToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_far, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DFarToHulls");

        if( P_D_far.Dimension(1) != 10 )
        {
            eprint("in DFarToHulls: P_D_far.Dimension(1) != 10. Aborting");
        }

		//ptr<Real> V_coords__  = V_coords.data();
		//ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_far__   = P_D_far.data();
		mut<Int>  buffer__    = buffer.data();
        
        if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[3*i+0] += P_D_far__[10*i+1];
				buffer__[3*i+1] += P_D_far__[10*i+2];
				buffer__[3*i+2] += P_D_far__[10*i+3];
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[3*i+0] = P_D_far__[10*i+1];
				buffer__[3*i+1] = P_D_far__[10*i+2];
				buffer__[3*i+2] = P_D_far__[10*i+3];
			}
		}

        ptoc(ClassName()+"::DFarToHulls");
        
    }

	}; // SimplicialMeshDetails<0,3,Real,Int>

//----------------------------------------------------------------------------------------------

    template<typename Real, typename Int>
    struct SimplicialMeshDetails<0,4,Real,Int>
    {
	private:

		const Int thread_count = 1;

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(std::max(static_cast<Int>(1),thread_count_))
		{}
	
		inline Int ThreadCount() const
		{
			return thread_count;
		}
	
		inline std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<0,4,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void ComputeNearFarData( 
		const Tensor2<Real,Int> & V_coords,
		const Tensor2<Int ,Int> & simplices,
			  Tensor2<Real,Int> & P_near,
			  Tensor2<Real,Int> & P_far
	) const
    {
        ptic(ClassName()+"::ComputeNearFarData");
        
        //Int size       = 1;
        //Int amb_dim    = 4;
        //Int dom_dim    = 0;

        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			ptr<Real> V_coords__      = V_coords.data();	
			ptr<Int>  simplices__     = simplices.data();

			Real hull    [1][4];

			Int simplex  [1];
			
			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {
				ptr<Real> near = P_near.data(i);                    
				ptr<Real> far  = P_far.data(i);   
            
				simplex[0] = simplices__[1*i +0];

				near[1] = hull[0][0] = V_coords__[4*simplex[0]+0];
				near[2] = hull[0][1] = V_coords__[4*simplex[0]+1];
				near[3] = hull[0][2] = V_coords__[4*simplex[0]+2];
				near[4] = hull[0][3] = V_coords__[4*simplex[0]+3];

				far[1] = static_cast<Real>(1.) * ( hull[0][0] );
				far[2] = static_cast<Real>(1.) * ( hull[0][1] );
				far[3] = static_cast<Real>(1.) * ( hull[0][2] );
				far[4] = static_cast<Real>(1.) * ( hull[0][3] );

				near[ 0] = far[ 0] = static_cast<Real>(1);
				near[ 5] = far[ 5] = static_cast<Real>(1);
				near[ 6] = far[ 6] = static_cast<Real>(0);
				near[ 7] = far[ 7] = static_cast<Real>(0);
				near[ 8] = far[ 8] = static_cast<Real>(0);
				near[ 9] = far[ 9] = static_cast<Real>(1);
				near[10] = far[10] = static_cast<Real>(0);
				near[11] = far[11] = static_cast<Real>(0);
				near[12] = far[12] = static_cast<Real>(1);
				near[13] = far[13] = static_cast<Real>(0);
				near[14] = far[14] = static_cast<Real>(1);

            } // for( Int i = i_begin; i < i_end; ++i )

		} // #pragma omp parallel num_threads( ThreadCount() )

        ptoc(ClassName()+"::ComputeNearFarData");
    }

    void ComputeNearFarDataOps( 
		const Tensor2<Real,Int> & V_coords,
        const Tensor2<Int ,Int> & simplices,
		      Tensor2<Real,Int> & P_coords,
		      Tensor3<Real,Int> & P_hull_coords,
		      Tensor2<Real,Int> & P_near,
		      Tensor2<Real,Int> & P_far,
		Sparse::MatrixCSR<Real,Int,Int> & DiffOp,
		Sparse::MatrixCSR<Real,Int,Int> & AvOp 
	) const
    {
        ptic(ClassName()+"::ComputeNearFarDataOps");
        
        //Int size       = 1;
        //Int amb_dim    = 4;
        //Int dom_dim    = 0;
        
        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			ptr<Int> AvOp_outer = AvOp.Outer().data();
			ptr<Int> AvOp_inner = AvOp.Inner().data();
			AvOp.Value().Fill(static_cast<Real>(1));

			mut<Int> DiffOp_outer = DiffOp.Outer().data();
			mut<Int> DiffOp_inner = DiffOp.Inner().data();
			DiffOp.Value().SetZero();

			ptr<Real> V_coords__      = V_coords.data();
			
			ptr<Int>  simplices__     = simplices.data();
			mut<Real> P_hull_coords__ = P_hull_coords.data();
		    mut<Real> P_coords__      = P_coords.data();

			Int simplex        [1];
			Int sorted_simplex [1];

			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {

				mut<Real> near = P_near.data(i);                    
				mut<Real> far  = P_far.data(i);

				simplex[0] = sorted_simplex[0] = simplices__[1*i +0];
                  
                // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
                std::sort( sorted_simplex, sorted_simplex + 1 );

				AvOp_outer[i+1] = (i+1) * 1;  
                      
				AvOp_inner[1*i+0] = sorted_simplex[0];

				DiffOp_outer[4*i+0] = (4 * i + 0) * 1;
				DiffOp_outer[4*i+1] = (4 * i + 1) * 1;
				DiffOp_outer[4*i+2] = (4 * i + 2) * 1;
				DiffOp_outer[4*i+3] = (4 * i + 3) * 1;

				DiffOp_inner[(i*4+0)*1+0] = sorted_simplex[0];
				DiffOp_inner[(i*4+1)*1+0] = sorted_simplex[0];
				DiffOp_inner[(i*4+2)*1+0] = sorted_simplex[0];
				DiffOp_inner[(i*4+3)*1+0] = sorted_simplex[0];

				near[1] = P_hull_coords__[4*i+0] = V_coords__[4*simplex[0]+0];
				near[2] = P_hull_coords__[4*i+1] = V_coords__[4*simplex[0]+1];
				near[3] = P_hull_coords__[4*i+2] = V_coords__[4*simplex[0]+2];
				near[4] = P_hull_coords__[4*i+3] = V_coords__[4*simplex[0]+3];

				far[1] = P_coords__[4*i+0] = 1. * ( P_hull_coords__[4*i+0] );
				far[2] = P_coords__[4*i+1] = 1. * ( P_hull_coords__[4*i+1] );
				far[3] = P_coords__[4*i+2] = 1. * ( P_hull_coords__[4*i+2] );
				far[4] = P_coords__[4*i+3] = 1. * ( P_hull_coords__[4*i+3] );

                near[ 0] = far[ 0] = static_cast<Real>(1);
				near[ 5] = far[ 5] = static_cast<Real>(1);
				near[ 6] = far[ 6] = static_cast<Real>(0);
				near[ 7] = far[ 7] = static_cast<Real>(0);
				near[ 8] = far[ 8] = static_cast<Real>(0);
				near[ 9] = far[ 9] = static_cast<Real>(1);
				near[10] = far[10] = static_cast<Real>(0);
				near[11] = far[11] = static_cast<Real>(0);
				near[12] = far[12] = static_cast<Real>(1);
				near[13] = far[13] = static_cast<Real>(0);
				near[14] = far[14] = static_cast<Real>(1);

            }
        }
        ptoc(ClassName()+"::ComputeNearFarDataOps");
    }

	void DNearToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_near, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DNearToHulls");

        if( P_D_near.Dimension(1) != 15 )
        {
            eprint("in DNearToHulls: P_D_near.Dimension(1) != 15. Aborting");
        }

		//ptr<Real> V_coords__  = V_coords.data();
		//ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_near__  = P_D_near.data();
		mut<Real> buffer__    = buffer.data();
        
        if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[4*i+0] += P_D_near__[15*i+1];
				buffer__[4*i+1] += P_D_near__[15*i+2];
				buffer__[4*i+2] += P_D_near__[15*i+3];
				buffer__[4*i+3] += P_D_near__[15*i+4];
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[4*i+0] = P_D_near__[15*i+1];
				buffer__[4*i+1] = P_D_near__[15*i+2];
				buffer__[4*i+2] = P_D_near__[15*i+3];
				buffer__[4*i+3] = P_D_near__[15*i+4];
			}
		}

        ptoc(ClassName()+"::DNearToHulls");
        
    }

	void DFarToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_far, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DFarToHulls");

        if( P_D_far.Dimension(1) != 15 )
        {
            eprint("in DFarToHulls: P_D_far.Dimension(1) != 15. Aborting");
        }

		//ptr<Real> V_coords__  = V_coords.data();
		//ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_far__   = P_D_far.data();
		mut<Int>  buffer__    = buffer.data();
        
        if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[4*i+0] += P_D_far__[15*i+1];
				buffer__[4*i+1] += P_D_far__[15*i+2];
				buffer__[4*i+2] += P_D_far__[15*i+3];
				buffer__[4*i+3] += P_D_far__[15*i+4];
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				buffer__[4*i+0] = P_D_far__[15*i+1];
				buffer__[4*i+1] = P_D_far__[15*i+2];
				buffer__[4*i+2] = P_D_far__[15*i+3];
				buffer__[4*i+3] = P_D_far__[15*i+4];
			}
		}

        ptoc(ClassName()+"::DFarToHulls");
        
    }

	}; // SimplicialMeshDetails<0,4,Real,Int>

//----------------------------------------------------------------------------------------------

    template<typename Real, typename Int>
    struct SimplicialMeshDetails<1,2,Real,Int>
    {
	private:

		const Int thread_count = 1;

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(std::max(static_cast<Int>(1),thread_count_))
		{}
	
		inline Int ThreadCount() const
		{
			return thread_count;
		}
	
		inline std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<1,2,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void ComputeNearFarData( 
		const Tensor2<Real,Int> & V_coords,
		const Tensor2<Int ,Int> & simplices,
			  Tensor2<Real,Int> & P_near,
			  Tensor2<Real,Int> & P_far
	) const
    {
        ptic(ClassName()+"::ComputeNearFarData");
        
        //Int size       = 2;
        //Int amb_dim    = 2;
        //Int dom_dim    = 1;

        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			ptr<Real> V_coords__      = V_coords.data();	
			const Int  * restrict const simplices__     = simplices.data();

			Real hull    [2][2];
			Real df      [2][1];
			Real dfdagger[1][2];
			Real g       [1][1];
			Real ginv    [1][1];

			Int simplex  [2];
			
			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {
				Real * restrict const near = P_near.data(i);                    
				Real * restrict const far  = P_far.data(i);   
            
				simplex[0] = simplices__[2*i +0];
				simplex[1] = simplices__[2*i +1];

				near[1] = hull[0][0] = V_coords__[2*simplex[0]+0];
				near[2] = hull[0][1] = V_coords__[2*simplex[0]+1];
				near[3] = hull[1][0] = V_coords__[2*simplex[1]+0];
				near[4] = hull[1][1] = V_coords__[2*simplex[1]+1];

				far[1] = static_cast<Real>(0.5) * ( hull[0][0] + hull[1][0] );
				far[2] = static_cast<Real>(0.5) * ( hull[0][1] + hull[1][1] );

				df[0][0] = hull[1][0] - hull[0][0];
				df[1][0] = hull[1][1] - hull[0][1];

				g[0][0] = df[0][0] * df[0][0] + df[1][0] * df[1][0];


				near[0] = far[0] = sqrt( fabs(g[0][0]) );

                ginv[0][0] = static_cast<Real>(1)/g[0][0];
                
                //  dfdagger = g^{-1} * df^T (1 x 2 matrix)
				dfdagger[0][0] = ginv[0][0] * df[0][0];
				dfdagger[0][1] = ginv[0][0] * df[1][0];
            
				near[5] = far[3]  = static_cast<Real>(1) - df[0][0] * dfdagger[0][0];
				near[6] = far[4]  =    - df[0][0] * dfdagger[0][1];
				near[7] = far[5]  = static_cast<Real>(1) - df[1][0] * dfdagger[0][1];

            } // for( Int i = i_begin; i < i_end; ++i )

		} // #pragma omp parallel num_threads( ThreadCount() )

        ptoc(ClassName()+"::ComputeNearFarData");
    }

    void ComputeNearFarDataOps( 
		const Tensor2<Real,Int> & V_coords,
        const Tensor2<Int ,Int> & simplices,
		      Tensor2<Real,Int> & P_coords,
		      Tensor3<Real,Int> & P_hull_coords,
		      Tensor2<Real,Int> & P_near,
		      Tensor2<Real,Int> & P_far,
		Sparse::MatrixCSR<Real,Int,Int> & DiffOp,
		Sparse::MatrixCSR<Real,Int,Int> & AvOp 
	) const
    {
        ptic(ClassName()+"::ComputeNearFarDataOps");
        
        //Int size       = 2;
        //Int amb_dim    = 2;
        //Int dom_dim    = 1;
        
        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			mut<Int>  AvOp_outer = AvOp.Outer().data();
			mut<Int>  AvOp_inner = AvOp.Inner().data();
			mut<Real> AvOp_value = AvOp.Values().data();

			mut<Int>  DiffOp_outer = DiffOp.Outer().data();
			mut<Int>  DiffOp_inner = DiffOp.Inner().data();
			mut<Real> DiffOp_value = DiffOp.Value().data();

			ptr<Real> V_coords__      = V_coords.data();
			
			ptr<Int>  simplices__     = simplices.data();
		    mut<Real> P_hull_coords__ = P_hull_coords.data();
			mut<Real> P_coords__      = P_coords.data();

			Real df       [2][1];
			Real dfdagger [1][2];
			Real g        [1][1];
			Real ginv     [1][1];

			Int simplex        [2];
			Int sorted_simplex [2];

			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {

				mut<Real> near = P_near.data(i);                    
				mut<Real> far  = P_far.data(i);

				simplex[0] = sorted_simplex[0] = simplices__[2*i +0];
				simplex[1] = sorted_simplex[1] = simplices__[2*i +1];
                  
                // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
                std::sort( sorted_simplex, sorted_simplex + 2 );

				AvOp_outer[i+1] = (i+1) * 2;  
                      
				AvOp_inner[2*i+0] = sorted_simplex[0];
				AvOp_inner[2*i+1] = sorted_simplex[1];

				AvOp_value[2*i+0] = 0.5;
				AvOp_value[2*i+1] = 0.5;

				DiffOp_outer[2*i+0] = (2 * i + 0) * 2;
				DiffOp_outer[2*i+1] = (2 * i + 1) * 2;

				DiffOp_inner[(i*2+0)*2+0] = sorted_simplex[0];
				DiffOp_inner[(i*2+1)*2+0] = sorted_simplex[0];
				DiffOp_inner[(i*2+0)*2+1] = sorted_simplex[1];
				DiffOp_inner[(i*2+1)*2+1] = sorted_simplex[1];

				near[1] = P_hull_coords__[4*i+0] = V_coords__[2*simplex[0]+0];
				near[2] = P_hull_coords__[4*i+1] = V_coords__[2*simplex[0]+1];
				near[3] = P_hull_coords__[4*i+2] = V_coords__[2*simplex[1]+0];
				near[4] = P_hull_coords__[4*i+3] = V_coords__[2*simplex[1]+1];

				far[1] = P_coords__[2*i+0] = 0.5 * ( P_hull_coords__[4*i+0] + P_hull_coords__[4*i+2] );
				far[2] = P_coords__[2*i+1] = 0.5 * ( P_hull_coords__[4*i+1] + P_hull_coords__[4*i+3] );

				df[0][0] = V_coords__[2*sorted_simplex[1]+0] - V_coords__[2*sorted_simplex[0]+0];
				df[1][0] = V_coords__[2*sorted_simplex[1]+1] - V_coords__[2*sorted_simplex[0]+1];

				g[0][0] = df[0][0] * df[0][0] + df[1][0] * df[1][0];

                near[0] = far[0] = sqrt( fabs(g[0][0]) );

                ginv[0][0] =  static_cast<Real>(1)/g[0][0];
                
                //  dfdagger = g^{-1} * df^T (1 x 2 matrix)
				dfdagger[0][0] = ginv[0][0] * df[0][0];
				dfdagger[0][1] = ginv[0][0] * df[1][0];

				near[5] = far[3] = static_cast<Real>(1.0) - df[0][0] * dfdagger[0][0];
				near[6] = far[4] =    - df[0][0] * dfdagger[0][1];
				near[7] = far[5] = static_cast<Real>(1.0) - df[1][0] * dfdagger[0][1];

                // derivative operator  (2 x 2 matrix)

                Real * Df = &DiffOp_value[ 4 * i ];

				Df[0] = - dfdagger[0][0];
				Df[1] =   dfdagger[0][0];
				Df[2] = - dfdagger[0][1];
				Df[3] =   dfdagger[0][1];

            }
        }
        ptoc(ClassName()+"::ComputeNearFarDataOps");
    }

	void DNearToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_near, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DNearToHulls");

        if( P_D_near.Dimension(1) != 8 )
        {
            eprint("in DNearToHulls: P_D_near.Dimension(1) != 8. Aborting");
        }

		ptr<Real> V_coords__  = V_coords.data();
		ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_near__  = P_D_near.data();
		mut<Real> buffer__    = buffer.data();
        
        if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[2*simplices__[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = V_coords__[2*simplices__[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = V_coords__[2*simplices__[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = V_coords__[2*simplices__[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = s4 + s9;
				const Real s11 = sqrt(s10);
				const Real s12 = 1/s11;
				const Real s13 = s10*s11;
				const Real s14 = 1/s13;
				const Real s15 = s3*s4;
				const Real s16 = s10*s10;
				const Real s17 = 1/s16;
				const Real s18 = 1/s10;
				const Real s19 = P_D_near__[8*i+0];
				const Real s20 = P_D_near__[8*i+1];
				const Real s21 = P_D_near__[8*i+3];
				const Real s22 = P_D_near__[8*i+4];
				const Real s23 = P_D_near__[8*i+6];
				const Real s24 = P_D_near__[8*i+2];
				const Real s25 = P_D_near__[8*i+5];
				const Real s26 = -(s18*s4);
				const Real s27 = 1 + s26;
				const Real s28 = P_D_near__[8*i+7];
				const Real s29 = s8*s9;
				const Real s30 = -(s18*s9);
				const Real s31 = 1 + s30;
				buffer__[4*i+0] += -(s12*s19*s3) - s12*s2*s21*s3 + s20*(s11 - s0*s12*s3) + s25*(-(s12*s27*s3) + s11*(-2*s15*s17 + 2*s18*s3)) - s12*s24*s3*s5 - s12*s22*s3*s7 + s23*(s12*s8 - s14*s4*s8) + s28*(-(s12*s3*s31) - 2*s14*s3*s9);
				buffer__[4*i+1] += -(s12*s19*s8) - s0*s12*s20*s8 - s12*s2*s21*s8 - s12*s22*s7*s8 + s25*(-(s12*s27*s8) - 2*s14*s4*s8) + s24*(s11 - s12*s5*s8) + s28*(-(s12*s31*s8) + s11*(-2*s17*s29 + 2*s18*s8)) + s23*(s12*s3 - s14*s3*s9);
				buffer__[4*i+2] += s12*s19*s3 + s0*s12*s20*s3 + s21*(s11 + s12*s2*s3) + s25*(s12*s27*s3 + s11*(2*s15*s17 - 2*s18*s3)) + s12*s24*s3*s5 + s12*s22*s3*s7 + s23*(-(s12*s8) + s14*s4*s8) + s28*(s12*s3*s31 + 2*s14*s3*s9);
				buffer__[4*i+3] += s12*s19*s8 + s0*s12*s20*s8 + s12*s2*s21*s8 + s12*s24*s5*s8 + s25*(s12*s27*s8 + 2*s14*s4*s8) + s22*(s11 + s12*s7*s8) + s28*(s12*s31*s8 + s11*(2*s17*s29 - 2*s18*s8)) + s23*(-(s12*s3) + s14*s3*s9);
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[2*simplices__[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = V_coords__[2*simplices__[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = V_coords__[2*simplices__[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = V_coords__[2*simplices__[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = s4 + s9;
				const Real s11 = sqrt(s10);
				const Real s12 = 1/s11;
				const Real s13 = s10*s11;
				const Real s14 = 1/s13;
				const Real s15 = s3*s4;
				const Real s16 = s10*s10;
				const Real s17 = 1/s16;
				const Real s18 = 1/s10;
				const Real s19 = P_D_near__[8*i+0];
				const Real s20 = P_D_near__[8*i+1];
				const Real s21 = P_D_near__[8*i+3];
				const Real s22 = P_D_near__[8*i+4];
				const Real s23 = P_D_near__[8*i+6];
				const Real s24 = P_D_near__[8*i+2];
				const Real s25 = P_D_near__[8*i+5];
				const Real s26 = -(s18*s4);
				const Real s27 = 1 + s26;
				const Real s28 = P_D_near__[8*i+7];
				const Real s29 = s8*s9;
				const Real s30 = -(s18*s9);
				const Real s31 = 1 + s30;
				buffer__[4*i+0] = -(s12*s19*s3) - s12*s2*s21*s3 + s20*(s11 - s0*s12*s3) + s25*(-(s12*s27*s3) + s11*(-2*s15*s17 + 2*s18*s3)) - s12*s24*s3*s5 - s12*s22*s3*s7 + s23*(s12*s8 - s14*s4*s8) + s28*(-(s12*s3*s31) - 2*s14*s3*s9);
				buffer__[4*i+1] = -(s12*s19*s8) - s0*s12*s20*s8 - s12*s2*s21*s8 - s12*s22*s7*s8 + s25*(-(s12*s27*s8) - 2*s14*s4*s8) + s24*(s11 - s12*s5*s8) + s28*(-(s12*s31*s8) + s11*(-2*s17*s29 + 2*s18*s8)) + s23*(s12*s3 - s14*s3*s9);
				buffer__[4*i+2] = s12*s19*s3 + s0*s12*s20*s3 + s21*(s11 + s12*s2*s3) + s25*(s12*s27*s3 + s11*(2*s15*s17 - 2*s18*s3)) + s12*s24*s3*s5 + s12*s22*s3*s7 + s23*(-(s12*s8) + s14*s4*s8) + s28*(s12*s3*s31 + 2*s14*s3*s9);
				buffer__[4*i+3] = s12*s19*s8 + s0*s12*s20*s8 + s12*s2*s21*s8 + s12*s24*s5*s8 + s25*(s12*s27*s8 + 2*s14*s4*s8) + s22*(s11 + s12*s7*s8) + s28*(s12*s31*s8 + s11*(2*s17*s29 - 2*s18*s8)) + s23*(-(s12*s3) + s14*s3*s9);
			}
		}

        ptoc(ClassName()+"::DNearToHulls");
        
    }

	void DFarToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_far, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DFarToHulls");

        if( P_D_far.Dimension(1) != 6 )
        {
            eprint("in DFarToHulls: P_D_far.Dimension(1) != 6. Aborting");
        }

		ptr<Real> V_coords__  = V_coords.data();
		ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_far__   = P_D_far.data();
		mut<Real> buffer__    = buffer.data();
        
        if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[2*simplices__[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = V_coords__[2*simplices__[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = V_coords__[2*simplices__[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = V_coords__[2*simplices__[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = s4 + s9;
				const Real s11 = sqrt(s10);
				const Real s12 = 1/s11;
				const Real s13 = s10*s11;
				const Real s14 = 1/s13;
				const Real s15 = s3*s4;
				const Real s16 = s10*s10;
				const Real s17 = 1/s16;
				const Real s18 = 1/s10;
				const Real s19 = P_D_far__[6*i+0];
				const Real s20 = P_D_far__[6*i+1];
				const Real s21 = s0 + s2;
				const Real s22 = P_D_far__[6*i+4];
				const Real s23 = P_D_far__[6*i+2];
				const Real s24 = s5 + s7;
				const Real s25 = s11/2.;
				const Real s26 = P_D_far__[6*i+3];
				const Real s27 = -(s18*s4);
				const Real s28 = 1 + s27;
				const Real s29 = P_D_far__[6*i+5];
				const Real s30 = s8*s9;
				const Real s31 = -(s18*s9);
				const Real s32 = 1 + s31;
				buffer__[4*i+0] += -(s12*s19*s3) - (s12*s23*s24*s3)/2. + s20*(s25 - (s12*s21*s3)/2.) + s26*(-(s12*s28*s3) + s11*(-2*s15*s17 + 2*s18*s3)) + s22*(s12*s8 - s14*s4*s8) + s29*(-(s12*s3*s32) - 2*s14*s3*s9);
				buffer__[4*i+1] += -(s12*s19*s8) - (s12*s20*s21*s8)/2. + s23*(s25 - (s12*s24*s8)/2.) + s26*(-(s12*s28*s8) - 2*s14*s4*s8) + s29*(-(s12*s32*s8) + s11*(-2*s17*s30 + 2*s18*s8)) + s22*(s12*s3 - s14*s3*s9);
				buffer__[4*i+2] += s12*s19*s3 + (s12*s23*s24*s3)/2. + s20*(s25 + (s12*s21*s3)/2.) + s26*(s12*s28*s3 + s11*(2*s15*s17 - 2*s18*s3)) + s22*(-(s12*s8) + s14*s4*s8) + s29*(s12*s3*s32 + 2*s14*s3*s9);
				buffer__[4*i+3] += s12*s19*s8 + (s12*s20*s21*s8)/2. + s23*(s25 + (s12*s24*s8)/2.) + s26*(s12*s28*s8 + 2*s14*s4*s8) + s29*(s12*s32*s8 + s11*(2*s17*s30 - 2*s18*s8)) + s22*(-(s12*s3) + s14*s3*s9);
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[2*simplices__[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = V_coords__[2*simplices__[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = V_coords__[2*simplices__[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = V_coords__[2*simplices__[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = s4 + s9;
				const Real s11 = sqrt(s10);
				const Real s12 = 1/s11;
				const Real s13 = s10*s11;
				const Real s14 = 1/s13;
				const Real s15 = s3*s4;
				const Real s16 = s10*s10;
				const Real s17 = 1/s16;
				const Real s18 = 1/s10;
				const Real s19 = P_D_far__[6*i+0];
				const Real s20 = P_D_far__[6*i+1];
				const Real s21 = s0 + s2;
				const Real s22 = P_D_far__[6*i+4];
				const Real s23 = P_D_far__[6*i+2];
				const Real s24 = s5 + s7;
				const Real s25 = s11/2.;
				const Real s26 = P_D_far__[6*i+3];
				const Real s27 = -(s18*s4);
				const Real s28 = 1 + s27;
				const Real s29 = P_D_far__[6*i+5];
				const Real s30 = s8*s9;
				const Real s31 = -(s18*s9);
				const Real s32 = 1 + s31;
				buffer__[4*i+0] = -(s12*s19*s3) - (s12*s23*s24*s3)/2. + s20*(s25 - (s12*s21*s3)/2.) + s26*(-(s12*s28*s3) + s11*(-2*s15*s17 + 2*s18*s3)) + s22*(s12*s8 - s14*s4*s8) + s29*(-(s12*s3*s32) - 2*s14*s3*s9);
				buffer__[4*i+1] = -(s12*s19*s8) - (s12*s20*s21*s8)/2. + s23*(s25 - (s12*s24*s8)/2.) + s26*(-(s12*s28*s8) - 2*s14*s4*s8) + s29*(-(s12*s32*s8) + s11*(-2*s17*s30 + 2*s18*s8)) + s22*(s12*s3 - s14*s3*s9);
				buffer__[4*i+2] = s12*s19*s3 + (s12*s23*s24*s3)/2. + s20*(s25 + (s12*s21*s3)/2.) + s26*(s12*s28*s3 + s11*(2*s15*s17 - 2*s18*s3)) + s22*(-(s12*s8) + s14*s4*s8) + s29*(s12*s3*s32 + 2*s14*s3*s9);
				buffer__[4*i+3] = s12*s19*s8 + (s12*s20*s21*s8)/2. + s23*(s25 + (s12*s24*s8)/2.) + s26*(s12*s28*s8 + 2*s14*s4*s8) + s29*(s12*s32*s8 + s11*(2*s17*s30 - 2*s18*s8)) + s22*(-(s12*s3) + s14*s3*s9);
			}
		}

        ptoc(ClassName()+"::DFarToHulls");
        
    }

	}; // SimplicialMeshDetails<1,2,Real,Int>

//----------------------------------------------------------------------------------------------

    template<typename Real, typename Int>
    struct SimplicialMeshDetails<1,3,Real,Int>
    {
	private:

		const Int thread_count = 1;

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(std::max(static_cast<Int>(1),thread_count_))
		{}
	
		inline Int ThreadCount() const
		{
			return thread_count;
		}
	
		inline std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<1,3,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void ComputeNearFarData( 
		const Tensor2<Real,Int> & V_coords,
		const Tensor2<Int ,Int> & simplices,
			  Tensor2<Real,Int> & P_near,
			  Tensor2<Real,Int> & P_far
	) const
    {
        ptic(ClassName()+"::ComputeNearFarData");
        
        //Int size       = 2;
        //Int amb_dim    = 3;
        //Int dom_dim    = 1;

        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			ptr<Real> V_coords__      = V_coords.data();	
			const Int  * restrict const simplices__     = simplices.data();

			Real hull    [2][3];
			Real df      [3][1];
			Real dfdagger[1][3];
			Real g       [1][1];
			Real ginv    [1][1];

			Int simplex  [2];
			
			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {
				Real * restrict const near = P_near.data(i);                    
				Real * restrict const far  = P_far.data(i);   
            
				simplex[0] = simplices__[2*i +0];
				simplex[1] = simplices__[2*i +1];

				near[1] = hull[0][0] = V_coords__[3*simplex[0]+0];
				near[2] = hull[0][1] = V_coords__[3*simplex[0]+1];
				near[3] = hull[0][2] = V_coords__[3*simplex[0]+2];
				near[4] = hull[1][0] = V_coords__[3*simplex[1]+0];
				near[5] = hull[1][1] = V_coords__[3*simplex[1]+1];
				near[6] = hull[1][2] = V_coords__[3*simplex[1]+2];

				far[1] = static_cast<Real>(0.5) * ( hull[0][0] + hull[1][0] );
				far[2] = static_cast<Real>(0.5) * ( hull[0][1] + hull[1][1] );
				far[3] = static_cast<Real>(0.5) * ( hull[0][2] + hull[1][2] );

				df[0][0] = hull[1][0] - hull[0][0];
				df[1][0] = hull[1][1] - hull[0][1];
				df[2][0] = hull[1][2] - hull[0][2];

				g[0][0] = df[0][0] * df[0][0] + df[1][0] * df[1][0] + df[2][0] * df[2][0];


				near[0] = far[0] = sqrt( fabs(g[0][0]) );

                ginv[0][0] = static_cast<Real>(1)/g[0][0];
                
                //  dfdagger = g^{-1} * df^T (1 x 3 matrix)
				dfdagger[0][0] = ginv[0][0] * df[0][0];
				dfdagger[0][1] = ginv[0][0] * df[1][0];
				dfdagger[0][2] = ginv[0][0] * df[2][0];
            
				near[ 7] = far[ 4]  = static_cast<Real>(1) - df[0][0] * dfdagger[0][0];
				near[ 8] = far[ 5]  =    - df[0][0] * dfdagger[0][1];
				near[ 9] = far[ 6]  =    - df[0][0] * dfdagger[0][2];
				near[10] = far[ 7]  = static_cast<Real>(1) - df[1][0] * dfdagger[0][1];
				near[11] = far[ 8]  =    - df[1][0] * dfdagger[0][2];
				near[12] = far[ 9]  = static_cast<Real>(1) - df[2][0] * dfdagger[0][2];

            } // for( Int i = i_begin; i < i_end; ++i )

		} // #pragma omp parallel num_threads( ThreadCount() )

        ptoc(ClassName()+"::ComputeNearFarData");
    }

    void ComputeNearFarDataOps( 
		const Tensor2<Real,Int> & V_coords,
        const Tensor2<Int ,Int> & simplices,
		      Tensor2<Real,Int> & P_coords,
		      Tensor3<Real,Int> & P_hull_coords,
		      Tensor2<Real,Int> & P_near,
		      Tensor2<Real,Int> & P_far,
		Sparse::MatrixCSR<Real,Int,Int> & DiffOp,
		Sparse::MatrixCSR<Real,Int,Int> & AvOp 
	) const
    {
        ptic(ClassName()+"::ComputeNearFarDataOps");
        
        //Int size       = 2;
        //Int amb_dim    = 3;
        //Int dom_dim    = 1;
        
        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			mut<Int>  AvOp_outer = AvOp.Outer().data();
			mut<Int>  AvOp_inner = AvOp.Inner().data();
			mut<Real> AvOp_value = AvOp.Values().data();

			mut<Int>  DiffOp_outer = DiffOp.Outer().data();
			mut<Int>  DiffOp_inner = DiffOp.Inner().data();
			mut<Real> DiffOp_value = DiffOp.Value().data();

			ptr<Real> V_coords__      = V_coords.data();
			
			ptr<Int>  simplices__     = simplices.data();
		    mut<Real> P_hull_coords__ = P_hull_coords.data();
			mut<Real> P_coords__      = P_coords.data();

			Real df       [3][1];
			Real dfdagger [1][3];
			Real g        [1][1];
			Real ginv     [1][1];

			Int simplex        [2];
			Int sorted_simplex [2];

			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {

				mut<Real> near = P_near.data(i);                    
				mut<Real> far  = P_far.data(i);

				simplex[0] = sorted_simplex[0] = simplices__[2*i +0];
				simplex[1] = sorted_simplex[1] = simplices__[2*i +1];
                  
                // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
                std::sort( sorted_simplex, sorted_simplex + 2 );

				AvOp_outer[i+1] = (i+1) * 2;  
                      
				AvOp_inner[2*i+0] = sorted_simplex[0];
				AvOp_inner[2*i+1] = sorted_simplex[1];

				AvOp_value[2*i+0] = 0.5;
				AvOp_value[2*i+1] = 0.5;

				DiffOp_outer[3*i+0] = (3 * i + 0) * 2;
				DiffOp_outer[3*i+1] = (3 * i + 1) * 2;
				DiffOp_outer[3*i+2] = (3 * i + 2) * 2;

				DiffOp_inner[(i*3+0)*2+0] = sorted_simplex[0];
				DiffOp_inner[(i*3+1)*2+0] = sorted_simplex[0];
				DiffOp_inner[(i*3+2)*2+0] = sorted_simplex[0];
				DiffOp_inner[(i*3+0)*2+1] = sorted_simplex[1];
				DiffOp_inner[(i*3+1)*2+1] = sorted_simplex[1];
				DiffOp_inner[(i*3+2)*2+1] = sorted_simplex[1];

				near[1] = P_hull_coords__[6*i+0] = V_coords__[3*simplex[0]+0];
				near[2] = P_hull_coords__[6*i+1] = V_coords__[3*simplex[0]+1];
				near[3] = P_hull_coords__[6*i+2] = V_coords__[3*simplex[0]+2];
				near[4] = P_hull_coords__[6*i+3] = V_coords__[3*simplex[1]+0];
				near[5] = P_hull_coords__[6*i+4] = V_coords__[3*simplex[1]+1];
				near[6] = P_hull_coords__[6*i+5] = V_coords__[3*simplex[1]+2];

				far[1] = P_coords__[3*i+0] = 0.5 * ( P_hull_coords__[6*i+0] + P_hull_coords__[6*i+3] );
				far[2] = P_coords__[3*i+1] = 0.5 * ( P_hull_coords__[6*i+1] + P_hull_coords__[6*i+4] );
				far[3] = P_coords__[3*i+2] = 0.5 * ( P_hull_coords__[6*i+2] + P_hull_coords__[6*i+5] );

				df[0][0] = V_coords__[3*sorted_simplex[1]+0] - V_coords__[3*sorted_simplex[0]+0];
				df[1][0] = V_coords__[3*sorted_simplex[1]+1] - V_coords__[3*sorted_simplex[0]+1];
				df[2][0] = V_coords__[3*sorted_simplex[1]+2] - V_coords__[3*sorted_simplex[0]+2];

				g[0][0] = df[0][0] * df[0][0] + df[1][0] * df[1][0] + df[2][0] * df[2][0];

                near[0] = far[0] = sqrt( fabs(g[0][0]) );

                ginv[0][0] =  static_cast<Real>(1)/g[0][0];
                
                //  dfdagger = g^{-1} * df^T (1 x 3 matrix)
				dfdagger[0][0] = ginv[0][0] * df[0][0];
				dfdagger[0][1] = ginv[0][0] * df[1][0];
				dfdagger[0][2] = ginv[0][0] * df[2][0];

				near[ 7] = far[ 4] = static_cast<Real>(1.0) - df[0][0] * dfdagger[0][0];
				near[ 8] = far[ 5] =    - df[0][0] * dfdagger[0][1];
				near[ 9] = far[ 6] =    - df[0][0] * dfdagger[0][2];
				near[10] = far[ 7] = static_cast<Real>(1.0) - df[1][0] * dfdagger[0][1];
				near[11] = far[ 8] =    - df[1][0] * dfdagger[0][2];
				near[12] = far[ 9] = static_cast<Real>(1.0) - df[2][0] * dfdagger[0][2];

                // derivative operator  (3 x 2 matrix)

                Real * Df = &DiffOp_value[ 6 * i ];

				Df[ 0] = - dfdagger[0][0];
				Df[ 1] =   dfdagger[0][0];
				Df[ 2] = - dfdagger[0][1];
				Df[ 3] =   dfdagger[0][1];
				Df[ 4] = - dfdagger[0][2];
				Df[ 5] =   dfdagger[0][2];

            }
        }
        ptoc(ClassName()+"::ComputeNearFarDataOps");
    }

	void DNearToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_near, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DNearToHulls");

        if( P_D_near.Dimension(1) != 13 )
        {
            eprint("in DNearToHulls: P_D_near.Dimension(1) != 13. Aborting");
        }

		ptr<Real> V_coords__  = V_coords.data();
		ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_near__  = P_D_near.data();
		mut<Real> buffer__    = buffer.data();
        
        if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[3*simplices__[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = V_coords__[3*simplices__[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = V_coords__[3*simplices__[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = V_coords__[3*simplices__[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = V_coords__[3*simplices__[2*i+0]+2];
				const Real s11 = -s10;
				const Real s12 = V_coords__[3*simplices__[2*i+1]+2];
				const Real s13 = s11 + s12;
				const Real s14 = s13*s13;
				const Real s15 = s14 + s4 + s9;
				const Real s16 = sqrt(s15);
				const Real s17 = s15*s16;
				const Real s18 = 1/s17;
				const Real s19 = 1/s16;
				const Real s20 = s3*s4;
				const Real s21 = s15*s15;
				const Real s22 = 1/s21;
				const Real s23 = 1/s15;
				const Real s24 = P_D_near__[13*i+9];
				const Real s25 = P_D_near__[13*i+0];
				const Real s26 = P_D_near__[13*i+1];
				const Real s27 = P_D_near__[13*i+3];
				const Real s28 = P_D_near__[13*i+4];
				const Real s29 = P_D_near__[13*i+5];
				const Real s30 = P_D_near__[13*i+6];
				const Real s31 = P_D_near__[13*i+8];
				const Real s32 = P_D_near__[13*i+11];
				const Real s33 = s13*s19;
				const Real s34 = P_D_near__[13*i+2];
				const Real s35 = P_D_near__[13*i+7];
				const Real s36 = -(s23*s4);
				const Real s37 = 1 + s36;
				const Real s38 = P_D_near__[13*i+10];
				const Real s39 = s8*s9;
				const Real s40 = -(s23*s9);
				const Real s41 = 1 + s40;
				const Real s42 = P_D_near__[13*i+12];
				const Real s43 = -(s14*s23);
				const Real s44 = 1 + s43;
				const Real s45 = s19*s3;
				const Real s46 = s19*s8;
				const Real s47 = s13*s14;
				const Real s48 = -(s13*s19);
				const Real s49 = -(s19*s3);
				const Real s50 = -(s19*s8);
				buffer__[6*i+0] += -(s19*s25*s3) - s10*s19*s27*s3 - s19*s2*s28*s3 + s26*(s16 - s0*s19*s3) - s12*s19*s3*s30 + s35*(s16*(-2*s20*s22 + 2*s23*s3) - s19*s3*s37) + s24*(s33 - s13*s18*s4) + s42*(-2*s14*s18*s3 - s19*s3*s44) - s19*s3*s34*s5 - s19*s29*s3*s7 - s13*s18*s3*s32*s8 + s31*(s46 - s18*s4*s8) + s38*(-(s19*s3*s41) - 2*s18*s3*s9);
				buffer__[6*i+1] += -(s19*s25*s8) - s0*s19*s26*s8 - s10*s19*s27*s8 - s19*s2*s28*s8 - s13*s18*s24*s3*s8 - s12*s19*s30*s8 - s19*s29*s7*s8 + s35*(-(s19*s37*s8) - 2*s18*s4*s8) + s42*(-2*s14*s18*s8 - s19*s44*s8) + s34*(s16 - s19*s5*s8) + s38*(-(s19*s41*s8) + s16*(-2*s22*s39 + 2*s23*s8)) + s32*(s33 - s13*s18*s9) + s31*(s45 - s18*s3*s9);
				buffer__[6*i+2] += -(s13*s19*s25) - s0*s13*s19*s26 + (s16 - s10*s13*s19)*s27 - s13*s19*s2*s28 - s12*s13*s19*s30 + s35*(-(s13*s19*s37) - 2*s13*s18*s4) + s24*(-(s14*s18*s3) + s45) + s42*(-(s13*s19*s44) + s16*(2*s13*s23 - 2*s22*s47)) - s13*s19*s34*s5 - s13*s19*s29*s7 - s13*s18*s3*s31*s8 + s32*(s46 - s14*s18*s8) + s38*(-(s13*s19*s41) - 2*s13*s18*s9);
				buffer__[6*i+3] += s19*s25*s3 + s0*s19*s26*s3 + s10*s19*s27*s3 + s28*(s16 + s19*s2*s3) + s12*s19*s3*s30 + s35*(s16*(2*s20*s22 - 2*s23*s3) + s19*s3*s37) + s42*(2*s14*s18*s3 + s19*s3*s44) + s24*(s13*s18*s4 + s48) + s19*s3*s34*s5 + s19*s29*s3*s7 + s13*s18*s3*s32*s8 + s31*(s50 + s18*s4*s8) + s38*(s19*s3*s41 + 2*s18*s3*s9);
				buffer__[6*i+4] += s19*s25*s8 + s0*s19*s26*s8 + s10*s19*s27*s8 + s19*s2*s28*s8 + s13*s18*s24*s3*s8 + s12*s19*s30*s8 + s19*s34*s5*s8 + s35*(s19*s37*s8 + 2*s18*s4*s8) + s42*(2*s14*s18*s8 + s19*s44*s8) + s29*(s16 + s19*s7*s8) + s38*(s19*s41*s8 + s16*(2*s22*s39 - 2*s23*s8)) + s32*(s48 + s13*s18*s9) + s31*(s49 + s18*s3*s9);
				buffer__[6*i+5] += s13*s19*s25 + s0*s13*s19*s26 + s10*s13*s19*s27 + s13*s19*s2*s28 + (s16 + s12*s13*s19)*s30 + s35*(s13*s19*s37 + 2*s13*s18*s4) + s42*(s13*s19*s44 + s16*(-2*s13*s23 + 2*s22*s47)) + s24*(s14*s18*s3 + s49) + s13*s19*s34*s5 + s13*s19*s29*s7 + s13*s18*s3*s31*s8 + s32*(s50 + s14*s18*s8) + s38*(s13*s19*s41 + 2*s13*s18*s9);
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[3*simplices__[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = V_coords__[3*simplices__[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = V_coords__[3*simplices__[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = V_coords__[3*simplices__[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = V_coords__[3*simplices__[2*i+0]+2];
				const Real s11 = -s10;
				const Real s12 = V_coords__[3*simplices__[2*i+1]+2];
				const Real s13 = s11 + s12;
				const Real s14 = s13*s13;
				const Real s15 = s14 + s4 + s9;
				const Real s16 = sqrt(s15);
				const Real s17 = s15*s16;
				const Real s18 = 1/s17;
				const Real s19 = 1/s16;
				const Real s20 = s3*s4;
				const Real s21 = s15*s15;
				const Real s22 = 1/s21;
				const Real s23 = 1/s15;
				const Real s24 = P_D_near__[13*i+9];
				const Real s25 = P_D_near__[13*i+0];
				const Real s26 = P_D_near__[13*i+1];
				const Real s27 = P_D_near__[13*i+3];
				const Real s28 = P_D_near__[13*i+4];
				const Real s29 = P_D_near__[13*i+5];
				const Real s30 = P_D_near__[13*i+6];
				const Real s31 = P_D_near__[13*i+8];
				const Real s32 = P_D_near__[13*i+11];
				const Real s33 = s13*s19;
				const Real s34 = P_D_near__[13*i+2];
				const Real s35 = P_D_near__[13*i+7];
				const Real s36 = -(s23*s4);
				const Real s37 = 1 + s36;
				const Real s38 = P_D_near__[13*i+10];
				const Real s39 = s8*s9;
				const Real s40 = -(s23*s9);
				const Real s41 = 1 + s40;
				const Real s42 = P_D_near__[13*i+12];
				const Real s43 = -(s14*s23);
				const Real s44 = 1 + s43;
				const Real s45 = s19*s3;
				const Real s46 = s19*s8;
				const Real s47 = s13*s14;
				const Real s48 = -(s13*s19);
				const Real s49 = -(s19*s3);
				const Real s50 = -(s19*s8);
				buffer__[6*i+0] = -(s19*s25*s3) - s10*s19*s27*s3 - s19*s2*s28*s3 + s26*(s16 - s0*s19*s3) - s12*s19*s3*s30 + s35*(s16*(-2*s20*s22 + 2*s23*s3) - s19*s3*s37) + s24*(s33 - s13*s18*s4) + s42*(-2*s14*s18*s3 - s19*s3*s44) - s19*s3*s34*s5 - s19*s29*s3*s7 - s13*s18*s3*s32*s8 + s31*(s46 - s18*s4*s8) + s38*(-(s19*s3*s41) - 2*s18*s3*s9);
				buffer__[6*i+1] = -(s19*s25*s8) - s0*s19*s26*s8 - s10*s19*s27*s8 - s19*s2*s28*s8 - s13*s18*s24*s3*s8 - s12*s19*s30*s8 - s19*s29*s7*s8 + s35*(-(s19*s37*s8) - 2*s18*s4*s8) + s42*(-2*s14*s18*s8 - s19*s44*s8) + s34*(s16 - s19*s5*s8) + s38*(-(s19*s41*s8) + s16*(-2*s22*s39 + 2*s23*s8)) + s32*(s33 - s13*s18*s9) + s31*(s45 - s18*s3*s9);
				buffer__[6*i+2] = -(s13*s19*s25) - s0*s13*s19*s26 + (s16 - s10*s13*s19)*s27 - s13*s19*s2*s28 - s12*s13*s19*s30 + s35*(-(s13*s19*s37) - 2*s13*s18*s4) + s24*(-(s14*s18*s3) + s45) + s42*(-(s13*s19*s44) + s16*(2*s13*s23 - 2*s22*s47)) - s13*s19*s34*s5 - s13*s19*s29*s7 - s13*s18*s3*s31*s8 + s32*(s46 - s14*s18*s8) + s38*(-(s13*s19*s41) - 2*s13*s18*s9);
				buffer__[6*i+3] = s19*s25*s3 + s0*s19*s26*s3 + s10*s19*s27*s3 + s28*(s16 + s19*s2*s3) + s12*s19*s3*s30 + s35*(s16*(2*s20*s22 - 2*s23*s3) + s19*s3*s37) + s42*(2*s14*s18*s3 + s19*s3*s44) + s24*(s13*s18*s4 + s48) + s19*s3*s34*s5 + s19*s29*s3*s7 + s13*s18*s3*s32*s8 + s31*(s50 + s18*s4*s8) + s38*(s19*s3*s41 + 2*s18*s3*s9);
				buffer__[6*i+4] = s19*s25*s8 + s0*s19*s26*s8 + s10*s19*s27*s8 + s19*s2*s28*s8 + s13*s18*s24*s3*s8 + s12*s19*s30*s8 + s19*s34*s5*s8 + s35*(s19*s37*s8 + 2*s18*s4*s8) + s42*(2*s14*s18*s8 + s19*s44*s8) + s29*(s16 + s19*s7*s8) + s38*(s19*s41*s8 + s16*(2*s22*s39 - 2*s23*s8)) + s32*(s48 + s13*s18*s9) + s31*(s49 + s18*s3*s9);
				buffer__[6*i+5] = s13*s19*s25 + s0*s13*s19*s26 + s10*s13*s19*s27 + s13*s19*s2*s28 + (s16 + s12*s13*s19)*s30 + s35*(s13*s19*s37 + 2*s13*s18*s4) + s42*(s13*s19*s44 + s16*(-2*s13*s23 + 2*s22*s47)) + s24*(s14*s18*s3 + s49) + s13*s19*s34*s5 + s13*s19*s29*s7 + s13*s18*s3*s31*s8 + s32*(s50 + s14*s18*s8) + s38*(s13*s19*s41 + 2*s13*s18*s9);
			}
		}

        ptoc(ClassName()+"::DNearToHulls");
        
    }

	void DFarToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_far, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DFarToHulls");

        if( P_D_far.Dimension(1) != 10 )
        {
            eprint("in DFarToHulls: P_D_far.Dimension(1) != 10. Aborting");
        }

		ptr<Real> V_coords__  = V_coords.data();
		ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_far__   = P_D_far.data();
		mut<Real> buffer__    = buffer.data();
        
        if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[3*simplices__[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = V_coords__[3*simplices__[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = V_coords__[3*simplices__[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = V_coords__[3*simplices__[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = V_coords__[3*simplices__[2*i+0]+2];
				const Real s11 = -s10;
				const Real s12 = V_coords__[3*simplices__[2*i+1]+2];
				const Real s13 = s11 + s12;
				const Real s14 = s13*s13;
				const Real s15 = s14 + s4 + s9;
				const Real s16 = sqrt(s15);
				const Real s17 = s15*s16;
				const Real s18 = 1/s17;
				const Real s19 = 1/s16;
				const Real s20 = s3*s4;
				const Real s21 = s15*s15;
				const Real s22 = 1/s21;
				const Real s23 = 1/s15;
				const Real s24 = P_D_far__[10*i+6];
				const Real s25 = P_D_far__[10*i+0];
				const Real s26 = P_D_far__[10*i+1];
				const Real s27 = s0 + s2;
				const Real s28 = P_D_far__[10*i+3];
				const Real s29 = s10 + s12;
				const Real s30 = P_D_far__[10*i+5];
				const Real s31 = P_D_far__[10*i+8];
				const Real s32 = s13*s19;
				const Real s33 = P_D_far__[10*i+2];
				const Real s34 = s5 + s7;
				const Real s35 = s16/2.;
				const Real s36 = P_D_far__[10*i+4];
				const Real s37 = -(s23*s4);
				const Real s38 = 1 + s37;
				const Real s39 = P_D_far__[10*i+7];
				const Real s40 = s8*s9;
				const Real s41 = -(s23*s9);
				const Real s42 = 1 + s41;
				const Real s43 = P_D_far__[10*i+9];
				const Real s44 = -(s14*s23);
				const Real s45 = 1 + s44;
				const Real s46 = s19*s3;
				const Real s47 = s19*s8;
				const Real s48 = s13*s14;
				const Real s49 = -(s13*s19);
				const Real s50 = -(s19*s3);
				const Real s51 = -(s19*s8);
				buffer__[6*i+0] += -(s19*s25*s3) - (s19*s28*s29*s3)/2. - (s19*s3*s33*s34)/2. + s26*(-0.5*(s19*s27*s3) + s35) + s36*(s16*(-2*s20*s22 + 2*s23*s3) - s19*s3*s38) + s24*(s32 - s13*s18*s4) + s43*(-2*s14*s18*s3 - s19*s3*s45) - s13*s18*s3*s31*s8 + s30*(s47 - s18*s4*s8) + s39*(-(s19*s3*s42) - 2*s18*s3*s9);
				buffer__[6*i+1] += -(s19*s25*s8) - (s19*s26*s27*s8)/2. - (s19*s28*s29*s8)/2. - s13*s18*s24*s3*s8 + s33*(s35 - (s19*s34*s8)/2.) + s36*(-(s19*s38*s8) - 2*s18*s4*s8) + s43*(-2*s14*s18*s8 - s19*s45*s8) + s39*(-(s19*s42*s8) + s16*(-2*s22*s40 + 2*s23*s8)) + s31*(s32 - s13*s18*s9) + s30*(s46 - s18*s3*s9);
				buffer__[6*i+2] += -(s13*s19*s25) - (s13*s19*s26*s27)/2. - (s13*s19*s33*s34)/2. + s28*(-0.5*(s13*s19*s29) + s35) + s36*(-(s13*s19*s38) - 2*s13*s18*s4) + s24*(-(s14*s18*s3) + s46) + s43*(-(s13*s19*s45) + s16*(2*s13*s23 - 2*s22*s48)) - s13*s18*s3*s30*s8 + s31*(s47 - s14*s18*s8) + s39*(-(s13*s19*s42) - 2*s13*s18*s9);
				buffer__[6*i+3] += s19*s25*s3 + (s19*s28*s29*s3)/2. + (s19*s3*s33*s34)/2. + s26*((s19*s27*s3)/2. + s35) + s36*(s16*(2*s20*s22 - 2*s23*s3) + s19*s3*s38) + s43*(2*s14*s18*s3 + s19*s3*s45) + s24*(s13*s18*s4 + s49) + s13*s18*s3*s31*s8 + s30*(s51 + s18*s4*s8) + s39*(s19*s3*s42 + 2*s18*s3*s9);
				buffer__[6*i+4] += s19*s25*s8 + (s19*s26*s27*s8)/2. + (s19*s28*s29*s8)/2. + s13*s18*s24*s3*s8 + s33*(s35 + (s19*s34*s8)/2.) + s36*(s19*s38*s8 + 2*s18*s4*s8) + s43*(2*s14*s18*s8 + s19*s45*s8) + s39*(s19*s42*s8 + s16*(2*s22*s40 - 2*s23*s8)) + s31*(s49 + s13*s18*s9) + s30*(s50 + s18*s3*s9);
				buffer__[6*i+5] += s13*s19*s25 + (s13*s19*s26*s27)/2. + (s13*s19*s33*s34)/2. + s28*((s13*s19*s29)/2. + s35) + s36*(s13*s19*s38 + 2*s13*s18*s4) + s43*(s13*s19*s45 + s16*(-2*s13*s23 + 2*s22*s48)) + s24*(s14*s18*s3 + s50) + s13*s18*s3*s30*s8 + s31*(s51 + s14*s18*s8) + s39*(s13*s19*s42 + 2*s13*s18*s9);
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[3*simplices__[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = V_coords__[3*simplices__[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = V_coords__[3*simplices__[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = V_coords__[3*simplices__[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = V_coords__[3*simplices__[2*i+0]+2];
				const Real s11 = -s10;
				const Real s12 = V_coords__[3*simplices__[2*i+1]+2];
				const Real s13 = s11 + s12;
				const Real s14 = s13*s13;
				const Real s15 = s14 + s4 + s9;
				const Real s16 = sqrt(s15);
				const Real s17 = s15*s16;
				const Real s18 = 1/s17;
				const Real s19 = 1/s16;
				const Real s20 = s3*s4;
				const Real s21 = s15*s15;
				const Real s22 = 1/s21;
				const Real s23 = 1/s15;
				const Real s24 = P_D_far__[10*i+6];
				const Real s25 = P_D_far__[10*i+0];
				const Real s26 = P_D_far__[10*i+1];
				const Real s27 = s0 + s2;
				const Real s28 = P_D_far__[10*i+3];
				const Real s29 = s10 + s12;
				const Real s30 = P_D_far__[10*i+5];
				const Real s31 = P_D_far__[10*i+8];
				const Real s32 = s13*s19;
				const Real s33 = P_D_far__[10*i+2];
				const Real s34 = s5 + s7;
				const Real s35 = s16/2.;
				const Real s36 = P_D_far__[10*i+4];
				const Real s37 = -(s23*s4);
				const Real s38 = 1 + s37;
				const Real s39 = P_D_far__[10*i+7];
				const Real s40 = s8*s9;
				const Real s41 = -(s23*s9);
				const Real s42 = 1 + s41;
				const Real s43 = P_D_far__[10*i+9];
				const Real s44 = -(s14*s23);
				const Real s45 = 1 + s44;
				const Real s46 = s19*s3;
				const Real s47 = s19*s8;
				const Real s48 = s13*s14;
				const Real s49 = -(s13*s19);
				const Real s50 = -(s19*s3);
				const Real s51 = -(s19*s8);
				buffer__[6*i+0] = -(s19*s25*s3) - (s19*s28*s29*s3)/2. - (s19*s3*s33*s34)/2. + s26*(-0.5*(s19*s27*s3) + s35) + s36*(s16*(-2*s20*s22 + 2*s23*s3) - s19*s3*s38) + s24*(s32 - s13*s18*s4) + s43*(-2*s14*s18*s3 - s19*s3*s45) - s13*s18*s3*s31*s8 + s30*(s47 - s18*s4*s8) + s39*(-(s19*s3*s42) - 2*s18*s3*s9);
				buffer__[6*i+1] = -(s19*s25*s8) - (s19*s26*s27*s8)/2. - (s19*s28*s29*s8)/2. - s13*s18*s24*s3*s8 + s33*(s35 - (s19*s34*s8)/2.) + s36*(-(s19*s38*s8) - 2*s18*s4*s8) + s43*(-2*s14*s18*s8 - s19*s45*s8) + s39*(-(s19*s42*s8) + s16*(-2*s22*s40 + 2*s23*s8)) + s31*(s32 - s13*s18*s9) + s30*(s46 - s18*s3*s9);
				buffer__[6*i+2] = -(s13*s19*s25) - (s13*s19*s26*s27)/2. - (s13*s19*s33*s34)/2. + s28*(-0.5*(s13*s19*s29) + s35) + s36*(-(s13*s19*s38) - 2*s13*s18*s4) + s24*(-(s14*s18*s3) + s46) + s43*(-(s13*s19*s45) + s16*(2*s13*s23 - 2*s22*s48)) - s13*s18*s3*s30*s8 + s31*(s47 - s14*s18*s8) + s39*(-(s13*s19*s42) - 2*s13*s18*s9);
				buffer__[6*i+3] = s19*s25*s3 + (s19*s28*s29*s3)/2. + (s19*s3*s33*s34)/2. + s26*((s19*s27*s3)/2. + s35) + s36*(s16*(2*s20*s22 - 2*s23*s3) + s19*s3*s38) + s43*(2*s14*s18*s3 + s19*s3*s45) + s24*(s13*s18*s4 + s49) + s13*s18*s3*s31*s8 + s30*(s51 + s18*s4*s8) + s39*(s19*s3*s42 + 2*s18*s3*s9);
				buffer__[6*i+4] = s19*s25*s8 + (s19*s26*s27*s8)/2. + (s19*s28*s29*s8)/2. + s13*s18*s24*s3*s8 + s33*(s35 + (s19*s34*s8)/2.) + s36*(s19*s38*s8 + 2*s18*s4*s8) + s43*(2*s14*s18*s8 + s19*s45*s8) + s39*(s19*s42*s8 + s16*(2*s22*s40 - 2*s23*s8)) + s31*(s49 + s13*s18*s9) + s30*(s50 + s18*s3*s9);
				buffer__[6*i+5] = s13*s19*s25 + (s13*s19*s26*s27)/2. + (s13*s19*s33*s34)/2. + s28*((s13*s19*s29)/2. + s35) + s36*(s13*s19*s38 + 2*s13*s18*s4) + s43*(s13*s19*s45 + s16*(-2*s13*s23 + 2*s22*s48)) + s24*(s14*s18*s3 + s50) + s13*s18*s3*s30*s8 + s31*(s51 + s14*s18*s8) + s39*(s13*s19*s42 + 2*s13*s18*s9);
			}
		}

        ptoc(ClassName()+"::DFarToHulls");
        
    }

	}; // SimplicialMeshDetails<1,3,Real,Int>

//----------------------------------------------------------------------------------------------

    template<typename Real, typename Int>
    struct SimplicialMeshDetails<1,4,Real,Int>
    {
	private:

		const Int thread_count = 1;

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(std::max(static_cast<Int>(1),thread_count_))
		{}
	
		inline Int ThreadCount() const
		{
			return thread_count;
		}
	
		inline std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<1,4,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void ComputeNearFarData( 
		const Tensor2<Real,Int> & V_coords,
		const Tensor2<Int ,Int> & simplices,
			  Tensor2<Real,Int> & P_near,
			  Tensor2<Real,Int> & P_far
	) const
    {
        ptic(ClassName()+"::ComputeNearFarData");
        
        //Int size       = 2;
        //Int amb_dim    = 4;
        //Int dom_dim    = 1;

        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			ptr<Real> V_coords__      = V_coords.data();	
			const Int  * restrict const simplices__     = simplices.data();

			Real hull    [2][4];
			Real df      [4][1];
			Real dfdagger[1][4];
			Real g       [1][1];
			Real ginv    [1][1];

			Int simplex  [2];
			
			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {
				Real * restrict const near = P_near.data(i);                    
				Real * restrict const far  = P_far.data(i);   
            
				simplex[0] = simplices__[2*i +0];
				simplex[1] = simplices__[2*i +1];

				near[1] = hull[0][0] = V_coords__[4*simplex[0]+0];
				near[2] = hull[0][1] = V_coords__[4*simplex[0]+1];
				near[3] = hull[0][2] = V_coords__[4*simplex[0]+2];
				near[4] = hull[0][3] = V_coords__[4*simplex[0]+3];
				near[5] = hull[1][0] = V_coords__[4*simplex[1]+0];
				near[6] = hull[1][1] = V_coords__[4*simplex[1]+1];
				near[7] = hull[1][2] = V_coords__[4*simplex[1]+2];
				near[8] = hull[1][3] = V_coords__[4*simplex[1]+3];

				far[1] = static_cast<Real>(0.5) * ( hull[0][0] + hull[1][0] );
				far[2] = static_cast<Real>(0.5) * ( hull[0][1] + hull[1][1] );
				far[3] = static_cast<Real>(0.5) * ( hull[0][2] + hull[1][2] );
				far[4] = static_cast<Real>(0.5) * ( hull[0][3] + hull[1][3] );

				df[0][0] = hull[1][0] - hull[0][0];
				df[1][0] = hull[1][1] - hull[0][1];
				df[2][0] = hull[1][2] - hull[0][2];
				df[3][0] = hull[1][3] - hull[0][3];

				g[0][0] = df[0][0] * df[0][0] + df[1][0] * df[1][0] + df[2][0] * df[2][0] + df[3][0] * df[3][0];


				near[0] = far[0] = sqrt( fabs(g[0][0]) );

                ginv[0][0] = static_cast<Real>(1)/g[0][0];
                
                //  dfdagger = g^{-1} * df^T (1 x 4 matrix)
				dfdagger[0][0] = ginv[0][0] * df[0][0];
				dfdagger[0][1] = ginv[0][0] * df[1][0];
				dfdagger[0][2] = ginv[0][0] * df[2][0];
				dfdagger[0][3] = ginv[0][0] * df[3][0];
            
				near[ 9] = far[ 5]  = static_cast<Real>(1) - df[0][0] * dfdagger[0][0];
				near[10] = far[ 6]  =    - df[0][0] * dfdagger[0][1];
				near[11] = far[ 7]  =    - df[0][0] * dfdagger[0][2];
				near[12] = far[ 8]  =    - df[0][0] * dfdagger[0][3];
				near[13] = far[ 9]  = static_cast<Real>(1) - df[1][0] * dfdagger[0][1];
				near[14] = far[10]  =    - df[1][0] * dfdagger[0][2];
				near[15] = far[11]  =    - df[1][0] * dfdagger[0][3];
				near[16] = far[12]  = static_cast<Real>(1) - df[2][0] * dfdagger[0][2];
				near[17] = far[13]  =    - df[2][0] * dfdagger[0][3];
				near[18] = far[14]  = static_cast<Real>(1) - df[3][0] * dfdagger[0][3];

            } // for( Int i = i_begin; i < i_end; ++i )

		} // #pragma omp parallel num_threads( ThreadCount() )

        ptoc(ClassName()+"::ComputeNearFarData");
    }

    void ComputeNearFarDataOps( 
		const Tensor2<Real,Int> & V_coords,
        const Tensor2<Int ,Int> & simplices,
		      Tensor2<Real,Int> & P_coords,
		      Tensor3<Real,Int> & P_hull_coords,
		      Tensor2<Real,Int> & P_near,
		      Tensor2<Real,Int> & P_far,
		Sparse::MatrixCSR<Real,Int,Int> & DiffOp,
		Sparse::MatrixCSR<Real,Int,Int> & AvOp 
	) const
    {
        ptic(ClassName()+"::ComputeNearFarDataOps");
        
        //Int size       = 2;
        //Int amb_dim    = 4;
        //Int dom_dim    = 1;
        
        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			mut<Int>  AvOp_outer = AvOp.Outer().data();
			mut<Int>  AvOp_inner = AvOp.Inner().data();
			mut<Real> AvOp_value = AvOp.Values().data();

			mut<Int>  DiffOp_outer = DiffOp.Outer().data();
			mut<Int>  DiffOp_inner = DiffOp.Inner().data();
			mut<Real> DiffOp_value = DiffOp.Value().data();

			ptr<Real> V_coords__      = V_coords.data();
			
			ptr<Int>  simplices__     = simplices.data();
		    mut<Real> P_hull_coords__ = P_hull_coords.data();
			mut<Real> P_coords__      = P_coords.data();

			Real df       [4][1];
			Real dfdagger [1][4];
			Real g        [1][1];
			Real ginv     [1][1];

			Int simplex        [2];
			Int sorted_simplex [2];

			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {

				mut<Real> near = P_near.data(i);                    
				mut<Real> far  = P_far.data(i);

				simplex[0] = sorted_simplex[0] = simplices__[2*i +0];
				simplex[1] = sorted_simplex[1] = simplices__[2*i +1];
                  
                // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
                std::sort( sorted_simplex, sorted_simplex + 2 );

				AvOp_outer[i+1] = (i+1) * 2;  
                      
				AvOp_inner[2*i+0] = sorted_simplex[0];
				AvOp_inner[2*i+1] = sorted_simplex[1];

				AvOp_value[2*i+0] = 0.5;
				AvOp_value[2*i+1] = 0.5;

				DiffOp_outer[4*i+0] = (4 * i + 0) * 2;
				DiffOp_outer[4*i+1] = (4 * i + 1) * 2;
				DiffOp_outer[4*i+2] = (4 * i + 2) * 2;
				DiffOp_outer[4*i+3] = (4 * i + 3) * 2;

				DiffOp_inner[(i*4+0)*2+0] = sorted_simplex[0];
				DiffOp_inner[(i*4+1)*2+0] = sorted_simplex[0];
				DiffOp_inner[(i*4+2)*2+0] = sorted_simplex[0];
				DiffOp_inner[(i*4+3)*2+0] = sorted_simplex[0];
				DiffOp_inner[(i*4+0)*2+1] = sorted_simplex[1];
				DiffOp_inner[(i*4+1)*2+1] = sorted_simplex[1];
				DiffOp_inner[(i*4+2)*2+1] = sorted_simplex[1];
				DiffOp_inner[(i*4+3)*2+1] = sorted_simplex[1];

				near[1] = P_hull_coords__[8*i+0] = V_coords__[4*simplex[0]+0];
				near[2] = P_hull_coords__[8*i+1] = V_coords__[4*simplex[0]+1];
				near[3] = P_hull_coords__[8*i+2] = V_coords__[4*simplex[0]+2];
				near[4] = P_hull_coords__[8*i+3] = V_coords__[4*simplex[0]+3];
				near[5] = P_hull_coords__[8*i+4] = V_coords__[4*simplex[1]+0];
				near[6] = P_hull_coords__[8*i+5] = V_coords__[4*simplex[1]+1];
				near[7] = P_hull_coords__[8*i+6] = V_coords__[4*simplex[1]+2];
				near[8] = P_hull_coords__[8*i+7] = V_coords__[4*simplex[1]+3];

				far[1] = P_coords__[4*i+0] = 0.5 * ( P_hull_coords__[8*i+0] + P_hull_coords__[8*i+4] );
				far[2] = P_coords__[4*i+1] = 0.5 * ( P_hull_coords__[8*i+1] + P_hull_coords__[8*i+5] );
				far[3] = P_coords__[4*i+2] = 0.5 * ( P_hull_coords__[8*i+2] + P_hull_coords__[8*i+6] );
				far[4] = P_coords__[4*i+3] = 0.5 * ( P_hull_coords__[8*i+3] + P_hull_coords__[8*i+7] );

				df[0][0] = V_coords__[4*sorted_simplex[1]+0] - V_coords__[4*sorted_simplex[0]+0];
				df[1][0] = V_coords__[4*sorted_simplex[1]+1] - V_coords__[4*sorted_simplex[0]+1];
				df[2][0] = V_coords__[4*sorted_simplex[1]+2] - V_coords__[4*sorted_simplex[0]+2];
				df[3][0] = V_coords__[4*sorted_simplex[1]+3] - V_coords__[4*sorted_simplex[0]+3];

				g[0][0] = df[0][0] * df[0][0] + df[1][0] * df[1][0] + df[2][0] * df[2][0] + df[3][0] * df[3][0];

                near[0] = far[0] = sqrt( fabs(g[0][0]) );

                ginv[0][0] =  static_cast<Real>(1)/g[0][0];
                
                //  dfdagger = g^{-1} * df^T (1 x 4 matrix)
				dfdagger[0][0] = ginv[0][0] * df[0][0];
				dfdagger[0][1] = ginv[0][0] * df[1][0];
				dfdagger[0][2] = ginv[0][0] * df[2][0];
				dfdagger[0][3] = ginv[0][0] * df[3][0];

				near[ 9] = far[ 5] = static_cast<Real>(1.0) - df[0][0] * dfdagger[0][0];
				near[10] = far[ 6] =    - df[0][0] * dfdagger[0][1];
				near[11] = far[ 7] =    - df[0][0] * dfdagger[0][2];
				near[12] = far[ 8] =    - df[0][0] * dfdagger[0][3];
				near[13] = far[ 9] = static_cast<Real>(1.0) - df[1][0] * dfdagger[0][1];
				near[14] = far[10] =    - df[1][0] * dfdagger[0][2];
				near[15] = far[11] =    - df[1][0] * dfdagger[0][3];
				near[16] = far[12] = static_cast<Real>(1.0) - df[2][0] * dfdagger[0][2];
				near[17] = far[13] =    - df[2][0] * dfdagger[0][3];
				near[18] = far[14] = static_cast<Real>(1.0) - df[3][0] * dfdagger[0][3];

                // derivative operator  (4 x 2 matrix)

                Real * Df = &DiffOp_value[ 8 * i ];

				Df[ 0] = - dfdagger[0][0];
				Df[ 1] =   dfdagger[0][0];
				Df[ 2] = - dfdagger[0][1];
				Df[ 3] =   dfdagger[0][1];
				Df[ 4] = - dfdagger[0][2];
				Df[ 5] =   dfdagger[0][2];
				Df[ 6] = - dfdagger[0][3];
				Df[ 7] =   dfdagger[0][3];

            }
        }
        ptoc(ClassName()+"::ComputeNearFarDataOps");
    }

	void DNearToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_near, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DNearToHulls");

        if( P_D_near.Dimension(1) != 19 )
        {
            eprint("in DNearToHulls: P_D_near.Dimension(1) != 19. Aborting");
        }

		ptr<Real> V_coords__  = V_coords.data();
		ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_near__  = P_D_near.data();
		mut<Real> buffer__    = buffer.data();
        
        if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[4*simplices__[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = V_coords__[4*simplices__[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = V_coords__[4*simplices__[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = V_coords__[4*simplices__[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = V_coords__[4*simplices__[2*i+0]+2];
				const Real s11 = -s10;
				const Real s12 = V_coords__[4*simplices__[2*i+1]+2];
				const Real s13 = s11 + s12;
				const Real s14 = s13*s13;
				const Real s15 = V_coords__[4*simplices__[2*i+0]+3];
				const Real s16 = -s15;
				const Real s17 = V_coords__[4*simplices__[2*i+1]+3];
				const Real s18 = s16 + s17;
				const Real s19 = s18*s18;
				const Real s20 = s14 + s19 + s4 + s9;
				const Real s21 = sqrt(s20);
				const Real s22 = s20*s21;
				const Real s23 = 1/s22;
				const Real s24 = 1/s21;
				const Real s25 = s3*s4;
				const Real s26 = s20*s20;
				const Real s27 = 1/s26;
				const Real s28 = 1/s20;
				const Real s29 = P_D_near__[19*i+11];
				const Real s30 = P_D_near__[19*i+12];
				const Real s31 = P_D_near__[19*i+17];
				const Real s32 = P_D_near__[19*i+0];
				const Real s33 = P_D_near__[19*i+1];
				const Real s34 = P_D_near__[19*i+3];
				const Real s35 = P_D_near__[19*i+4];
				const Real s36 = P_D_near__[19*i+5];
				const Real s37 = P_D_near__[19*i+6];
				const Real s38 = P_D_near__[19*i+7];
				const Real s39 = P_D_near__[19*i+8];
				const Real s40 = P_D_near__[19*i+10];
				const Real s41 = P_D_near__[19*i+14];
				const Real s42 = s13*s24;
				const Real s43 = P_D_near__[19*i+15];
				const Real s44 = s18*s24;
				const Real s45 = P_D_near__[19*i+2];
				const Real s46 = P_D_near__[19*i+9];
				const Real s47 = -(s28*s4);
				const Real s48 = 1 + s47;
				const Real s49 = P_D_near__[19*i+13];
				const Real s50 = s8*s9;
				const Real s51 = -(s28*s9);
				const Real s52 = 1 + s51;
				const Real s53 = P_D_near__[19*i+16];
				const Real s54 = -(s14*s28);
				const Real s55 = 1 + s54;
				const Real s56 = P_D_near__[19*i+18];
				const Real s57 = -(s19*s28);
				const Real s58 = 1 + s57;
				const Real s59 = s24*s3;
				const Real s60 = s24*s8;
				const Real s61 = s13*s14;
				const Real s62 = s18*s19;
				const Real s63 = -(s13*s24);
				const Real s64 = -(s18*s24);
				const Real s65 = -(s24*s3);
				const Real s66 = -(s24*s8);
				buffer__[8*i+0] += -(s13*s18*s23*s3*s31) - s24*s3*s32 + (s21 - s0*s24*s3)*s33 - s10*s24*s3*s34 - s15*s24*s3*s35 - s2*s24*s3*s36 - s12*s24*s3*s38 - s17*s24*s3*s39 + s29*(-(s13*s23*s4) + s42) + s30*(-(s18*s23*s4) + s44) + s46*(s21*(-2*s25*s27 + 2*s28*s3) - s24*s3*s48) - s24*s3*s45*s5 + s53*(-2*s14*s23*s3 - s24*s3*s55) + s56*(-2*s19*s23*s3 - s24*s3*s58) - s24*s3*s37*s7 - s13*s23*s3*s41*s8 - s18*s23*s3*s43*s8 + s40*(s60 - s23*s4*s8) + s49*(-(s24*s3*s52) - 2*s23*s3*s9);
				buffer__[8*i+1] += -(s13*s23*s29*s3*s8) - s18*s23*s3*s30*s8 - s13*s18*s23*s31*s8 - s24*s32*s8 - s0*s24*s33*s8 - s10*s24*s34*s8 - s15*s24*s35*s8 - s2*s24*s36*s8 - s12*s24*s38*s8 - s17*s24*s39*s8 - s24*s37*s7*s8 + s46*(-2*s23*s4*s8 - s24*s48*s8) + s45*(s21 - s24*s5*s8) + s53*(-2*s14*s23*s8 - s24*s55*s8) + s56*(-2*s19*s23*s8 - s24*s58*s8) + s49*(-(s24*s52*s8) + s21*(-2*s27*s50 + 2*s28*s8)) + s41*(s42 - s13*s23*s9) + s43*(s44 - s18*s23*s9) + s40*(s59 - s23*s3*s9);
				buffer__[8*i+2] += -(s13*s18*s23*s3*s30) - s13*s24*s32 - s0*s13*s24*s33 + (s21 - s10*s13*s24)*s34 - s13*s15*s24*s35 - s13*s2*s24*s36 - s12*s13*s24*s38 - s13*s17*s24*s39 + s31*(-(s14*s18*s23) + s44) + s46*(-2*s13*s23*s4 - s13*s24*s48) - s13*s24*s45*s5 + s56*(-2*s13*s19*s23 - s13*s24*s58) + s29*(-(s14*s23*s3) + s59) + s53*(-(s13*s24*s55) + s21*(2*s13*s28 - 2*s27*s61)) - s13*s24*s37*s7 - s13*s23*s3*s40*s8 - s13*s18*s23*s43*s8 + s41*(s60 - s14*s23*s8) + s49*(-(s13*s24*s52) - 2*s13*s23*s9);
				buffer__[8*i+3] += -(s13*s18*s23*s29*s3) - s18*s24*s32 - s0*s18*s24*s33 - s10*s18*s24*s34 + (s21 - s15*s18*s24)*s35 - s18*s2*s24*s36 - s12*s18*s24*s38 - s17*s18*s24*s39 + s31*(-(s13*s19*s23) + s42) + s46*(-2*s18*s23*s4 - s18*s24*s48) - s18*s24*s45*s5 + s53*(-2*s14*s18*s23 - s18*s24*s55) + s30*(-(s19*s23*s3) + s59) + s56*(-(s18*s24*s58) + s21*(2*s18*s28 - 2*s27*s62)) - s18*s24*s37*s7 - s18*s23*s3*s40*s8 - s13*s18*s23*s41*s8 + s43*(s60 - s19*s23*s8) + s49*(-(s18*s24*s52) - 2*s18*s23*s9);
				buffer__[8*i+4] += s13*s18*s23*s3*s31 + s24*s3*s32 + s0*s24*s3*s33 + s10*s24*s3*s34 + s15*s24*s3*s35 + (s21 + s2*s24*s3)*s36 + s12*s24*s3*s38 + s17*s24*s3*s39 + s46*(s21*(2*s25*s27 - 2*s28*s3) + s24*s3*s48) + s24*s3*s45*s5 + s53*(2*s14*s23*s3 + s24*s3*s55) + s56*(2*s19*s23*s3 + s24*s3*s58) + s29*(s13*s23*s4 + s63) + s30*(s18*s23*s4 + s64) + s24*s3*s37*s7 + s13*s23*s3*s41*s8 + s18*s23*s3*s43*s8 + s40*(s66 + s23*s4*s8) + s49*(s24*s3*s52 + 2*s23*s3*s9);
				buffer__[8*i+5] += s13*s23*s29*s3*s8 + s18*s23*s3*s30*s8 + s13*s18*s23*s31*s8 + s24*s32*s8 + s0*s24*s33*s8 + s10*s24*s34*s8 + s15*s24*s35*s8 + s2*s24*s36*s8 + s12*s24*s38*s8 + s17*s24*s39*s8 + s24*s45*s5*s8 + s46*(2*s23*s4*s8 + s24*s48*s8) + s53*(2*s14*s23*s8 + s24*s55*s8) + s56*(2*s19*s23*s8 + s24*s58*s8) + s37*(s21 + s24*s7*s8) + s49*(s24*s52*s8 + s21*(2*s27*s50 - 2*s28*s8)) + s41*(s63 + s13*s23*s9) + s43*(s64 + s18*s23*s9) + s40*(s65 + s23*s3*s9);
				buffer__[8*i+6] += s13*s18*s23*s3*s30 + s13*s24*s32 + s0*s13*s24*s33 + s10*s13*s24*s34 + s13*s15*s24*s35 + s13*s2*s24*s36 + (s21 + s12*s13*s24)*s38 + s13*s17*s24*s39 + s46*(2*s13*s23*s4 + s13*s24*s48) + s13*s24*s45*s5 + s56*(2*s13*s19*s23 + s13*s24*s58) + s53*(s13*s24*s55 + s21*(-2*s13*s28 + 2*s27*s61)) + s31*(s14*s18*s23 + s64) + s29*(s14*s23*s3 + s65) + s13*s24*s37*s7 + s13*s23*s3*s40*s8 + s13*s18*s23*s43*s8 + s41*(s66 + s14*s23*s8) + s49*(s13*s24*s52 + 2*s13*s23*s9);
				buffer__[8*i+7] += s13*s18*s23*s29*s3 + s18*s24*s32 + s0*s18*s24*s33 + s10*s18*s24*s34 + s15*s18*s24*s35 + s18*s2*s24*s36 + s12*s18*s24*s38 + (s21 + s17*s18*s24)*s39 + s46*(2*s18*s23*s4 + s18*s24*s48) + s18*s24*s45*s5 + s53*(2*s14*s18*s23 + s18*s24*s55) + s56*(s18*s24*s58 + s21*(-2*s18*s28 + 2*s27*s62)) + s31*(s13*s19*s23 + s63) + s30*(s19*s23*s3 + s65) + s18*s24*s37*s7 + s18*s23*s3*s40*s8 + s13*s18*s23*s41*s8 + s43*(s66 + s19*s23*s8) + s49*(s18*s24*s52 + 2*s18*s23*s9);
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[4*simplices__[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = V_coords__[4*simplices__[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = V_coords__[4*simplices__[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = V_coords__[4*simplices__[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = V_coords__[4*simplices__[2*i+0]+2];
				const Real s11 = -s10;
				const Real s12 = V_coords__[4*simplices__[2*i+1]+2];
				const Real s13 = s11 + s12;
				const Real s14 = s13*s13;
				const Real s15 = V_coords__[4*simplices__[2*i+0]+3];
				const Real s16 = -s15;
				const Real s17 = V_coords__[4*simplices__[2*i+1]+3];
				const Real s18 = s16 + s17;
				const Real s19 = s18*s18;
				const Real s20 = s14 + s19 + s4 + s9;
				const Real s21 = sqrt(s20);
				const Real s22 = s20*s21;
				const Real s23 = 1/s22;
				const Real s24 = 1/s21;
				const Real s25 = s3*s4;
				const Real s26 = s20*s20;
				const Real s27 = 1/s26;
				const Real s28 = 1/s20;
				const Real s29 = P_D_near__[19*i+11];
				const Real s30 = P_D_near__[19*i+12];
				const Real s31 = P_D_near__[19*i+17];
				const Real s32 = P_D_near__[19*i+0];
				const Real s33 = P_D_near__[19*i+1];
				const Real s34 = P_D_near__[19*i+3];
				const Real s35 = P_D_near__[19*i+4];
				const Real s36 = P_D_near__[19*i+5];
				const Real s37 = P_D_near__[19*i+6];
				const Real s38 = P_D_near__[19*i+7];
				const Real s39 = P_D_near__[19*i+8];
				const Real s40 = P_D_near__[19*i+10];
				const Real s41 = P_D_near__[19*i+14];
				const Real s42 = s13*s24;
				const Real s43 = P_D_near__[19*i+15];
				const Real s44 = s18*s24;
				const Real s45 = P_D_near__[19*i+2];
				const Real s46 = P_D_near__[19*i+9];
				const Real s47 = -(s28*s4);
				const Real s48 = 1 + s47;
				const Real s49 = P_D_near__[19*i+13];
				const Real s50 = s8*s9;
				const Real s51 = -(s28*s9);
				const Real s52 = 1 + s51;
				const Real s53 = P_D_near__[19*i+16];
				const Real s54 = -(s14*s28);
				const Real s55 = 1 + s54;
				const Real s56 = P_D_near__[19*i+18];
				const Real s57 = -(s19*s28);
				const Real s58 = 1 + s57;
				const Real s59 = s24*s3;
				const Real s60 = s24*s8;
				const Real s61 = s13*s14;
				const Real s62 = s18*s19;
				const Real s63 = -(s13*s24);
				const Real s64 = -(s18*s24);
				const Real s65 = -(s24*s3);
				const Real s66 = -(s24*s8);
				buffer__[8*i+0] = -(s13*s18*s23*s3*s31) - s24*s3*s32 + (s21 - s0*s24*s3)*s33 - s10*s24*s3*s34 - s15*s24*s3*s35 - s2*s24*s3*s36 - s12*s24*s3*s38 - s17*s24*s3*s39 + s29*(-(s13*s23*s4) + s42) + s30*(-(s18*s23*s4) + s44) + s46*(s21*(-2*s25*s27 + 2*s28*s3) - s24*s3*s48) - s24*s3*s45*s5 + s53*(-2*s14*s23*s3 - s24*s3*s55) + s56*(-2*s19*s23*s3 - s24*s3*s58) - s24*s3*s37*s7 - s13*s23*s3*s41*s8 - s18*s23*s3*s43*s8 + s40*(s60 - s23*s4*s8) + s49*(-(s24*s3*s52) - 2*s23*s3*s9);
				buffer__[8*i+1] = -(s13*s23*s29*s3*s8) - s18*s23*s3*s30*s8 - s13*s18*s23*s31*s8 - s24*s32*s8 - s0*s24*s33*s8 - s10*s24*s34*s8 - s15*s24*s35*s8 - s2*s24*s36*s8 - s12*s24*s38*s8 - s17*s24*s39*s8 - s24*s37*s7*s8 + s46*(-2*s23*s4*s8 - s24*s48*s8) + s45*(s21 - s24*s5*s8) + s53*(-2*s14*s23*s8 - s24*s55*s8) + s56*(-2*s19*s23*s8 - s24*s58*s8) + s49*(-(s24*s52*s8) + s21*(-2*s27*s50 + 2*s28*s8)) + s41*(s42 - s13*s23*s9) + s43*(s44 - s18*s23*s9) + s40*(s59 - s23*s3*s9);
				buffer__[8*i+2] = -(s13*s18*s23*s3*s30) - s13*s24*s32 - s0*s13*s24*s33 + (s21 - s10*s13*s24)*s34 - s13*s15*s24*s35 - s13*s2*s24*s36 - s12*s13*s24*s38 - s13*s17*s24*s39 + s31*(-(s14*s18*s23) + s44) + s46*(-2*s13*s23*s4 - s13*s24*s48) - s13*s24*s45*s5 + s56*(-2*s13*s19*s23 - s13*s24*s58) + s29*(-(s14*s23*s3) + s59) + s53*(-(s13*s24*s55) + s21*(2*s13*s28 - 2*s27*s61)) - s13*s24*s37*s7 - s13*s23*s3*s40*s8 - s13*s18*s23*s43*s8 + s41*(s60 - s14*s23*s8) + s49*(-(s13*s24*s52) - 2*s13*s23*s9);
				buffer__[8*i+3] = -(s13*s18*s23*s29*s3) - s18*s24*s32 - s0*s18*s24*s33 - s10*s18*s24*s34 + (s21 - s15*s18*s24)*s35 - s18*s2*s24*s36 - s12*s18*s24*s38 - s17*s18*s24*s39 + s31*(-(s13*s19*s23) + s42) + s46*(-2*s18*s23*s4 - s18*s24*s48) - s18*s24*s45*s5 + s53*(-2*s14*s18*s23 - s18*s24*s55) + s30*(-(s19*s23*s3) + s59) + s56*(-(s18*s24*s58) + s21*(2*s18*s28 - 2*s27*s62)) - s18*s24*s37*s7 - s18*s23*s3*s40*s8 - s13*s18*s23*s41*s8 + s43*(s60 - s19*s23*s8) + s49*(-(s18*s24*s52) - 2*s18*s23*s9);
				buffer__[8*i+4] = s13*s18*s23*s3*s31 + s24*s3*s32 + s0*s24*s3*s33 + s10*s24*s3*s34 + s15*s24*s3*s35 + (s21 + s2*s24*s3)*s36 + s12*s24*s3*s38 + s17*s24*s3*s39 + s46*(s21*(2*s25*s27 - 2*s28*s3) + s24*s3*s48) + s24*s3*s45*s5 + s53*(2*s14*s23*s3 + s24*s3*s55) + s56*(2*s19*s23*s3 + s24*s3*s58) + s29*(s13*s23*s4 + s63) + s30*(s18*s23*s4 + s64) + s24*s3*s37*s7 + s13*s23*s3*s41*s8 + s18*s23*s3*s43*s8 + s40*(s66 + s23*s4*s8) + s49*(s24*s3*s52 + 2*s23*s3*s9);
				buffer__[8*i+5] = s13*s23*s29*s3*s8 + s18*s23*s3*s30*s8 + s13*s18*s23*s31*s8 + s24*s32*s8 + s0*s24*s33*s8 + s10*s24*s34*s8 + s15*s24*s35*s8 + s2*s24*s36*s8 + s12*s24*s38*s8 + s17*s24*s39*s8 + s24*s45*s5*s8 + s46*(2*s23*s4*s8 + s24*s48*s8) + s53*(2*s14*s23*s8 + s24*s55*s8) + s56*(2*s19*s23*s8 + s24*s58*s8) + s37*(s21 + s24*s7*s8) + s49*(s24*s52*s8 + s21*(2*s27*s50 - 2*s28*s8)) + s41*(s63 + s13*s23*s9) + s43*(s64 + s18*s23*s9) + s40*(s65 + s23*s3*s9);
				buffer__[8*i+6] = s13*s18*s23*s3*s30 + s13*s24*s32 + s0*s13*s24*s33 + s10*s13*s24*s34 + s13*s15*s24*s35 + s13*s2*s24*s36 + (s21 + s12*s13*s24)*s38 + s13*s17*s24*s39 + s46*(2*s13*s23*s4 + s13*s24*s48) + s13*s24*s45*s5 + s56*(2*s13*s19*s23 + s13*s24*s58) + s53*(s13*s24*s55 + s21*(-2*s13*s28 + 2*s27*s61)) + s31*(s14*s18*s23 + s64) + s29*(s14*s23*s3 + s65) + s13*s24*s37*s7 + s13*s23*s3*s40*s8 + s13*s18*s23*s43*s8 + s41*(s66 + s14*s23*s8) + s49*(s13*s24*s52 + 2*s13*s23*s9);
				buffer__[8*i+7] = s13*s18*s23*s29*s3 + s18*s24*s32 + s0*s18*s24*s33 + s10*s18*s24*s34 + s15*s18*s24*s35 + s18*s2*s24*s36 + s12*s18*s24*s38 + (s21 + s17*s18*s24)*s39 + s46*(2*s18*s23*s4 + s18*s24*s48) + s18*s24*s45*s5 + s53*(2*s14*s18*s23 + s18*s24*s55) + s56*(s18*s24*s58 + s21*(-2*s18*s28 + 2*s27*s62)) + s31*(s13*s19*s23 + s63) + s30*(s19*s23*s3 + s65) + s18*s24*s37*s7 + s18*s23*s3*s40*s8 + s13*s18*s23*s41*s8 + s43*(s66 + s19*s23*s8) + s49*(s18*s24*s52 + 2*s18*s23*s9);
			}
		}

        ptoc(ClassName()+"::DNearToHulls");
        
    }

	void DFarToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_far, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DFarToHulls");

        if( P_D_far.Dimension(1) != 15 )
        {
            eprint("in DFarToHulls: P_D_far.Dimension(1) != 15. Aborting");
        }

		ptr<Real> V_coords__  = V_coords.data();
		ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_far__   = P_D_far.data();
		mut<Real> buffer__    = buffer.data();
        
        if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[4*simplices__[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = V_coords__[4*simplices__[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = V_coords__[4*simplices__[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = V_coords__[4*simplices__[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = V_coords__[4*simplices__[2*i+0]+2];
				const Real s11 = -s10;
				const Real s12 = V_coords__[4*simplices__[2*i+1]+2];
				const Real s13 = s11 + s12;
				const Real s14 = s13*s13;
				const Real s15 = V_coords__[4*simplices__[2*i+0]+3];
				const Real s16 = -s15;
				const Real s17 = V_coords__[4*simplices__[2*i+1]+3];
				const Real s18 = s16 + s17;
				const Real s19 = s18*s18;
				const Real s20 = s14 + s19 + s4 + s9;
				const Real s21 = sqrt(s20);
				const Real s22 = s20*s21;
				const Real s23 = 1/s22;
				const Real s24 = 1/s21;
				const Real s25 = s3*s4;
				const Real s26 = s20*s20;
				const Real s27 = 1/s26;
				const Real s28 = 1/s20;
				const Real s29 = P_D_far__[15*i+7];
				const Real s30 = P_D_far__[15*i+8];
				const Real s31 = P_D_far__[15*i+13];
				const Real s32 = P_D_far__[15*i+0];
				const Real s33 = P_D_far__[15*i+1];
				const Real s34 = s0 + s2;
				const Real s35 = P_D_far__[15*i+3];
				const Real s36 = s10 + s12;
				const Real s37 = P_D_far__[15*i+4];
				const Real s38 = s15 + s17;
				const Real s39 = P_D_far__[15*i+6];
				const Real s40 = P_D_far__[15*i+10];
				const Real s41 = s13*s24;
				const Real s42 = P_D_far__[15*i+11];
				const Real s43 = s18*s24;
				const Real s44 = P_D_far__[15*i+2];
				const Real s45 = s5 + s7;
				const Real s46 = s21/2.;
				const Real s47 = P_D_far__[15*i+5];
				const Real s48 = -(s28*s4);
				const Real s49 = 1 + s48;
				const Real s50 = P_D_far__[15*i+9];
				const Real s51 = s8*s9;
				const Real s52 = -(s28*s9);
				const Real s53 = 1 + s52;
				const Real s54 = P_D_far__[15*i+12];
				const Real s55 = -(s14*s28);
				const Real s56 = 1 + s55;
				const Real s57 = P_D_far__[15*i+14];
				const Real s58 = -(s19*s28);
				const Real s59 = 1 + s58;
				const Real s60 = s24*s3;
				const Real s61 = s24*s8;
				const Real s62 = s13*s14;
				const Real s63 = s18*s19;
				const Real s64 = -(s13*s24);
				const Real s65 = -(s18*s24);
				const Real s66 = -(s24*s3);
				const Real s67 = -(s24*s8);
				buffer__[8*i+0] += -(s13*s18*s23*s3*s31) - s24*s3*s32 - (s24*s3*s35*s36)/2. - (s24*s3*s37*s38)/2. + s29*(-(s13*s23*s4) + s41) + s30*(-(s18*s23*s4) + s43) - (s24*s3*s44*s45)/2. + s33*(-0.5*(s24*s3*s34) + s46) + s47*(s21*(-2*s25*s27 + 2*s28*s3) - s24*s3*s49) + s54*(-2*s14*s23*s3 - s24*s3*s56) + s57*(-2*s19*s23*s3 - s24*s3*s59) - s13*s23*s3*s40*s8 - s18*s23*s3*s42*s8 + s39*(s61 - s23*s4*s8) + s50*(-(s24*s3*s53) - 2*s23*s3*s9);
				buffer__[8*i+1] += -(s13*s23*s29*s3*s8) - s18*s23*s3*s30*s8 - s13*s18*s23*s31*s8 - s24*s32*s8 - (s24*s33*s34*s8)/2. - (s24*s35*s36*s8)/2. - (s24*s37*s38*s8)/2. + s44*(s46 - (s24*s45*s8)/2.) + s47*(-2*s23*s4*s8 - s24*s49*s8) + s54*(-2*s14*s23*s8 - s24*s56*s8) + s57*(-2*s19*s23*s8 - s24*s59*s8) + s50*(-(s24*s53*s8) + s21*(-2*s27*s51 + 2*s28*s8)) + s40*(s41 - s13*s23*s9) + s42*(s43 - s18*s23*s9) + s39*(s60 - s23*s3*s9);
				buffer__[8*i+2] += -(s13*s18*s23*s3*s30) - s13*s24*s32 - (s13*s24*s33*s34)/2. - (s13*s24*s37*s38)/2. + s31*(-(s14*s18*s23) + s43) - (s13*s24*s44*s45)/2. + s35*(-0.5*(s13*s24*s36) + s46) + s47*(-2*s13*s23*s4 - s13*s24*s49) + s57*(-2*s13*s19*s23 - s13*s24*s59) + s29*(-(s14*s23*s3) + s60) + s54*(-(s13*s24*s56) + s21*(2*s13*s28 - 2*s27*s62)) - s13*s23*s3*s39*s8 - s13*s18*s23*s42*s8 + s40*(s61 - s14*s23*s8) + s50*(-(s13*s24*s53) - 2*s13*s23*s9);
				buffer__[8*i+3] += -(s13*s18*s23*s29*s3) - s18*s24*s32 - (s18*s24*s33*s34)/2. - (s18*s24*s35*s36)/2. + s31*(-(s13*s19*s23) + s41) - (s18*s24*s44*s45)/2. + s37*(-0.5*(s18*s24*s38) + s46) + s47*(-2*s18*s23*s4 - s18*s24*s49) + s54*(-2*s14*s18*s23 - s18*s24*s56) + s30*(-(s19*s23*s3) + s60) + s57*(-(s18*s24*s59) + s21*(2*s18*s28 - 2*s27*s63)) - s18*s23*s3*s39*s8 - s13*s18*s23*s40*s8 + s42*(s61 - s19*s23*s8) + s50*(-(s18*s24*s53) - 2*s18*s23*s9);
				buffer__[8*i+4] += s13*s18*s23*s3*s31 + s24*s3*s32 + (s24*s3*s35*s36)/2. + (s24*s3*s37*s38)/2. + (s24*s3*s44*s45)/2. + s33*((s24*s3*s34)/2. + s46) + s47*(s21*(2*s25*s27 - 2*s28*s3) + s24*s3*s49) + s54*(2*s14*s23*s3 + s24*s3*s56) + s57*(2*s19*s23*s3 + s24*s3*s59) + s29*(s13*s23*s4 + s64) + s30*(s18*s23*s4 + s65) + s13*s23*s3*s40*s8 + s18*s23*s3*s42*s8 + s39*(s67 + s23*s4*s8) + s50*(s24*s3*s53 + 2*s23*s3*s9);
				buffer__[8*i+5] += s13*s23*s29*s3*s8 + s18*s23*s3*s30*s8 + s13*s18*s23*s31*s8 + s24*s32*s8 + (s24*s33*s34*s8)/2. + (s24*s35*s36*s8)/2. + (s24*s37*s38*s8)/2. + s44*(s46 + (s24*s45*s8)/2.) + s47*(2*s23*s4*s8 + s24*s49*s8) + s54*(2*s14*s23*s8 + s24*s56*s8) + s57*(2*s19*s23*s8 + s24*s59*s8) + s50*(s24*s53*s8 + s21*(2*s27*s51 - 2*s28*s8)) + s40*(s64 + s13*s23*s9) + s42*(s65 + s18*s23*s9) + s39*(s66 + s23*s3*s9);
				buffer__[8*i+6] += s13*s18*s23*s3*s30 + s13*s24*s32 + (s13*s24*s33*s34)/2. + (s13*s24*s37*s38)/2. + (s13*s24*s44*s45)/2. + s35*((s13*s24*s36)/2. + s46) + s47*(2*s13*s23*s4 + s13*s24*s49) + s57*(2*s13*s19*s23 + s13*s24*s59) + s54*(s13*s24*s56 + s21*(-2*s13*s28 + 2*s27*s62)) + s31*(s14*s18*s23 + s65) + s29*(s14*s23*s3 + s66) + s13*s23*s3*s39*s8 + s13*s18*s23*s42*s8 + s40*(s67 + s14*s23*s8) + s50*(s13*s24*s53 + 2*s13*s23*s9);
				buffer__[8*i+7] += s13*s18*s23*s29*s3 + s18*s24*s32 + (s18*s24*s33*s34)/2. + (s18*s24*s35*s36)/2. + (s18*s24*s44*s45)/2. + s37*((s18*s24*s38)/2. + s46) + s47*(2*s18*s23*s4 + s18*s24*s49) + s54*(2*s14*s18*s23 + s18*s24*s56) + s57*(s18*s24*s59 + s21*(-2*s18*s28 + 2*s27*s63)) + s31*(s13*s19*s23 + s64) + s30*(s19*s23*s3 + s66) + s18*s23*s3*s39*s8 + s13*s18*s23*s40*s8 + s42*(s67 + s19*s23*s8) + s50*(s18*s24*s53 + 2*s18*s23*s9);
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[4*simplices__[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = V_coords__[4*simplices__[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = V_coords__[4*simplices__[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = V_coords__[4*simplices__[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = V_coords__[4*simplices__[2*i+0]+2];
				const Real s11 = -s10;
				const Real s12 = V_coords__[4*simplices__[2*i+1]+2];
				const Real s13 = s11 + s12;
				const Real s14 = s13*s13;
				const Real s15 = V_coords__[4*simplices__[2*i+0]+3];
				const Real s16 = -s15;
				const Real s17 = V_coords__[4*simplices__[2*i+1]+3];
				const Real s18 = s16 + s17;
				const Real s19 = s18*s18;
				const Real s20 = s14 + s19 + s4 + s9;
				const Real s21 = sqrt(s20);
				const Real s22 = s20*s21;
				const Real s23 = 1/s22;
				const Real s24 = 1/s21;
				const Real s25 = s3*s4;
				const Real s26 = s20*s20;
				const Real s27 = 1/s26;
				const Real s28 = 1/s20;
				const Real s29 = P_D_far__[15*i+7];
				const Real s30 = P_D_far__[15*i+8];
				const Real s31 = P_D_far__[15*i+13];
				const Real s32 = P_D_far__[15*i+0];
				const Real s33 = P_D_far__[15*i+1];
				const Real s34 = s0 + s2;
				const Real s35 = P_D_far__[15*i+3];
				const Real s36 = s10 + s12;
				const Real s37 = P_D_far__[15*i+4];
				const Real s38 = s15 + s17;
				const Real s39 = P_D_far__[15*i+6];
				const Real s40 = P_D_far__[15*i+10];
				const Real s41 = s13*s24;
				const Real s42 = P_D_far__[15*i+11];
				const Real s43 = s18*s24;
				const Real s44 = P_D_far__[15*i+2];
				const Real s45 = s5 + s7;
				const Real s46 = s21/2.;
				const Real s47 = P_D_far__[15*i+5];
				const Real s48 = -(s28*s4);
				const Real s49 = 1 + s48;
				const Real s50 = P_D_far__[15*i+9];
				const Real s51 = s8*s9;
				const Real s52 = -(s28*s9);
				const Real s53 = 1 + s52;
				const Real s54 = P_D_far__[15*i+12];
				const Real s55 = -(s14*s28);
				const Real s56 = 1 + s55;
				const Real s57 = P_D_far__[15*i+14];
				const Real s58 = -(s19*s28);
				const Real s59 = 1 + s58;
				const Real s60 = s24*s3;
				const Real s61 = s24*s8;
				const Real s62 = s13*s14;
				const Real s63 = s18*s19;
				const Real s64 = -(s13*s24);
				const Real s65 = -(s18*s24);
				const Real s66 = -(s24*s3);
				const Real s67 = -(s24*s8);
				buffer__[8*i+0] = -(s13*s18*s23*s3*s31) - s24*s3*s32 - (s24*s3*s35*s36)/2. - (s24*s3*s37*s38)/2. + s29*(-(s13*s23*s4) + s41) + s30*(-(s18*s23*s4) + s43) - (s24*s3*s44*s45)/2. + s33*(-0.5*(s24*s3*s34) + s46) + s47*(s21*(-2*s25*s27 + 2*s28*s3) - s24*s3*s49) + s54*(-2*s14*s23*s3 - s24*s3*s56) + s57*(-2*s19*s23*s3 - s24*s3*s59) - s13*s23*s3*s40*s8 - s18*s23*s3*s42*s8 + s39*(s61 - s23*s4*s8) + s50*(-(s24*s3*s53) - 2*s23*s3*s9);
				buffer__[8*i+1] = -(s13*s23*s29*s3*s8) - s18*s23*s3*s30*s8 - s13*s18*s23*s31*s8 - s24*s32*s8 - (s24*s33*s34*s8)/2. - (s24*s35*s36*s8)/2. - (s24*s37*s38*s8)/2. + s44*(s46 - (s24*s45*s8)/2.) + s47*(-2*s23*s4*s8 - s24*s49*s8) + s54*(-2*s14*s23*s8 - s24*s56*s8) + s57*(-2*s19*s23*s8 - s24*s59*s8) + s50*(-(s24*s53*s8) + s21*(-2*s27*s51 + 2*s28*s8)) + s40*(s41 - s13*s23*s9) + s42*(s43 - s18*s23*s9) + s39*(s60 - s23*s3*s9);
				buffer__[8*i+2] = -(s13*s18*s23*s3*s30) - s13*s24*s32 - (s13*s24*s33*s34)/2. - (s13*s24*s37*s38)/2. + s31*(-(s14*s18*s23) + s43) - (s13*s24*s44*s45)/2. + s35*(-0.5*(s13*s24*s36) + s46) + s47*(-2*s13*s23*s4 - s13*s24*s49) + s57*(-2*s13*s19*s23 - s13*s24*s59) + s29*(-(s14*s23*s3) + s60) + s54*(-(s13*s24*s56) + s21*(2*s13*s28 - 2*s27*s62)) - s13*s23*s3*s39*s8 - s13*s18*s23*s42*s8 + s40*(s61 - s14*s23*s8) + s50*(-(s13*s24*s53) - 2*s13*s23*s9);
				buffer__[8*i+3] = -(s13*s18*s23*s29*s3) - s18*s24*s32 - (s18*s24*s33*s34)/2. - (s18*s24*s35*s36)/2. + s31*(-(s13*s19*s23) + s41) - (s18*s24*s44*s45)/2. + s37*(-0.5*(s18*s24*s38) + s46) + s47*(-2*s18*s23*s4 - s18*s24*s49) + s54*(-2*s14*s18*s23 - s18*s24*s56) + s30*(-(s19*s23*s3) + s60) + s57*(-(s18*s24*s59) + s21*(2*s18*s28 - 2*s27*s63)) - s18*s23*s3*s39*s8 - s13*s18*s23*s40*s8 + s42*(s61 - s19*s23*s8) + s50*(-(s18*s24*s53) - 2*s18*s23*s9);
				buffer__[8*i+4] = s13*s18*s23*s3*s31 + s24*s3*s32 + (s24*s3*s35*s36)/2. + (s24*s3*s37*s38)/2. + (s24*s3*s44*s45)/2. + s33*((s24*s3*s34)/2. + s46) + s47*(s21*(2*s25*s27 - 2*s28*s3) + s24*s3*s49) + s54*(2*s14*s23*s3 + s24*s3*s56) + s57*(2*s19*s23*s3 + s24*s3*s59) + s29*(s13*s23*s4 + s64) + s30*(s18*s23*s4 + s65) + s13*s23*s3*s40*s8 + s18*s23*s3*s42*s8 + s39*(s67 + s23*s4*s8) + s50*(s24*s3*s53 + 2*s23*s3*s9);
				buffer__[8*i+5] = s13*s23*s29*s3*s8 + s18*s23*s3*s30*s8 + s13*s18*s23*s31*s8 + s24*s32*s8 + (s24*s33*s34*s8)/2. + (s24*s35*s36*s8)/2. + (s24*s37*s38*s8)/2. + s44*(s46 + (s24*s45*s8)/2.) + s47*(2*s23*s4*s8 + s24*s49*s8) + s54*(2*s14*s23*s8 + s24*s56*s8) + s57*(2*s19*s23*s8 + s24*s59*s8) + s50*(s24*s53*s8 + s21*(2*s27*s51 - 2*s28*s8)) + s40*(s64 + s13*s23*s9) + s42*(s65 + s18*s23*s9) + s39*(s66 + s23*s3*s9);
				buffer__[8*i+6] = s13*s18*s23*s3*s30 + s13*s24*s32 + (s13*s24*s33*s34)/2. + (s13*s24*s37*s38)/2. + (s13*s24*s44*s45)/2. + s35*((s13*s24*s36)/2. + s46) + s47*(2*s13*s23*s4 + s13*s24*s49) + s57*(2*s13*s19*s23 + s13*s24*s59) + s54*(s13*s24*s56 + s21*(-2*s13*s28 + 2*s27*s62)) + s31*(s14*s18*s23 + s65) + s29*(s14*s23*s3 + s66) + s13*s23*s3*s39*s8 + s13*s18*s23*s42*s8 + s40*(s67 + s14*s23*s8) + s50*(s13*s24*s53 + 2*s13*s23*s9);
				buffer__[8*i+7] = s13*s18*s23*s29*s3 + s18*s24*s32 + (s18*s24*s33*s34)/2. + (s18*s24*s35*s36)/2. + (s18*s24*s44*s45)/2. + s37*((s18*s24*s38)/2. + s46) + s47*(2*s18*s23*s4 + s18*s24*s49) + s54*(2*s14*s18*s23 + s18*s24*s56) + s57*(s18*s24*s59 + s21*(-2*s18*s28 + 2*s27*s63)) + s31*(s13*s19*s23 + s64) + s30*(s19*s23*s3 + s66) + s18*s23*s3*s39*s8 + s13*s18*s23*s40*s8 + s42*(s67 + s19*s23*s8) + s50*(s18*s24*s53 + 2*s18*s23*s9);
			}
		}

        ptoc(ClassName()+"::DFarToHulls");
        
    }

	}; // SimplicialMeshDetails<1,4,Real,Int>

//----------------------------------------------------------------------------------------------

    template<typename Real, typename Int>
    struct SimplicialMeshDetails<2,3,Real,Int>
    {
	private:

		const Int thread_count = 1;

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(std::max(static_cast<Int>(1),thread_count_))
		{}
	
		inline Int ThreadCount() const
		{
			return thread_count;
		}
	
		inline std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<2,3,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void ComputeNearFarData( 
		const Tensor2<Real,Int> & V_coords,
		const Tensor2<Int ,Int> & simplices,
			  Tensor2<Real,Int> & P_near,
			  Tensor2<Real,Int> & P_far
	) const
    {
        ptic(ClassName()+"::ComputeNearFarData");
        
        //Int size       = 3;
        //Int amb_dim    = 3;
        //Int dom_dim    = 2;

        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			ptr<Real> V_coords__      = V_coords.data();	
			ptr<Int>  simplices__     = simplices.data();

			Real hull    [3][3];
			Real df      [3][2];
			Real dfdagger[2][3];
			Real g       [2][2];
			Real ginv    [2][2];

			Int simplex  [3];
			
			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {
				Real * restrict const near = P_near.data(i);                    
				Real * restrict const far  = P_far.data(i);   
            
				simplex[0] = simplices__[3*i +0];
				simplex[1] = simplices__[3*i +1];
				simplex[2] = simplices__[3*i +2];

				near[1] = hull[0][0] = V_coords__[3*simplex[0]+0];
				near[2] = hull[0][1] = V_coords__[3*simplex[0]+1];
				near[3] = hull[0][2] = V_coords__[3*simplex[0]+2];
				near[4] = hull[1][0] = V_coords__[3*simplex[1]+0];
				near[5] = hull[1][1] = V_coords__[3*simplex[1]+1];
				near[6] = hull[1][2] = V_coords__[3*simplex[1]+2];
				near[7] = hull[2][0] = V_coords__[3*simplex[2]+0];
				near[8] = hull[2][1] = V_coords__[3*simplex[2]+1];
				near[9] = hull[2][2] = V_coords__[3*simplex[2]+2];

				far[1] = static_cast<Real>(0.3333333333333333) * ( hull[0][0] + hull[1][0] + hull[2][0] );
				far[2] = static_cast<Real>(0.3333333333333333) * ( hull[0][1] + hull[1][1] + hull[2][1] );
				far[3] = static_cast<Real>(0.3333333333333333) * ( hull[0][2] + hull[1][2] + hull[2][2] );

				df[0][0] = hull[1][0] - hull[0][0];
				df[0][1] = hull[2][0] - hull[0][0];
				df[1][0] = hull[1][1] - hull[0][1];
				df[1][1] = hull[2][1] - hull[0][1];
				df[2][0] = hull[1][2] - hull[0][2];
				df[2][1] = hull[2][2] - hull[0][2];

				g[0][0] = df[0][0] * df[0][0] + df[1][0] * df[1][0] + df[2][0] * df[2][0];
				g[0][1] = df[0][0] * df[0][1] + df[1][0] * df[1][1] + df[2][0] * df[2][1];
				g[1][0] = df[0][1] * df[0][0] + df[1][1] * df[1][0] + df[2][1] * df[2][0];
				g[1][1] = df[0][1] * df[0][1] + df[1][1] * df[1][1] + df[2][1] * df[2][1];

                Real det = g[0][0] * g[1][1] - g[0][1] * g[0][1];

				near[0] = far[0] = sqrt( fabs(det) ) * static_cast<Real>(0.5);

                Real invdet = static_cast<Real>(1)/det;
                ginv[0][0] =  g[1][1] * invdet;
                ginv[0][1] = -g[0][1] * invdet;
                ginv[1][1] =  g[0][0] * invdet;
                
                //  dfdagger = g^{-1} * df^T (2 x 3 matrix)
				dfdagger[0][0] = ginv[0][0] * df[0][0] + ginv[0][1] * df[0][1];
				dfdagger[0][1] = ginv[0][0] * df[1][0] + ginv[0][1] * df[1][1];
				dfdagger[0][2] = ginv[0][0] * df[2][0] + ginv[0][1] * df[2][1];
				dfdagger[1][0] = ginv[0][1] * df[0][0] + ginv[1][1] * df[0][1];
				dfdagger[1][1] = ginv[0][1] * df[1][0] + ginv[1][1] * df[1][1];
				dfdagger[1][2] = ginv[0][1] * df[2][0] + ginv[1][1] * df[2][1];
            
				near[10] = far[ 4]  = static_cast<Real>(1) - df[0][0] * dfdagger[0][0] - df[0][1] * dfdagger[1][0];
				near[11] = far[ 5]  =    - df[0][0] * dfdagger[0][1] - df[0][1] * dfdagger[1][1];
				near[12] = far[ 6]  =    - df[0][0] * dfdagger[0][2] - df[0][1] * dfdagger[1][2];
				near[13] = far[ 7]  = static_cast<Real>(1) - df[1][0] * dfdagger[0][1] - df[1][1] * dfdagger[1][1];
				near[14] = far[ 8]  =    - df[1][0] * dfdagger[0][2] - df[1][1] * dfdagger[1][2];
				near[15] = far[ 9]  = static_cast<Real>(1) - df[2][0] * dfdagger[0][2] - df[2][1] * dfdagger[1][2];

            } // for( Int i = i_begin; i < i_end; ++i )

        } // #pragma omp parallel num_threads( ThreadCount() )

        ptoc(ClassName()+"::ComputeNearFarData");
    }

	void ComputeNearFarDataOps(
		const Tensor2<Real,Int> & V_coords,
        const Tensor2<Int ,Int> & simplices,
		      Tensor2<Real,Int> & P_coords,
		      Tensor3<Real,Int> & P_hull_coords,
		      Tensor2<Real,Int> & P_near,
		      Tensor2<Real,Int> & P_far,
		Sparse::MatrixCSR<Real,Int,Int> & DiffOp,
		Sparse::MatrixCSR<Real,Int,Int> & AvOp 
	) const
    {
        ptic(ClassName()+"::ComputeNearFarDataOps");

        //Int size       = 3;
        //Int amb_dim    = 3;
        //Int dom_dim    = 2;

        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );

        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			mut<Int>  AvOp_outer = AvOp.Outer().data();
			mut<Int>  AvOp_inner = AvOp.Inner().data();
			mut<Real> AvOp_value = AvOp.Values().data();

			mut<Int>  DiffOp_outer = DiffOp.Outer().data();
			mut<Int>  DiffOp_inner = DiffOp.Inner().data();
			mut<Real> DiffOp_value = DiffOp.Values().data();

			ptr<Real> V_coords__      = V_coords.data();
			
			ptr<Int>  simplices__     = simplices.data();
		    mut<Real> P_hull_coords__ = P_hull_coords.data();
			mut<Real> P_coords__      = P_coords.data();

			Real df       [3][2];
			Real dfdagger [2][3];
			Real g        [2][2];
			Real ginv     [2][2];

			Int simplex        [3];
            Int sorted_simplex [3];

			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {

				mut<Real> near = P_near.data(i);                    
				mut<Real> far  = P_far.data(i);

				simplex[0] = sorted_simplex[0] = simplices__[3*i +0];
				simplex[1] = sorted_simplex[1] = simplices__[3*i +1];
				simplex[2] = sorted_simplex[2] = simplices__[3*i +2];
                  
                // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
                std::sort( sorted_simplex, sorted_simplex + 3 );

				AvOp_outer[i+1] = (i+1) * 3;                      
				AvOp_inner[3*i+0] = sorted_simplex[0];
				AvOp_inner[3*i+1] = sorted_simplex[1];
				AvOp_inner[3*i+2] = sorted_simplex[2];

				AvOp_value[3*i+0] = 0.3333333333333333;
				AvOp_value[3*i+1] = 0.3333333333333333;
				AvOp_value[3*i+2] = 0.3333333333333333;

				DiffOp_outer[3*i+0] = (3 * i + 0) * 3;
				DiffOp_outer[3*i+1] = (3 * i + 1) * 3;
				DiffOp_outer[3*i+2] = (3 * i + 2) * 3;

				DiffOp_inner[(i * 3 + 0) * 3 + 0 ] = sorted_simplex[0];
				DiffOp_inner[(i * 3 + 0) * 3 + 1 ] = sorted_simplex[1];
				DiffOp_inner[(i * 3 + 0) * 3 + 2 ] = sorted_simplex[2];
				DiffOp_inner[(i * 3 + 1) * 3 + 0 ] = sorted_simplex[0];
				DiffOp_inner[(i * 3 + 1) * 3 + 1 ] = sorted_simplex[1];
				DiffOp_inner[(i * 3 + 1) * 3 + 2 ] = sorted_simplex[2];
				DiffOp_inner[(i * 3 + 2) * 3 + 0 ] = sorted_simplex[0];
				DiffOp_inner[(i * 3 + 2) * 3 + 1 ] = sorted_simplex[1];
				DiffOp_inner[(i * 3 + 2) * 3 + 2 ] = sorted_simplex[2];

				near[1] = P_hull_coords__[9*i+0] = V_coords__[3*simplex[0]+0];
				near[2] = P_hull_coords__[9*i+1] = V_coords__[3*simplex[0]+1];
				near[3] = P_hull_coords__[9*i+2] = V_coords__[3*simplex[0]+2];
				near[4] = P_hull_coords__[9*i+3] = V_coords__[3*simplex[1]+0];
				near[5] = P_hull_coords__[9*i+4] = V_coords__[3*simplex[1]+1];
				near[6] = P_hull_coords__[9*i+5] = V_coords__[3*simplex[1]+2];
				near[7] = P_hull_coords__[9*i+6] = V_coords__[3*simplex[2]+0];
				near[8] = P_hull_coords__[9*i+7] = V_coords__[3*simplex[2]+1];
				near[9] = P_hull_coords__[9*i+8] = V_coords__[3*simplex[2]+2];

				far[1] = P_coords__[3*i+0] = 0.3333333333333333 * ( P_hull_coords__[9*i+0] + P_hull_coords__[9*i+3] + P_hull_coords__[9*i+6] );
				far[2] = P_coords__[3*i+1] = 0.3333333333333333 * ( P_hull_coords__[9*i+1] + P_hull_coords__[9*i+4] + P_hull_coords__[9*i+7] );
				far[3] = P_coords__[3*i+2] = 0.3333333333333333 * ( P_hull_coords__[9*i+2] + P_hull_coords__[9*i+5] + P_hull_coords__[9*i+8] );

				df[0][0] = V_coords__[3*sorted_simplex[1]+0] - V_coords__[3*sorted_simplex[0]+0];
				df[0][1] = V_coords__[3*sorted_simplex[2]+0] - V_coords__[3*sorted_simplex[0]+0];
				df[1][0] = V_coords__[3*sorted_simplex[1]+1] - V_coords__[3*sorted_simplex[0]+1];
				df[1][1] = V_coords__[3*sorted_simplex[2]+1] - V_coords__[3*sorted_simplex[0]+1];
				df[2][0] = V_coords__[3*sorted_simplex[1]+2] - V_coords__[3*sorted_simplex[0]+2];
				df[2][1] = V_coords__[3*sorted_simplex[2]+2] - V_coords__[3*sorted_simplex[0]+2];

				g[0][0] = df[0][0] * df[0][0] + df[1][0] * df[1][0] + df[2][0] * df[2][0];
				g[0][1] = df[0][0] * df[0][1] + df[1][0] * df[1][1] + df[2][0] * df[2][1];
				g[1][0] = df[0][1] * df[0][0] + df[1][1] * df[1][0] + df[2][1] * df[2][0];
				g[1][1] = df[0][1] * df[0][1] + df[1][1] * df[1][1] + df[2][1] * df[2][1];

                Real det = g[0][0] * g[1][1] - g[0][1] * g[0][1];

                near[0] = far[0] = sqrt( fabs(det) ) * static_cast<Real>(0.5);

                Real invdet = static_cast<Real>(1.0)/det;
                ginv[0][0] =  g[1][1] * invdet;
                ginv[0][1] = -g[0][1] * invdet;
                ginv[1][1] =  g[0][0] * invdet;
                
                //  dfdagger = g^{-1} * df^T (2 x 3 matrix)
				dfdagger[0][0] = ginv[0][0] * df[0][0] + ginv[0][1] * df[0][1];
				dfdagger[0][1] = ginv[0][0] * df[1][0] + ginv[0][1] * df[1][1];
				dfdagger[0][2] = ginv[0][0] * df[2][0] + ginv[0][1] * df[2][1];
				dfdagger[1][0] = ginv[0][1] * df[0][0] + ginv[1][1] * df[0][1];
				dfdagger[1][1] = ginv[0][1] * df[1][0] + ginv[1][1] * df[1][1];
				dfdagger[1][2] = ginv[0][1] * df[2][0] + ginv[1][1] * df[2][1];
            
				near[10] = far[ 4]  = static_cast<Real>(1.0) - df[0][0] * dfdagger[0][0] - df[0][1] * dfdagger[1][0];
				near[11] = far[ 5]  =    - df[0][0] * dfdagger[0][1] - df[0][1] * dfdagger[1][1];
				near[12] = far[ 6]  =    - df[0][0] * dfdagger[0][2] - df[0][1] * dfdagger[1][2];
				near[13] = far[ 7]  = static_cast<Real>(1.0) - df[1][0] * dfdagger[0][1] - df[1][1] * dfdagger[1][1];
				near[14] = far[ 8]  =    - df[1][0] * dfdagger[0][2] - df[1][1] * dfdagger[1][2];
				near[15] = far[ 9]  = static_cast<Real>(1.0) - df[2][0] * dfdagger[0][2] - df[2][1] * dfdagger[1][2];

                // derivative operator  (3 x 3 matrix)

                mut<Real> Df = &DiffOp_value[ 9 * i ];

				Df[ 0] = - dfdagger[0][0] - dfdagger[1][0];
				Df[ 1] =   dfdagger[0][0];
				Df[ 2] =   dfdagger[1][0];
				Df[ 3] = - dfdagger[0][1] - dfdagger[1][1];
				Df[ 4] =   dfdagger[0][1];
				Df[ 5] =   dfdagger[1][1];
				Df[ 6] = - dfdagger[0][2] - dfdagger[1][2];
				Df[ 7] =   dfdagger[0][2];
				Df[ 8] =   dfdagger[1][2];

            }
        }
        ptoc(ClassName()+"::ComputeNearFarDataOps");
    }

    void DNearToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_near, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DNearToHulls");

        if( P_D_near.Dimension(1) != 16 )
        {
            eprint("in DNearToHulls: P_D_near.Dimension(1) != 16. Aborting");
        }
        
		ptr<Real> V_coords__  = V_coords.data();
        ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_near__  = P_D_near.data();
        mut<Real> buffer__    = buffer.data();

		if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[3*simplices__[3*i+0]+2];
				const Real s1 = V_coords__[3*simplices__[3*i+1]+1];
				const Real s2 = -(s0*s1);
				const Real s3 = V_coords__[3*simplices__[3*i+0]+1];
				const Real s4 = V_coords__[3*simplices__[3*i+1]+2];
				const Real s5 = s3*s4;
				const Real s6 = V_coords__[3*simplices__[3*i+2]+1];
				const Real s7 = s0*s6;
				const Real s8 = -(s4*s6);
				const Real s9 = V_coords__[3*simplices__[3*i+2]+2];
				const Real s10 = -(s3*s9);
				const Real s11 = s1*s9;
				const Real s12 = s10 + s11 + s2 + s5 + s7 + s8;
				const Real s13 = s12*s12;
				const Real s14 = V_coords__[3*simplices__[3*i+2]+0];
				const Real s15 = V_coords__[3*simplices__[3*i+0]+0];
				const Real s16 = V_coords__[3*simplices__[3*i+1]+0];
				const Real s17 = -(s16*s3);
				const Real s18 = s1*s15;
				const Real s19 = s14*s3;
				const Real s20 = -(s1*s14);
				const Real s21 = -(s15*s6);
				const Real s22 = s16*s6;
				const Real s23 = s17 + s18 + s19 + s20 + s21 + s22;
				const Real s24 = s23*s23;
				const Real s25 = s0*s16;
				const Real s26 = -(s15*s4);
				const Real s27 = -(s0*s14);
				const Real s28 = s14*s4;
				const Real s29 = s15*s9;
				const Real s30 = -(s16*s9);
				const Real s31 = s25 + s26 + s27 + s28 + s29 + s30;
				const Real s32 = s31*s31;
				const Real s33 = s13 + s24 + s32;
				const Real s34 = sqrt(s33);
				const Real s35 = s33*s34;
				const Real s36 = 1/s35;
				const Real s37 = -s6;
				const Real s38 = s1 + s37;
				const Real s39 = 2*s23*s38;
				const Real s40 = -s4;
				const Real s41 = s40 + s9;
				const Real s42 = 2*s31*s41;
				const Real s43 = s39 + s42;
				const Real s44 = 1/s34;
				const Real s45 = P_D_near__[16*i+13];
				const Real s46 = P_D_near__[16*i+0];
				const Real s47 = -s16;
				const Real s48 = s14 + s47;
				const Real s49 = 2*s23*s48;
				const Real s50 = -s9;
				const Real s51 = s4 + s50;
				const Real s52 = 2*s12*s51;
				const Real s53 = s49 + s52;
				const Real s54 = P_D_near__[16*i+1];
				const Real s55 = P_D_near__[16*i+3];
				const Real s56 = P_D_near__[16*i+4];
				const Real s57 = P_D_near__[16*i+5];
				const Real s58 = P_D_near__[16*i+6];
				const Real s59 = P_D_near__[16*i+7];
				const Real s60 = P_D_near__[16*i+8];
				const Real s61 = P_D_near__[16*i+9];
				const Real s62 = P_D_near__[16*i+15];
				const Real s63 = P_D_near__[16*i+14];
				const Real s64 = P_D_near__[16*i+11];
				const Real s65 = P_D_near__[16*i+12];
				const Real s66 = P_D_near__[16*i+10];
				const Real s67 = P_D_near__[16*i+2];
				const Real s68 = s34/2.;
				const Real s69 = -s14;
				const Real s70 = s16 + s69;
				const Real s71 = 2*s31*s70;
				const Real s72 = -s1;
				const Real s73 = s6 + s72;
				const Real s74 = 2*s12*s73;
				const Real s75 = s71 + s74;
				const Real s76 = -s3;
				const Real s77 = s6 + s76;
				const Real s78 = 2*s23*s77;
				const Real s79 = s0 + s50;
				const Real s80 = 2*s31*s79;
				const Real s81 = s78 + s80;
				const Real s82 = s15 + s69;
				const Real s83 = 2*s23*s82;
				const Real s84 = -s0;
				const Real s85 = s84 + s9;
				const Real s86 = 2*s12*s85;
				const Real s87 = s83 + s86;
				const Real s88 = -s15;
				const Real s89 = s14 + s88;
				const Real s90 = 2*s31*s89;
				const Real s91 = s3 + s37;
				const Real s92 = 2*s12*s91;
				const Real s93 = s90 + s92;
				const Real s94 = s3 + s72;
				const Real s95 = 2*s23*s94;
				const Real s96 = s4 + s84;
				const Real s97 = 2*s31*s96;
				const Real s98 = s95 + s97;
				const Real s99 = s16 + s88;
				const Real s100 = 2*s23*s99;
				const Real s101 = s0 + s40;
				const Real s102 = 2*s101*s12;
				const Real s103 = s100 + s102;
				const Real s104 = s15 + s47;
				const Real s105 = 2*s104*s31;
				const Real s106 = s1 + s76;
				const Real s107 = 2*s106*s12;
				const Real s108 = s105 + s107;
				buffer__[9*i+0] += (-0.25*(s32*s36*s43) + s31*s41*s44)*s45 + (s43*s44*s46)/4. + (s0*s43*s44*s55)/4. + (s16*s43*s44*s56)/4. + (s1*s43*s44*s57)/4. + (s4*s43*s44*s58)/4. + (s14*s43*s44*s59)/4. + (s43*s44*s6*s60)/4. + (-0.25*(s24*s36*s43) + s23*s38*s44)*s62 + (-0.25*(s23*s31*s36*s43) + (s31*s38*s44)/2. + (s23*s41*s44)/2.)*s63 + (-0.25*(s12*s31*s36*s43) + (s12*s41*s44)/2.)*s64 + (-0.25*(s12*s23*s36*s43) + (s12*s38*s44)/2.)*s65 - (s13*s36*s43*s66)/4. + (s3*s43*s44*s67)/4. + s54*((s15*s43*s44)/4. + s68) + (s43*s44*s61*s9)/4.;
				buffer__[9*i+1] += -0.25*(s32*s36*s45*s53) + (s44*s46*s53)/4. + (s15*s44*s53*s54)/4. + (s0*s44*s53*s55)/4. + (s16*s44*s53*s56)/4. + (s1*s44*s53*s57)/4. + (s4*s44*s53*s58)/4. + (s14*s44*s53*s59)/4. + (s44*s53*s6*s60)/4. + (s23*s44*s48 - (s24*s36*s53)/4.)*s62 + ((s31*s44*s48)/2. - (s23*s31*s36*s53)/4.)*s63 + ((s31*s44*s51)/2. - (s12*s31*s36*s53)/4.)*s64 + ((s12*s44*s48)/2. + (s23*s44*s51)/2. - (s12*s23*s36*s53)/4.)*s65 + (s12*s44*s51 - (s13*s36*s53)/4.)*s66 + s67*((s3*s44*s53)/4. + s68) + (s44*s53*s61*s9)/4.;
				buffer__[9*i+2] += (s44*s46*s75)/4. + (s15*s44*s54*s75)/4. + (s16*s44*s56*s75)/4. + (s1*s44*s57*s75)/4. + (s4*s44*s58*s75)/4. + (s14*s44*s59*s75)/4. + (s44*s6*s60*s75)/4. - (s24*s36*s62*s75)/4. + (s3*s44*s67*s75)/4. + s66*(s12*s44*s73 - (s13*s36*s75)/4.) + s65*((s23*s44*s73)/2. - (s12*s23*s36*s75)/4.) + s64*((s12*s44*s70)/2. + (s31*s44*s73)/2. - (s12*s31*s36*s75)/4.) + s63*((s23*s44*s70)/2. - (s23*s31*s36*s75)/4.) + s45*(s31*s44*s70 - (s32*s36*s75)/4.) + s55*(s68 + (s0*s44*s75)/4.) + (s44*s61*s75*s9)/4.;
				buffer__[9*i+3] += (s44*s46*s81)/4. + (s15*s44*s54*s81)/4. + (s0*s44*s55*s81)/4. + (s1*s44*s57*s81)/4. + (s4*s44*s58*s81)/4. + (s14*s44*s59*s81)/4. + (s44*s6*s60*s81)/4. - (s13*s36*s66*s81)/4. + (s3*s44*s67*s81)/4. + s65*((s12*s44*s77)/2. - (s12*s23*s36*s81)/4.) + s62*(s23*s44*s77 - (s24*s36*s81)/4.) + s64*((s12*s44*s79)/2. - (s12*s31*s36*s81)/4.) + s63*((s31*s44*s77)/2. + (s23*s44*s79)/2. - (s23*s31*s36*s81)/4.) + s45*(s31*s44*s79 - (s32*s36*s81)/4.) + s56*(s68 + (s16*s44*s81)/4.) + (s44*s61*s81*s9)/4.;
				buffer__[9*i+4] += -0.25*(s32*s36*s45*s87) + (s44*s46*s87)/4. + (s15*s44*s54*s87)/4. + (s0*s44*s55*s87)/4. + (s16*s44*s56*s87)/4. + (s4*s44*s58*s87)/4. + (s14*s44*s59*s87)/4. + (s44*s6*s60*s87)/4. + (s3*s44*s67*s87)/4. + s66*(s12*s44*s85 - (s13*s36*s87)/4.) + s65*((s12*s44*s82)/2. + (s23*s44*s85)/2. - (s12*s23*s36*s87)/4.) + s62*(s23*s44*s82 - (s24*s36*s87)/4.) + s64*((s31*s44*s85)/2. - (s12*s31*s36*s87)/4.) + s63*((s31*s44*s82)/2. - (s23*s31*s36*s87)/4.) + s57*(s68 + (s1*s44*s87)/4.) + (s44*s61*s87*s9)/4.;
				buffer__[9*i+5] += (s44*s46*s93)/4. + (s15*s44*s54*s93)/4. + (s0*s44*s55*s93)/4. + (s16*s44*s56*s93)/4. + (s1*s44*s57*s93)/4. + (s14*s44*s59*s93)/4. + (s44*s6*s60*s93)/4. - (s24*s36*s62*s93)/4. + (s3*s44*s67*s93)/4. + (s44*s61*s9*s93)/4. + s66*(s12*s44*s91 - (s13*s36*s93)/4.) + s65*((s23*s44*s91)/2. - (s12*s23*s36*s93)/4.) + s64*((s12*s44*s89)/2. + (s31*s44*s91)/2. - (s12*s31*s36*s93)/4.) + s63*((s23*s44*s89)/2. - (s23*s31*s36*s93)/4.) + s45*(s31*s44*s89 - (s32*s36*s93)/4.) + s58*(s68 + (s4*s44*s93)/4.);
				buffer__[9*i+6] += (s44*s46*s98)/4. + (s15*s44*s54*s98)/4. + (s0*s44*s55*s98)/4. + (s16*s44*s56*s98)/4. + (s1*s44*s57*s98)/4. + (s4*s44*s58*s98)/4. + (s44*s6*s60*s98)/4. - (s13*s36*s66*s98)/4. + (s3*s44*s67*s98)/4. + (s44*s61*s9*s98)/4. + s65*((s12*s44*s94)/2. - (s12*s23*s36*s98)/4.) + s62*(s23*s44*s94 - (s24*s36*s98)/4.) + s64*((s12*s44*s96)/2. - (s12*s31*s36*s98)/4.) + s63*((s31*s44*s94)/2. + (s23*s44*s96)/2. - (s23*s31*s36*s98)/4.) + s45*(s31*s44*s96 - (s32*s36*s98)/4.) + s59*(s68 + (s14*s44*s98)/4.);
				buffer__[9*i+7] += -0.25*(s103*s32*s36*s45) + (s103*s44*s46)/4. + (s103*s15*s44*s54)/4. + (s0*s103*s44*s55)/4. + (s103*s16*s44*s56)/4. + (s1*s103*s44*s57)/4. + (s103*s4*s44*s58)/4. + (s103*s14*s44*s59)/4. + (-0.25*(s103*s12*s31*s36) + (s101*s31*s44)/2.)*s64 + (-0.25*(s103*s13*s36) + s101*s12*s44)*s66 + (s103*s3*s44*s67)/4. + s60*((s103*s44*s6)/4. + s68) + (s103*s44*s61*s9)/4. + s65*(-0.25*(s103*s12*s23*s36) + (s101*s23*s44)/2. + (s12*s44*s99)/2.) + s62*(-0.25*(s103*s24*s36) + s23*s44*s99) + s63*(-0.25*(s103*s23*s31*s36) + (s31*s44*s99)/2.);
				buffer__[9*i+8] += (-0.25*(s108*s32*s36) + s104*s31*s44)*s45 + (s108*s44*s46)/4. + (s108*s15*s44*s54)/4. + (s0*s108*s44*s55)/4. + (s108*s16*s44*s56)/4. + (s1*s108*s44*s57)/4. + (s108*s4*s44*s58)/4. + (s108*s14*s44*s59)/4. + (s108*s44*s6*s60)/4. - (s108*s24*s36*s62)/4. + (-0.25*(s108*s23*s31*s36) + (s104*s23*s44)/2.)*s63 + (-0.25*(s108*s12*s31*s36) + (s104*s12*s44)/2. + (s106*s31*s44)/2.)*s64 + (-0.25*(s108*s12*s23*s36) + (s106*s23*s44)/2.)*s65 + (-0.25*(s108*s13*s36) + s106*s12*s44)*s66 + (s108*s3*s44*s67)/4. + s61*(s68 + (s108*s44*s9)/4.);
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[3*simplices__[3*i+0]+2];
				const Real s1 = V_coords__[3*simplices__[3*i+1]+1];
				const Real s2 = -(s0*s1);
				const Real s3 = V_coords__[3*simplices__[3*i+0]+1];
				const Real s4 = V_coords__[3*simplices__[3*i+1]+2];
				const Real s5 = s3*s4;
				const Real s6 = V_coords__[3*simplices__[3*i+2]+1];
				const Real s7 = s0*s6;
				const Real s8 = -(s4*s6);
				const Real s9 = V_coords__[3*simplices__[3*i+2]+2];
				const Real s10 = -(s3*s9);
				const Real s11 = s1*s9;
				const Real s12 = s10 + s11 + s2 + s5 + s7 + s8;
				const Real s13 = s12*s12;
				const Real s14 = V_coords__[3*simplices__[3*i+2]+0];
				const Real s15 = V_coords__[3*simplices__[3*i+0]+0];
				const Real s16 = V_coords__[3*simplices__[3*i+1]+0];
				const Real s17 = -(s16*s3);
				const Real s18 = s1*s15;
				const Real s19 = s14*s3;
				const Real s20 = -(s1*s14);
				const Real s21 = -(s15*s6);
				const Real s22 = s16*s6;
				const Real s23 = s17 + s18 + s19 + s20 + s21 + s22;
				const Real s24 = s23*s23;
				const Real s25 = s0*s16;
				const Real s26 = -(s15*s4);
				const Real s27 = -(s0*s14);
				const Real s28 = s14*s4;
				const Real s29 = s15*s9;
				const Real s30 = -(s16*s9);
				const Real s31 = s25 + s26 + s27 + s28 + s29 + s30;
				const Real s32 = s31*s31;
				const Real s33 = s13 + s24 + s32;
				const Real s34 = sqrt(s33);
				const Real s35 = s33*s34;
				const Real s36 = 1/s35;
				const Real s37 = -s6;
				const Real s38 = s1 + s37;
				const Real s39 = 2*s23*s38;
				const Real s40 = -s4;
				const Real s41 = s40 + s9;
				const Real s42 = 2*s31*s41;
				const Real s43 = s39 + s42;
				const Real s44 = 1/s34;
				const Real s45 = P_D_near__[16*i+13];
				const Real s46 = P_D_near__[16*i+0];
				const Real s47 = -s16;
				const Real s48 = s14 + s47;
				const Real s49 = 2*s23*s48;
				const Real s50 = -s9;
				const Real s51 = s4 + s50;
				const Real s52 = 2*s12*s51;
				const Real s53 = s49 + s52;
				const Real s54 = P_D_near__[16*i+1];
				const Real s55 = P_D_near__[16*i+3];
				const Real s56 = P_D_near__[16*i+4];
				const Real s57 = P_D_near__[16*i+5];
				const Real s58 = P_D_near__[16*i+6];
				const Real s59 = P_D_near__[16*i+7];
				const Real s60 = P_D_near__[16*i+8];
				const Real s61 = P_D_near__[16*i+9];
				const Real s62 = P_D_near__[16*i+15];
				const Real s63 = P_D_near__[16*i+14];
				const Real s64 = P_D_near__[16*i+11];
				const Real s65 = P_D_near__[16*i+12];
				const Real s66 = P_D_near__[16*i+10];
				const Real s67 = P_D_near__[16*i+2];
				const Real s68 = s34/2.;
				const Real s69 = -s14;
				const Real s70 = s16 + s69;
				const Real s71 = 2*s31*s70;
				const Real s72 = -s1;
				const Real s73 = s6 + s72;
				const Real s74 = 2*s12*s73;
				const Real s75 = s71 + s74;
				const Real s76 = -s3;
				const Real s77 = s6 + s76;
				const Real s78 = 2*s23*s77;
				const Real s79 = s0 + s50;
				const Real s80 = 2*s31*s79;
				const Real s81 = s78 + s80;
				const Real s82 = s15 + s69;
				const Real s83 = 2*s23*s82;
				const Real s84 = -s0;
				const Real s85 = s84 + s9;
				const Real s86 = 2*s12*s85;
				const Real s87 = s83 + s86;
				const Real s88 = -s15;
				const Real s89 = s14 + s88;
				const Real s90 = 2*s31*s89;
				const Real s91 = s3 + s37;
				const Real s92 = 2*s12*s91;
				const Real s93 = s90 + s92;
				const Real s94 = s3 + s72;
				const Real s95 = 2*s23*s94;
				const Real s96 = s4 + s84;
				const Real s97 = 2*s31*s96;
				const Real s98 = s95 + s97;
				const Real s99 = s16 + s88;
				const Real s100 = 2*s23*s99;
				const Real s101 = s0 + s40;
				const Real s102 = 2*s101*s12;
				const Real s103 = s100 + s102;
				const Real s104 = s15 + s47;
				const Real s105 = 2*s104*s31;
				const Real s106 = s1 + s76;
				const Real s107 = 2*s106*s12;
				const Real s108 = s105 + s107;
				buffer__[9*i+0] = (-0.25*(s32*s36*s43) + s31*s41*s44)*s45 + (s43*s44*s46)/4. + (s0*s43*s44*s55)/4. + (s16*s43*s44*s56)/4. + (s1*s43*s44*s57)/4. + (s4*s43*s44*s58)/4. + (s14*s43*s44*s59)/4. + (s43*s44*s6*s60)/4. + (-0.25*(s24*s36*s43) + s23*s38*s44)*s62 + (-0.25*(s23*s31*s36*s43) + (s31*s38*s44)/2. + (s23*s41*s44)/2.)*s63 + (-0.25*(s12*s31*s36*s43) + (s12*s41*s44)/2.)*s64 + (-0.25*(s12*s23*s36*s43) + (s12*s38*s44)/2.)*s65 - (s13*s36*s43*s66)/4. + (s3*s43*s44*s67)/4. + s54*((s15*s43*s44)/4. + s68) + (s43*s44*s61*s9)/4.;
				buffer__[9*i+1] = -0.25*(s32*s36*s45*s53) + (s44*s46*s53)/4. + (s15*s44*s53*s54)/4. + (s0*s44*s53*s55)/4. + (s16*s44*s53*s56)/4. + (s1*s44*s53*s57)/4. + (s4*s44*s53*s58)/4. + (s14*s44*s53*s59)/4. + (s44*s53*s6*s60)/4. + (s23*s44*s48 - (s24*s36*s53)/4.)*s62 + ((s31*s44*s48)/2. - (s23*s31*s36*s53)/4.)*s63 + ((s31*s44*s51)/2. - (s12*s31*s36*s53)/4.)*s64 + ((s12*s44*s48)/2. + (s23*s44*s51)/2. - (s12*s23*s36*s53)/4.)*s65 + (s12*s44*s51 - (s13*s36*s53)/4.)*s66 + s67*((s3*s44*s53)/4. + s68) + (s44*s53*s61*s9)/4.;
				buffer__[9*i+2] = (s44*s46*s75)/4. + (s15*s44*s54*s75)/4. + (s16*s44*s56*s75)/4. + (s1*s44*s57*s75)/4. + (s4*s44*s58*s75)/4. + (s14*s44*s59*s75)/4. + (s44*s6*s60*s75)/4. - (s24*s36*s62*s75)/4. + (s3*s44*s67*s75)/4. + s66*(s12*s44*s73 - (s13*s36*s75)/4.) + s65*((s23*s44*s73)/2. - (s12*s23*s36*s75)/4.) + s64*((s12*s44*s70)/2. + (s31*s44*s73)/2. - (s12*s31*s36*s75)/4.) + s63*((s23*s44*s70)/2. - (s23*s31*s36*s75)/4.) + s45*(s31*s44*s70 - (s32*s36*s75)/4.) + s55*(s68 + (s0*s44*s75)/4.) + (s44*s61*s75*s9)/4.;
				buffer__[9*i+3] = (s44*s46*s81)/4. + (s15*s44*s54*s81)/4. + (s0*s44*s55*s81)/4. + (s1*s44*s57*s81)/4. + (s4*s44*s58*s81)/4. + (s14*s44*s59*s81)/4. + (s44*s6*s60*s81)/4. - (s13*s36*s66*s81)/4. + (s3*s44*s67*s81)/4. + s65*((s12*s44*s77)/2. - (s12*s23*s36*s81)/4.) + s62*(s23*s44*s77 - (s24*s36*s81)/4.) + s64*((s12*s44*s79)/2. - (s12*s31*s36*s81)/4.) + s63*((s31*s44*s77)/2. + (s23*s44*s79)/2. - (s23*s31*s36*s81)/4.) + s45*(s31*s44*s79 - (s32*s36*s81)/4.) + s56*(s68 + (s16*s44*s81)/4.) + (s44*s61*s81*s9)/4.;
				buffer__[9*i+4] = -0.25*(s32*s36*s45*s87) + (s44*s46*s87)/4. + (s15*s44*s54*s87)/4. + (s0*s44*s55*s87)/4. + (s16*s44*s56*s87)/4. + (s4*s44*s58*s87)/4. + (s14*s44*s59*s87)/4. + (s44*s6*s60*s87)/4. + (s3*s44*s67*s87)/4. + s66*(s12*s44*s85 - (s13*s36*s87)/4.) + s65*((s12*s44*s82)/2. + (s23*s44*s85)/2. - (s12*s23*s36*s87)/4.) + s62*(s23*s44*s82 - (s24*s36*s87)/4.) + s64*((s31*s44*s85)/2. - (s12*s31*s36*s87)/4.) + s63*((s31*s44*s82)/2. - (s23*s31*s36*s87)/4.) + s57*(s68 + (s1*s44*s87)/4.) + (s44*s61*s87*s9)/4.;
				buffer__[9*i+5] = (s44*s46*s93)/4. + (s15*s44*s54*s93)/4. + (s0*s44*s55*s93)/4. + (s16*s44*s56*s93)/4. + (s1*s44*s57*s93)/4. + (s14*s44*s59*s93)/4. + (s44*s6*s60*s93)/4. - (s24*s36*s62*s93)/4. + (s3*s44*s67*s93)/4. + (s44*s61*s9*s93)/4. + s66*(s12*s44*s91 - (s13*s36*s93)/4.) + s65*((s23*s44*s91)/2. - (s12*s23*s36*s93)/4.) + s64*((s12*s44*s89)/2. + (s31*s44*s91)/2. - (s12*s31*s36*s93)/4.) + s63*((s23*s44*s89)/2. - (s23*s31*s36*s93)/4.) + s45*(s31*s44*s89 - (s32*s36*s93)/4.) + s58*(s68 + (s4*s44*s93)/4.);
				buffer__[9*i+6] = (s44*s46*s98)/4. + (s15*s44*s54*s98)/4. + (s0*s44*s55*s98)/4. + (s16*s44*s56*s98)/4. + (s1*s44*s57*s98)/4. + (s4*s44*s58*s98)/4. + (s44*s6*s60*s98)/4. - (s13*s36*s66*s98)/4. + (s3*s44*s67*s98)/4. + (s44*s61*s9*s98)/4. + s65*((s12*s44*s94)/2. - (s12*s23*s36*s98)/4.) + s62*(s23*s44*s94 - (s24*s36*s98)/4.) + s64*((s12*s44*s96)/2. - (s12*s31*s36*s98)/4.) + s63*((s31*s44*s94)/2. + (s23*s44*s96)/2. - (s23*s31*s36*s98)/4.) + s45*(s31*s44*s96 - (s32*s36*s98)/4.) + s59*(s68 + (s14*s44*s98)/4.);
				buffer__[9*i+7] = -0.25*(s103*s32*s36*s45) + (s103*s44*s46)/4. + (s103*s15*s44*s54)/4. + (s0*s103*s44*s55)/4. + (s103*s16*s44*s56)/4. + (s1*s103*s44*s57)/4. + (s103*s4*s44*s58)/4. + (s103*s14*s44*s59)/4. + (-0.25*(s103*s12*s31*s36) + (s101*s31*s44)/2.)*s64 + (-0.25*(s103*s13*s36) + s101*s12*s44)*s66 + (s103*s3*s44*s67)/4. + s60*((s103*s44*s6)/4. + s68) + (s103*s44*s61*s9)/4. + s65*(-0.25*(s103*s12*s23*s36) + (s101*s23*s44)/2. + (s12*s44*s99)/2.) + s62*(-0.25*(s103*s24*s36) + s23*s44*s99) + s63*(-0.25*(s103*s23*s31*s36) + (s31*s44*s99)/2.);
				buffer__[9*i+8] = (-0.25*(s108*s32*s36) + s104*s31*s44)*s45 + (s108*s44*s46)/4. + (s108*s15*s44*s54)/4. + (s0*s108*s44*s55)/4. + (s108*s16*s44*s56)/4. + (s1*s108*s44*s57)/4. + (s108*s4*s44*s58)/4. + (s108*s14*s44*s59)/4. + (s108*s44*s6*s60)/4. - (s108*s24*s36*s62)/4. + (-0.25*(s108*s23*s31*s36) + (s104*s23*s44)/2.)*s63 + (-0.25*(s108*s12*s31*s36) + (s104*s12*s44)/2. + (s106*s31*s44)/2.)*s64 + (-0.25*(s108*s12*s23*s36) + (s106*s23*s44)/2.)*s65 + (-0.25*(s108*s13*s36) + s106*s12*s44)*s66 + (s108*s3*s44*s67)/4. + s61*(s68 + (s108*s44*s9)/4.);
			}
		}

        ptoc(ClassName()+"::DNearToHulls");
        
    }

    void DFarToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_far, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DFarToHulls");

        if( P_D_far.Dimension(1) != 10 )
        {
            eprint("in DFarToHulls: P_D_far.Dimension(1) != 10. Aborting");
        }
        
		ptr<Real> V_coords__  = V_coords.data();
        ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_far__   = P_D_far.data();
              Real * restrict const buffer__    = buffer.data();

		if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[3*simplices__[3*i+0]+2];
				const Real s1 = V_coords__[3*simplices__[3*i+1]+1];
				const Real s2 = -(s0*s1);
				const Real s3 = V_coords__[3*simplices__[3*i+0]+1];
				const Real s4 = V_coords__[3*simplices__[3*i+1]+2];
				const Real s5 = s3*s4;
				const Real s6 = V_coords__[3*simplices__[3*i+2]+1];
				const Real s7 = s0*s6;
				const Real s8 = -(s4*s6);
				const Real s9 = V_coords__[3*simplices__[3*i+2]+2];
				const Real s10 = -(s3*s9);
				const Real s11 = s1*s9;
				const Real s12 = s10 + s11 + s2 + s5 + s7 + s8;
				const Real s13 = s12*s12;
				const Real s14 = V_coords__[3*simplices__[3*i+2]+0];
				const Real s15 = V_coords__[3*simplices__[3*i+0]+0];
				const Real s16 = V_coords__[3*simplices__[3*i+1]+0];
				const Real s17 = -(s16*s3);
				const Real s18 = s1*s15;
				const Real s19 = s14*s3;
				const Real s20 = -(s1*s14);
				const Real s21 = -(s15*s6);
				const Real s22 = s16*s6;
				const Real s23 = s17 + s18 + s19 + s20 + s21 + s22;
				const Real s24 = s23*s23;
				const Real s25 = s0*s16;
				const Real s26 = -(s15*s4);
				const Real s27 = -(s0*s14);
				const Real s28 = s14*s4;
				const Real s29 = s15*s9;
				const Real s30 = -(s16*s9);
				const Real s31 = s25 + s26 + s27 + s28 + s29 + s30;
				const Real s32 = s31*s31;
				const Real s33 = s13 + s24 + s32;
				const Real s34 = sqrt(s33);
				const Real s35 = s33*s34;
				const Real s36 = 1/s35;
				const Real s37 = -s6;
				const Real s38 = s1 + s37;
				const Real s39 = 2*s23*s38;
				const Real s40 = -s4;
				const Real s41 = s40 + s9;
				const Real s42 = 2*s31*s41;
				const Real s43 = s39 + s42;
				const Real s44 = 1/s34;
				const Real s45 = P_D_far__[10*i+7];
				const Real s46 = P_D_far__[10*i+0];
				const Real s47 = -s16;
				const Real s48 = s14 + s47;
				const Real s49 = 2*s23*s48;
				const Real s50 = -s9;
				const Real s51 = s4 + s50;
				const Real s52 = 2*s12*s51;
				const Real s53 = s49 + s52;
				const Real s54 = P_D_far__[10*i+1];
				const Real s55 = s14 + s15 + s16;
				const Real s56 = P_D_far__[10*i+3];
				const Real s57 = s0 + s4 + s9;
				const Real s58 = P_D_far__[10*i+9];
				const Real s59 = P_D_far__[10*i+8];
				const Real s60 = P_D_far__[10*i+5];
				const Real s61 = P_D_far__[10*i+6];
				const Real s62 = P_D_far__[10*i+4];
				const Real s63 = P_D_far__[10*i+2];
				const Real s64 = s1 + s3 + s6;
				const Real s65 = s34/6.;
				const Real s66 = -s14;
				const Real s67 = s16 + s66;
				const Real s68 = 2*s31*s67;
				const Real s69 = -s1;
				const Real s70 = s6 + s69;
				const Real s71 = 2*s12*s70;
				const Real s72 = s68 + s71;
				const Real s73 = -s3;
				const Real s74 = s6 + s73;
				const Real s75 = 2*s23*s74;
				const Real s76 = s0 + s50;
				const Real s77 = 2*s31*s76;
				const Real s78 = s75 + s77;
				const Real s79 = s15 + s66;
				const Real s80 = 2*s23*s79;
				const Real s81 = -s0;
				const Real s82 = s81 + s9;
				const Real s83 = 2*s12*s82;
				const Real s84 = s80 + s83;
				const Real s85 = -s15;
				const Real s86 = s14 + s85;
				const Real s87 = 2*s31*s86;
				const Real s88 = s3 + s37;
				const Real s89 = 2*s12*s88;
				const Real s90 = s87 + s89;
				const Real s91 = s3 + s69;
				const Real s92 = 2*s23*s91;
				const Real s93 = s4 + s81;
				const Real s94 = 2*s31*s93;
				const Real s95 = s92 + s94;
				const Real s96 = s16 + s85;
				const Real s97 = 2*s23*s96;
				const Real s98 = s0 + s40;
				const Real s99 = 2*s12*s98;
				const Real s100 = s97 + s99;
				const Real s101 = s15 + s47;
				const Real s102 = 2*s101*s31;
				const Real s103 = s1 + s73;
				const Real s104 = 2*s103*s12;
				const Real s105 = s102 + s104;
				buffer__[9*i+0] += (-0.25*(s32*s36*s43) + s31*s41*s44)*s45 + (s43*s44*s46)/4. + (s43*s44*s56*s57)/12. + (-0.25*(s24*s36*s43) + s23*s38*s44)*s58 + (-0.25*(s23*s31*s36*s43) + (s31*s38*s44)/2. + (s23*s41*s44)/2.)*s59 + (-0.25*(s12*s31*s36*s43) + (s12*s41*s44)/2.)*s60 + (-0.25*(s12*s23*s36*s43) + (s12*s38*s44)/2.)*s61 - (s13*s36*s43*s62)/4. + (s43*s44*s63*s64)/12. + s54*((s43*s44*s55)/12. + s65);
				buffer__[9*i+1] += -0.25*(s32*s36*s45*s53) + (s44*s46*s53)/4. + (s44*s53*s54*s55)/12. + (s44*s53*s56*s57)/12. + (s23*s44*s48 - (s24*s36*s53)/4.)*s58 + ((s31*s44*s48)/2. - (s23*s31*s36*s53)/4.)*s59 + ((s31*s44*s51)/2. - (s12*s31*s36*s53)/4.)*s60 + ((s12*s44*s48)/2. + (s23*s44*s51)/2. - (s12*s23*s36*s53)/4.)*s61 + (s12*s44*s51 - (s13*s36*s53)/4.)*s62 + s63*((s44*s53*s64)/12. + s65);
				buffer__[9*i+2] += (s44*s46*s72)/4. + (s44*s54*s55*s72)/12. - (s24*s36*s58*s72)/4. + (s44*s63*s64*s72)/12. + s62*(s12*s44*s70 - (s13*s36*s72)/4.) + s61*((s23*s44*s70)/2. - (s12*s23*s36*s72)/4.) + s60*((s12*s44*s67)/2. + (s31*s44*s70)/2. - (s12*s31*s36*s72)/4.) + s59*((s23*s44*s67)/2. - (s23*s31*s36*s72)/4.) + s45*(s31*s44*s67 - (s32*s36*s72)/4.) + s56*(s65 + (s44*s57*s72)/12.);
				buffer__[9*i+3] += (s44*s46*s78)/4. + (s44*s56*s57*s78)/12. - (s13*s36*s62*s78)/4. + (s44*s63*s64*s78)/12. + s61*((s12*s44*s74)/2. - (s12*s23*s36*s78)/4.) + s58*(s23*s44*s74 - (s24*s36*s78)/4.) + s60*((s12*s44*s76)/2. - (s12*s31*s36*s78)/4.) + s59*((s31*s44*s74)/2. + (s23*s44*s76)/2. - (s23*s31*s36*s78)/4.) + s45*(s31*s44*s76 - (s32*s36*s78)/4.) + s54*(s65 + (s44*s55*s78)/12.);
				buffer__[9*i+4] += -0.25*(s32*s36*s45*s84) + (s44*s46*s84)/4. + (s44*s54*s55*s84)/12. + (s44*s56*s57*s84)/12. + s62*(s12*s44*s82 - (s13*s36*s84)/4.) + s61*((s12*s44*s79)/2. + (s23*s44*s82)/2. - (s12*s23*s36*s84)/4.) + s58*(s23*s44*s79 - (s24*s36*s84)/4.) + s60*((s31*s44*s82)/2. - (s12*s31*s36*s84)/4.) + s59*((s31*s44*s79)/2. - (s23*s31*s36*s84)/4.) + s63*(s65 + (s44*s64*s84)/12.);
				buffer__[9*i+5] += (s44*s46*s90)/4. + (s44*s54*s55*s90)/12. - (s24*s36*s58*s90)/4. + (s44*s63*s64*s90)/12. + s62*(s12*s44*s88 - (s13*s36*s90)/4.) + s61*((s23*s44*s88)/2. - (s12*s23*s36*s90)/4.) + s60*((s12*s44*s86)/2. + (s31*s44*s88)/2. - (s12*s31*s36*s90)/4.) + s59*((s23*s44*s86)/2. - (s23*s31*s36*s90)/4.) + s45*(s31*s44*s86 - (s32*s36*s90)/4.) + s56*(s65 + (s44*s57*s90)/12.);
				buffer__[9*i+6] += (s44*s46*s95)/4. + (s44*s56*s57*s95)/12. - (s13*s36*s62*s95)/4. + (s44*s63*s64*s95)/12. + s61*((s12*s44*s91)/2. - (s12*s23*s36*s95)/4.) + s58*(s23*s44*s91 - (s24*s36*s95)/4.) + s60*((s12*s44*s93)/2. - (s12*s31*s36*s95)/4.) + s59*((s31*s44*s91)/2. + (s23*s44*s93)/2. - (s23*s31*s36*s95)/4.) + s45*(s31*s44*s93 - (s32*s36*s95)/4.) + s54*(s65 + (s44*s55*s95)/12.);
				buffer__[9*i+7] += -0.25*(s100*s32*s36*s45) + (s100*s44*s46)/4. + (s100*s44*s54*s55)/12. + (s100*s44*s56*s57)/12. + s63*((s100*s44*s64)/12. + s65) + s58*(-0.25*(s100*s24*s36) + s23*s44*s96) + s59*(-0.25*(s100*s23*s31*s36) + (s31*s44*s96)/2.) + s62*(-0.25*(s100*s13*s36) + s12*s44*s98) + s61*(-0.25*(s100*s12*s23*s36) + (s12*s44*s96)/2. + (s23*s44*s98)/2.) + s60*(-0.25*(s100*s12*s31*s36) + (s31*s44*s98)/2.);
				buffer__[9*i+8] += (-0.25*(s105*s32*s36) + s101*s31*s44)*s45 + (s105*s44*s46)/4. + (s105*s44*s54*s55)/12. - (s105*s24*s36*s58)/4. + (-0.25*(s105*s23*s31*s36) + (s101*s23*s44)/2.)*s59 + (-0.25*(s105*s12*s31*s36) + (s101*s12*s44)/2. + (s103*s31*s44)/2.)*s60 + (-0.25*(s105*s12*s23*s36) + (s103*s23*s44)/2.)*s61 + (-0.25*(s105*s13*s36) + s103*s12*s44)*s62 + (s105*s44*s63*s64)/12. + s56*((s105*s44*s57)/12. + s65);
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[3*simplices__[3*i+0]+2];
				const Real s1 = V_coords__[3*simplices__[3*i+1]+1];
				const Real s2 = -(s0*s1);
				const Real s3 = V_coords__[3*simplices__[3*i+0]+1];
				const Real s4 = V_coords__[3*simplices__[3*i+1]+2];
				const Real s5 = s3*s4;
				const Real s6 = V_coords__[3*simplices__[3*i+2]+1];
				const Real s7 = s0*s6;
				const Real s8 = -(s4*s6);
				const Real s9 = V_coords__[3*simplices__[3*i+2]+2];
				const Real s10 = -(s3*s9);
				const Real s11 = s1*s9;
				const Real s12 = s10 + s11 + s2 + s5 + s7 + s8;
				const Real s13 = s12*s12;
				const Real s14 = V_coords__[3*simplices__[3*i+2]+0];
				const Real s15 = V_coords__[3*simplices__[3*i+0]+0];
				const Real s16 = V_coords__[3*simplices__[3*i+1]+0];
				const Real s17 = -(s16*s3);
				const Real s18 = s1*s15;
				const Real s19 = s14*s3;
				const Real s20 = -(s1*s14);
				const Real s21 = -(s15*s6);
				const Real s22 = s16*s6;
				const Real s23 = s17 + s18 + s19 + s20 + s21 + s22;
				const Real s24 = s23*s23;
				const Real s25 = s0*s16;
				const Real s26 = -(s15*s4);
				const Real s27 = -(s0*s14);
				const Real s28 = s14*s4;
				const Real s29 = s15*s9;
				const Real s30 = -(s16*s9);
				const Real s31 = s25 + s26 + s27 + s28 + s29 + s30;
				const Real s32 = s31*s31;
				const Real s33 = s13 + s24 + s32;
				const Real s34 = sqrt(s33);
				const Real s35 = s33*s34;
				const Real s36 = 1/s35;
				const Real s37 = -s6;
				const Real s38 = s1 + s37;
				const Real s39 = 2*s23*s38;
				const Real s40 = -s4;
				const Real s41 = s40 + s9;
				const Real s42 = 2*s31*s41;
				const Real s43 = s39 + s42;
				const Real s44 = 1/s34;
				const Real s45 = P_D_far__[10*i+7];
				const Real s46 = P_D_far__[10*i+0];
				const Real s47 = -s16;
				const Real s48 = s14 + s47;
				const Real s49 = 2*s23*s48;
				const Real s50 = -s9;
				const Real s51 = s4 + s50;
				const Real s52 = 2*s12*s51;
				const Real s53 = s49 + s52;
				const Real s54 = P_D_far__[10*i+1];
				const Real s55 = s14 + s15 + s16;
				const Real s56 = P_D_far__[10*i+3];
				const Real s57 = s0 + s4 + s9;
				const Real s58 = P_D_far__[10*i+9];
				const Real s59 = P_D_far__[10*i+8];
				const Real s60 = P_D_far__[10*i+5];
				const Real s61 = P_D_far__[10*i+6];
				const Real s62 = P_D_far__[10*i+4];
				const Real s63 = P_D_far__[10*i+2];
				const Real s64 = s1 + s3 + s6;
				const Real s65 = s34/6.;
				const Real s66 = -s14;
				const Real s67 = s16 + s66;
				const Real s68 = 2*s31*s67;
				const Real s69 = -s1;
				const Real s70 = s6 + s69;
				const Real s71 = 2*s12*s70;
				const Real s72 = s68 + s71;
				const Real s73 = -s3;
				const Real s74 = s6 + s73;
				const Real s75 = 2*s23*s74;
				const Real s76 = s0 + s50;
				const Real s77 = 2*s31*s76;
				const Real s78 = s75 + s77;
				const Real s79 = s15 + s66;
				const Real s80 = 2*s23*s79;
				const Real s81 = -s0;
				const Real s82 = s81 + s9;
				const Real s83 = 2*s12*s82;
				const Real s84 = s80 + s83;
				const Real s85 = -s15;
				const Real s86 = s14 + s85;
				const Real s87 = 2*s31*s86;
				const Real s88 = s3 + s37;
				const Real s89 = 2*s12*s88;
				const Real s90 = s87 + s89;
				const Real s91 = s3 + s69;
				const Real s92 = 2*s23*s91;
				const Real s93 = s4 + s81;
				const Real s94 = 2*s31*s93;
				const Real s95 = s92 + s94;
				const Real s96 = s16 + s85;
				const Real s97 = 2*s23*s96;
				const Real s98 = s0 + s40;
				const Real s99 = 2*s12*s98;
				const Real s100 = s97 + s99;
				const Real s101 = s15 + s47;
				const Real s102 = 2*s101*s31;
				const Real s103 = s1 + s73;
				const Real s104 = 2*s103*s12;
				const Real s105 = s102 + s104;
				buffer__[9*i+0] = (-0.25*(s32*s36*s43) + s31*s41*s44)*s45 + (s43*s44*s46)/4. + (s43*s44*s56*s57)/12. + (-0.25*(s24*s36*s43) + s23*s38*s44)*s58 + (-0.25*(s23*s31*s36*s43) + (s31*s38*s44)/2. + (s23*s41*s44)/2.)*s59 + (-0.25*(s12*s31*s36*s43) + (s12*s41*s44)/2.)*s60 + (-0.25*(s12*s23*s36*s43) + (s12*s38*s44)/2.)*s61 - (s13*s36*s43*s62)/4. + (s43*s44*s63*s64)/12. + s54*((s43*s44*s55)/12. + s65);
				buffer__[9*i+1] = -0.25*(s32*s36*s45*s53) + (s44*s46*s53)/4. + (s44*s53*s54*s55)/12. + (s44*s53*s56*s57)/12. + (s23*s44*s48 - (s24*s36*s53)/4.)*s58 + ((s31*s44*s48)/2. - (s23*s31*s36*s53)/4.)*s59 + ((s31*s44*s51)/2. - (s12*s31*s36*s53)/4.)*s60 + ((s12*s44*s48)/2. + (s23*s44*s51)/2. - (s12*s23*s36*s53)/4.)*s61 + (s12*s44*s51 - (s13*s36*s53)/4.)*s62 + s63*((s44*s53*s64)/12. + s65);
				buffer__[9*i+2] = (s44*s46*s72)/4. + (s44*s54*s55*s72)/12. - (s24*s36*s58*s72)/4. + (s44*s63*s64*s72)/12. + s62*(s12*s44*s70 - (s13*s36*s72)/4.) + s61*((s23*s44*s70)/2. - (s12*s23*s36*s72)/4.) + s60*((s12*s44*s67)/2. + (s31*s44*s70)/2. - (s12*s31*s36*s72)/4.) + s59*((s23*s44*s67)/2. - (s23*s31*s36*s72)/4.) + s45*(s31*s44*s67 - (s32*s36*s72)/4.) + s56*(s65 + (s44*s57*s72)/12.);
				buffer__[9*i+3] = (s44*s46*s78)/4. + (s44*s56*s57*s78)/12. - (s13*s36*s62*s78)/4. + (s44*s63*s64*s78)/12. + s61*((s12*s44*s74)/2. - (s12*s23*s36*s78)/4.) + s58*(s23*s44*s74 - (s24*s36*s78)/4.) + s60*((s12*s44*s76)/2. - (s12*s31*s36*s78)/4.) + s59*((s31*s44*s74)/2. + (s23*s44*s76)/2. - (s23*s31*s36*s78)/4.) + s45*(s31*s44*s76 - (s32*s36*s78)/4.) + s54*(s65 + (s44*s55*s78)/12.);
				buffer__[9*i+4] = -0.25*(s32*s36*s45*s84) + (s44*s46*s84)/4. + (s44*s54*s55*s84)/12. + (s44*s56*s57*s84)/12. + s62*(s12*s44*s82 - (s13*s36*s84)/4.) + s61*((s12*s44*s79)/2. + (s23*s44*s82)/2. - (s12*s23*s36*s84)/4.) + s58*(s23*s44*s79 - (s24*s36*s84)/4.) + s60*((s31*s44*s82)/2. - (s12*s31*s36*s84)/4.) + s59*((s31*s44*s79)/2. - (s23*s31*s36*s84)/4.) + s63*(s65 + (s44*s64*s84)/12.);
				buffer__[9*i+5] = (s44*s46*s90)/4. + (s44*s54*s55*s90)/12. - (s24*s36*s58*s90)/4. + (s44*s63*s64*s90)/12. + s62*(s12*s44*s88 - (s13*s36*s90)/4.) + s61*((s23*s44*s88)/2. - (s12*s23*s36*s90)/4.) + s60*((s12*s44*s86)/2. + (s31*s44*s88)/2. - (s12*s31*s36*s90)/4.) + s59*((s23*s44*s86)/2. - (s23*s31*s36*s90)/4.) + s45*(s31*s44*s86 - (s32*s36*s90)/4.) + s56*(s65 + (s44*s57*s90)/12.);
				buffer__[9*i+6] = (s44*s46*s95)/4. + (s44*s56*s57*s95)/12. - (s13*s36*s62*s95)/4. + (s44*s63*s64*s95)/12. + s61*((s12*s44*s91)/2. - (s12*s23*s36*s95)/4.) + s58*(s23*s44*s91 - (s24*s36*s95)/4.) + s60*((s12*s44*s93)/2. - (s12*s31*s36*s95)/4.) + s59*((s31*s44*s91)/2. + (s23*s44*s93)/2. - (s23*s31*s36*s95)/4.) + s45*(s31*s44*s93 - (s32*s36*s95)/4.) + s54*(s65 + (s44*s55*s95)/12.);
				buffer__[9*i+7] = -0.25*(s100*s32*s36*s45) + (s100*s44*s46)/4. + (s100*s44*s54*s55)/12. + (s100*s44*s56*s57)/12. + s63*((s100*s44*s64)/12. + s65) + s58*(-0.25*(s100*s24*s36) + s23*s44*s96) + s59*(-0.25*(s100*s23*s31*s36) + (s31*s44*s96)/2.) + s62*(-0.25*(s100*s13*s36) + s12*s44*s98) + s61*(-0.25*(s100*s12*s23*s36) + (s12*s44*s96)/2. + (s23*s44*s98)/2.) + s60*(-0.25*(s100*s12*s31*s36) + (s31*s44*s98)/2.);
				buffer__[9*i+8] = (-0.25*(s105*s32*s36) + s101*s31*s44)*s45 + (s105*s44*s46)/4. + (s105*s44*s54*s55)/12. - (s105*s24*s36*s58)/4. + (-0.25*(s105*s23*s31*s36) + (s101*s23*s44)/2.)*s59 + (-0.25*(s105*s12*s31*s36) + (s101*s12*s44)/2. + (s103*s31*s44)/2.)*s60 + (-0.25*(s105*s12*s23*s36) + (s103*s23*s44)/2.)*s61 + (-0.25*(s105*s13*s36) + s103*s12*s44)*s62 + (s105*s44*s63*s64)/12. + s56*((s105*s44*s57)/12. + s65);
			}
		}

        ptoc(ClassName()+"::DFarToHulls");
        
    }

	}; // SimplicialMeshDetails<2,3,Real,Int>

//----------------------------------------------------------------------------------------------

    template<typename Real, typename Int>
    struct SimplicialMeshDetails<2,4,Real,Int>
    {
	private:

		const Int thread_count = 1;

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(std::max(static_cast<Int>(1),thread_count_))
		{}
	
		inline Int ThreadCount() const
		{
			return thread_count;
		}
	
		inline std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<2,4,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void ComputeNearFarData( 
		const Tensor2<Real,Int> & V_coords,
		const Tensor2<Int ,Int> & simplices,
			  Tensor2<Real,Int> & P_near,
			  Tensor2<Real,Int> & P_far
	) const
    {
        ptic(ClassName()+"::ComputeNearFarData");
        
        //Int size       = 3;
        //Int amb_dim    = 4;
        //Int dom_dim    = 2;

        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			ptr<Real> V_coords__      = V_coords.data();	
			ptr<Int>  simplices__     = simplices.data();

			Real hull    [3][4];
			Real df      [4][2];
			Real dfdagger[2][4];
			Real g       [2][2];
			Real ginv    [2][2];

			Int simplex  [3];
			
			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {
				Real * restrict const near = P_near.data(i);                    
				Real * restrict const far  = P_far.data(i);   
            
				simplex[0] = simplices__[3*i +0];
				simplex[1] = simplices__[3*i +1];
				simplex[2] = simplices__[3*i +2];

				near[1] = hull[0][0] = V_coords__[4*simplex[0]+0];
				near[2] = hull[0][1] = V_coords__[4*simplex[0]+1];
				near[3] = hull[0][2] = V_coords__[4*simplex[0]+2];
				near[4] = hull[0][3] = V_coords__[4*simplex[0]+3];
				near[5] = hull[1][0] = V_coords__[4*simplex[1]+0];
				near[6] = hull[1][1] = V_coords__[4*simplex[1]+1];
				near[7] = hull[1][2] = V_coords__[4*simplex[1]+2];
				near[8] = hull[1][3] = V_coords__[4*simplex[1]+3];
				near[9] = hull[2][0] = V_coords__[4*simplex[2]+0];
				near[10] = hull[2][1] = V_coords__[4*simplex[2]+1];
				near[11] = hull[2][2] = V_coords__[4*simplex[2]+2];
				near[12] = hull[2][3] = V_coords__[4*simplex[2]+3];

				far[1] = static_cast<Real>(0.3333333333333333) * ( hull[0][0] + hull[1][0] + hull[2][0] );
				far[2] = static_cast<Real>(0.3333333333333333) * ( hull[0][1] + hull[1][1] + hull[2][1] );
				far[3] = static_cast<Real>(0.3333333333333333) * ( hull[0][2] + hull[1][2] + hull[2][2] );
				far[4] = static_cast<Real>(0.3333333333333333) * ( hull[0][3] + hull[1][3] + hull[2][3] );

				df[0][0] = hull[1][0] - hull[0][0];
				df[0][1] = hull[2][0] - hull[0][0];
				df[1][0] = hull[1][1] - hull[0][1];
				df[1][1] = hull[2][1] - hull[0][1];
				df[2][0] = hull[1][2] - hull[0][2];
				df[2][1] = hull[2][2] - hull[0][2];
				df[3][0] = hull[1][3] - hull[0][3];
				df[3][1] = hull[2][3] - hull[0][3];

				g[0][0] = df[0][0] * df[0][0] + df[1][0] * df[1][0] + df[2][0] * df[2][0] + df[3][0] * df[3][0];
				g[0][1] = df[0][0] * df[0][1] + df[1][0] * df[1][1] + df[2][0] * df[2][1] + df[3][0] * df[3][1];
				g[1][0] = df[0][1] * df[0][0] + df[1][1] * df[1][0] + df[2][1] * df[2][0] + df[3][1] * df[3][0];
				g[1][1] = df[0][1] * df[0][1] + df[1][1] * df[1][1] + df[2][1] * df[2][1] + df[3][1] * df[3][1];

                Real det = g[0][0] * g[1][1] - g[0][1] * g[0][1];

				near[0] = far[0] = sqrt( fabs(det) ) * static_cast<Real>(0.5);

                Real invdet = static_cast<Real>(1)/det;
                ginv[0][0] =  g[1][1] * invdet;
                ginv[0][1] = -g[0][1] * invdet;
                ginv[1][1] =  g[0][0] * invdet;
                
                //  dfdagger = g^{-1} * df^T (2 x 4 matrix)
				dfdagger[0][0] = ginv[0][0] * df[0][0] + ginv[0][1] * df[0][1];
				dfdagger[0][1] = ginv[0][0] * df[1][0] + ginv[0][1] * df[1][1];
				dfdagger[0][2] = ginv[0][0] * df[2][0] + ginv[0][1] * df[2][1];
				dfdagger[0][3] = ginv[0][0] * df[3][0] + ginv[0][1] * df[3][1];
				dfdagger[1][0] = ginv[0][1] * df[0][0] + ginv[1][1] * df[0][1];
				dfdagger[1][1] = ginv[0][1] * df[1][0] + ginv[1][1] * df[1][1];
				dfdagger[1][2] = ginv[0][1] * df[2][0] + ginv[1][1] * df[2][1];
				dfdagger[1][3] = ginv[0][1] * df[3][0] + ginv[1][1] * df[3][1];
            
				near[13] = far[ 5]  = static_cast<Real>(1) - df[0][0] * dfdagger[0][0] - df[0][1] * dfdagger[1][0];
				near[14] = far[ 6]  =    - df[0][0] * dfdagger[0][1] - df[0][1] * dfdagger[1][1];
				near[15] = far[ 7]  =    - df[0][0] * dfdagger[0][2] - df[0][1] * dfdagger[1][2];
				near[16] = far[ 8]  =    - df[0][0] * dfdagger[0][3] - df[0][1] * dfdagger[1][3];
				near[17] = far[ 9]  = static_cast<Real>(1) - df[1][0] * dfdagger[0][1] - df[1][1] * dfdagger[1][1];
				near[18] = far[10]  =    - df[1][0] * dfdagger[0][2] - df[1][1] * dfdagger[1][2];
				near[19] = far[11]  =    - df[1][0] * dfdagger[0][3] - df[1][1] * dfdagger[1][3];
				near[20] = far[12]  = static_cast<Real>(1) - df[2][0] * dfdagger[0][2] - df[2][1] * dfdagger[1][2];
				near[21] = far[13]  =    - df[2][0] * dfdagger[0][3] - df[2][1] * dfdagger[1][3];
				near[22] = far[14]  = static_cast<Real>(1) - df[3][0] * dfdagger[0][3] - df[3][1] * dfdagger[1][3];

            } // for( Int i = i_begin; i < i_end; ++i )

        } // #pragma omp parallel num_threads( ThreadCount() )

        ptoc(ClassName()+"::ComputeNearFarData");
    }

	void ComputeNearFarDataOps(
		const Tensor2<Real,Int> & V_coords,
        const Tensor2<Int ,Int> & simplices,
		      Tensor2<Real,Int> & P_coords,
		      Tensor3<Real,Int> & P_hull_coords,
		      Tensor2<Real,Int> & P_near,
		      Tensor2<Real,Int> & P_far,
		Sparse::MatrixCSR<Real,Int,Int> & DiffOp,
		Sparse::MatrixCSR<Real,Int,Int> & AvOp 
	) const
    {
        ptic(ClassName()+"::ComputeNearFarDataOps");

        //Int size       = 3;
        //Int amb_dim    = 4;
        //Int dom_dim    = 2;

        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );

        #pragma omp parallel for num_threads( ThreadCount() )
        for( Int thread = 0; thread < ThreadCount(); ++thread )
        {
			mut<Int>  AvOp_outer = AvOp.Outer().data();
			mut<Int>  AvOp_inner = AvOp.Inner().data();
			mut<Real> AvOp_value = AvOp.Values().data();

			mut<Int>  DiffOp_outer = DiffOp.Outer().data();
			mut<Int>  DiffOp_inner = DiffOp.Inner().data();
			mut<Real> DiffOp_value = DiffOp.Values().data();

			ptr<Real> V_coords__      = V_coords.data();
			
			ptr<Int>  simplices__     = simplices.data();
		    mut<Real> P_hull_coords__ = P_hull_coords.data();
			mut<Real> P_coords__      = P_coords.data();

			Real df       [4][2];
			Real dfdagger [2][4];
			Real g        [2][2];
			Real ginv     [2][2];

			Int simplex        [3];
            Int sorted_simplex [3];

			const Int i_begin = job_ptr[thread];
			const Int i_end   = job_ptr[thread+1];

            for( Int i = i_begin; i < i_end; ++i )
            {

				mut<Real> near = P_near.data(i);                    
				mut<Real> far  = P_far.data(i);

				simplex[0] = sorted_simplex[0] = simplices__[3*i +0];
				simplex[1] = sorted_simplex[1] = simplices__[3*i +1];
				simplex[2] = sorted_simplex[2] = simplices__[3*i +2];
                  
                // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
                std::sort( sorted_simplex, sorted_simplex + 3 );

				AvOp_outer[i+1] = (i+1) * 3;                      
				AvOp_inner[3*i+0] = sorted_simplex[0];
				AvOp_inner[3*i+1] = sorted_simplex[1];
				AvOp_inner[3*i+2] = sorted_simplex[2];

				AvOp_value[3*i+0] = 0.3333333333333333;
				AvOp_value[3*i+1] = 0.3333333333333333;
				AvOp_value[3*i+2] = 0.3333333333333333;

				DiffOp_outer[4*i+0] = (4 * i + 0) * 3;
				DiffOp_outer[4*i+1] = (4 * i + 1) * 3;
				DiffOp_outer[4*i+2] = (4 * i + 2) * 3;
				DiffOp_outer[4*i+3] = (4 * i + 3) * 3;

				DiffOp_inner[(i * 4 + 0) * 3 + 0 ] = sorted_simplex[0];
				DiffOp_inner[(i * 4 + 0) * 3 + 1 ] = sorted_simplex[1];
				DiffOp_inner[(i * 4 + 0) * 3 + 2 ] = sorted_simplex[2];
				DiffOp_inner[(i * 4 + 1) * 3 + 0 ] = sorted_simplex[0];
				DiffOp_inner[(i * 4 + 1) * 3 + 1 ] = sorted_simplex[1];
				DiffOp_inner[(i * 4 + 1) * 3 + 2 ] = sorted_simplex[2];
				DiffOp_inner[(i * 4 + 2) * 3 + 0 ] = sorted_simplex[0];
				DiffOp_inner[(i * 4 + 2) * 3 + 1 ] = sorted_simplex[1];
				DiffOp_inner[(i * 4 + 2) * 3 + 2 ] = sorted_simplex[2];
				DiffOp_inner[(i * 4 + 3) * 3 + 0 ] = sorted_simplex[0];
				DiffOp_inner[(i * 4 + 3) * 3 + 1 ] = sorted_simplex[1];
				DiffOp_inner[(i * 4 + 3) * 3 + 2 ] = sorted_simplex[2];

				near[1] = P_hull_coords__[12*i+0] = V_coords__[4*simplex[0]+0];
				near[2] = P_hull_coords__[12*i+1] = V_coords__[4*simplex[0]+1];
				near[3] = P_hull_coords__[12*i+2] = V_coords__[4*simplex[0]+2];
				near[4] = P_hull_coords__[12*i+3] = V_coords__[4*simplex[0]+3];
				near[5] = P_hull_coords__[12*i+4] = V_coords__[4*simplex[1]+0];
				near[6] = P_hull_coords__[12*i+5] = V_coords__[4*simplex[1]+1];
				near[7] = P_hull_coords__[12*i+6] = V_coords__[4*simplex[1]+2];
				near[8] = P_hull_coords__[12*i+7] = V_coords__[4*simplex[1]+3];
				near[9] = P_hull_coords__[12*i+8] = V_coords__[4*simplex[2]+0];
				near[10] = P_hull_coords__[12*i+9] = V_coords__[4*simplex[2]+1];
				near[11] = P_hull_coords__[12*i+10] = V_coords__[4*simplex[2]+2];
				near[12] = P_hull_coords__[12*i+11] = V_coords__[4*simplex[2]+3];

				far[1] = P_coords__[4*i+0] = 0.3333333333333333 * ( P_hull_coords__[12*i+0] + P_hull_coords__[12*i+4] + P_hull_coords__[12*i+8] );
				far[2] = P_coords__[4*i+1] = 0.3333333333333333 * ( P_hull_coords__[12*i+1] + P_hull_coords__[12*i+5] + P_hull_coords__[12*i+9] );
				far[3] = P_coords__[4*i+2] = 0.3333333333333333 * ( P_hull_coords__[12*i+2] + P_hull_coords__[12*i+6] + P_hull_coords__[12*i+10] );
				far[4] = P_coords__[4*i+3] = 0.3333333333333333 * ( P_hull_coords__[12*i+3] + P_hull_coords__[12*i+7] + P_hull_coords__[12*i+11] );

				df[0][0] = V_coords__[4*sorted_simplex[1]+0] - V_coords__[4*sorted_simplex[0]+0];
				df[0][1] = V_coords__[4*sorted_simplex[2]+0] - V_coords__[4*sorted_simplex[0]+0];
				df[1][0] = V_coords__[4*sorted_simplex[1]+1] - V_coords__[4*sorted_simplex[0]+1];
				df[1][1] = V_coords__[4*sorted_simplex[2]+1] - V_coords__[4*sorted_simplex[0]+1];
				df[2][0] = V_coords__[4*sorted_simplex[1]+2] - V_coords__[4*sorted_simplex[0]+2];
				df[2][1] = V_coords__[4*sorted_simplex[2]+2] - V_coords__[4*sorted_simplex[0]+2];
				df[3][0] = V_coords__[4*sorted_simplex[1]+3] - V_coords__[4*sorted_simplex[0]+3];
				df[3][1] = V_coords__[4*sorted_simplex[2]+3] - V_coords__[4*sorted_simplex[0]+3];

				g[0][0] = df[0][0] * df[0][0] + df[1][0] * df[1][0] + df[2][0] * df[2][0] + df[3][0] * df[3][0];
				g[0][1] = df[0][0] * df[0][1] + df[1][0] * df[1][1] + df[2][0] * df[2][1] + df[3][0] * df[3][1];
				g[1][0] = df[0][1] * df[0][0] + df[1][1] * df[1][0] + df[2][1] * df[2][0] + df[3][1] * df[3][0];
				g[1][1] = df[0][1] * df[0][1] + df[1][1] * df[1][1] + df[2][1] * df[2][1] + df[3][1] * df[3][1];

                Real det = g[0][0] * g[1][1] - g[0][1] * g[0][1];

                near[0] = far[0] = sqrt( fabs(det) ) * static_cast<Real>(0.5);

                Real invdet = static_cast<Real>(1.0)/det;
                ginv[0][0] =  g[1][1] * invdet;
                ginv[0][1] = -g[0][1] * invdet;
                ginv[1][1] =  g[0][0] * invdet;
                
                //  dfdagger = g^{-1} * df^T (2 x 4 matrix)
				dfdagger[0][0] = ginv[0][0] * df[0][0] + ginv[0][1] * df[0][1];
				dfdagger[0][1] = ginv[0][0] * df[1][0] + ginv[0][1] * df[1][1];
				dfdagger[0][2] = ginv[0][0] * df[2][0] + ginv[0][1] * df[2][1];
				dfdagger[0][3] = ginv[0][0] * df[3][0] + ginv[0][1] * df[3][1];
				dfdagger[1][0] = ginv[0][1] * df[0][0] + ginv[1][1] * df[0][1];
				dfdagger[1][1] = ginv[0][1] * df[1][0] + ginv[1][1] * df[1][1];
				dfdagger[1][2] = ginv[0][1] * df[2][0] + ginv[1][1] * df[2][1];
				dfdagger[1][3] = ginv[0][1] * df[3][0] + ginv[1][1] * df[3][1];
            
				near[13] = far[ 5]  = static_cast<Real>(1.0) - df[0][0] * dfdagger[0][0] - df[0][1] * dfdagger[1][0];
				near[14] = far[ 6]  =    - df[0][0] * dfdagger[0][1] - df[0][1] * dfdagger[1][1];
				near[15] = far[ 7]  =    - df[0][0] * dfdagger[0][2] - df[0][1] * dfdagger[1][2];
				near[16] = far[ 8]  =    - df[0][0] * dfdagger[0][3] - df[0][1] * dfdagger[1][3];
				near[17] = far[ 9]  = static_cast<Real>(1.0) - df[1][0] * dfdagger[0][1] - df[1][1] * dfdagger[1][1];
				near[18] = far[10]  =    - df[1][0] * dfdagger[0][2] - df[1][1] * dfdagger[1][2];
				near[19] = far[11]  =    - df[1][0] * dfdagger[0][3] - df[1][1] * dfdagger[1][3];
				near[20] = far[12]  = static_cast<Real>(1.0) - df[2][0] * dfdagger[0][2] - df[2][1] * dfdagger[1][2];
				near[21] = far[13]  =    - df[2][0] * dfdagger[0][3] - df[2][1] * dfdagger[1][3];
				near[22] = far[14]  = static_cast<Real>(1.0) - df[3][0] * dfdagger[0][3] - df[3][1] * dfdagger[1][3];

                // derivative operator  (4 x 3 matrix)

                mut<Real> Df = &DiffOp_value[ 12 * i ];

				Df[ 0] = - dfdagger[0][0] - dfdagger[1][0];
				Df[ 1] =   dfdagger[0][0];
				Df[ 2] =   dfdagger[1][0];
				Df[ 3] = - dfdagger[0][1] - dfdagger[1][1];
				Df[ 4] =   dfdagger[0][1];
				Df[ 5] =   dfdagger[1][1];
				Df[ 6] = - dfdagger[0][2] - dfdagger[1][2];
				Df[ 7] =   dfdagger[0][2];
				Df[ 8] =   dfdagger[1][2];
				Df[ 9] = - dfdagger[0][3] - dfdagger[1][3];
				Df[10] =   dfdagger[0][3];
				Df[11] =   dfdagger[1][3];

            }
        }
        ptoc(ClassName()+"::ComputeNearFarDataOps");
    }

    void DNearToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_near, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DNearToHulls");

        if( P_D_near.Dimension(1) != 23 )
        {
            eprint("in DNearToHulls: P_D_near.Dimension(1) != 23. Aborting");
        }
        
		ptr<Real> V_coords__  = V_coords.data();
        ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_near__  = P_D_near.data();
        mut<Real> buffer__    = buffer.data();

		if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[4*simplices__[3*i+1]+1];
				const Real s1 = s0*s0;
				const Real s2 = V_coords__[4*simplices__[3*i+1]+0];
				const Real s3 = V_coords__[4*simplices__[3*i+0]+0];
				const Real s4 = V_coords__[4*simplices__[3*i+1]+2];
				const Real s5 = s4*s4;
				const Real s6 = V_coords__[4*simplices__[3*i+1]+3];
				const Real s7 = s6*s6;
				const Real s8 = V_coords__[4*simplices__[3*i+0]+1];
				const Real s9 = V_coords__[4*simplices__[3*i+2]+0];
				const Real s10 = V_coords__[4*simplices__[3*i+0]+2];
				const Real s11 = V_coords__[4*simplices__[3*i+0]+3];
				const Real s12 = V_coords__[4*simplices__[3*i+2]+1];
				const Real s13 = s12*s12;
				const Real s14 = V_coords__[4*simplices__[3*i+2]+2];
				const Real s15 = s14*s14;
				const Real s16 = V_coords__[4*simplices__[3*i+2]+3];
				const Real s17 = s16*s16;
				const Real s18 = s8*s8;
				const Real s19 = s2*s2;
				const Real s20 = s18*s19;
				const Real s21 = s10*s10;
				const Real s22 = s19*s21;
				const Real s23 = s11*s11;
				const Real s24 = s19*s23;
				const Real s25 = -2*s0*s2*s3*s8;
				const Real s26 = s3*s3;
				const Real s27 = s1*s26;
				const Real s28 = s1*s21;
				const Real s29 = s1*s23;
				const Real s30 = -2*s10*s2*s3*s4;
				const Real s31 = -2*s0*s10*s4*s8;
				const Real s32 = s26*s5;
				const Real s33 = s18*s5;
				const Real s34 = s23*s5;
				const Real s35 = -2*s11*s2*s3*s6;
				const Real s36 = -2*s0*s11*s6*s8;
				const Real s37 = -2*s10*s11*s4*s6;
				const Real s38 = s26*s7;
				const Real s39 = s18*s7;
				const Real s40 = s21*s7;
				const Real s41 = -2*s18*s2*s9;
				const Real s42 = -2*s2*s21*s9;
				const Real s43 = -2*s2*s23*s9;
				const Real s44 = 2*s0*s3*s8*s9;
				const Real s45 = 2*s0*s2*s8*s9;
				const Real s46 = -2*s1*s3*s9;
				const Real s47 = 2*s10*s3*s4*s9;
				const Real s48 = 2*s10*s2*s4*s9;
				const Real s49 = -2*s3*s5*s9;
				const Real s50 = 2*s11*s3*s6*s9;
				const Real s51 = 2*s11*s2*s6*s9;
				const Real s52 = -2*s3*s7*s9;
				const Real s53 = s9*s9;
				const Real s54 = s18*s53;
				const Real s55 = s21*s53;
				const Real s56 = s23*s53;
				const Real s57 = -2*s0*s53*s8;
				const Real s58 = s1*s53;
				const Real s59 = -2*s10*s4*s53;
				const Real s60 = s5*s53;
				const Real s61 = -2*s11*s53*s6;
				const Real s62 = s53*s7;
				const Real s63 = 2*s12*s2*s3*s8;
				const Real s64 = -2*s12*s19*s8;
				const Real s65 = -2*s0*s12*s26;
				const Real s66 = -2*s0*s12*s21;
				const Real s67 = -2*s0*s12*s23;
				const Real s68 = 2*s0*s12*s2*s3;
				const Real s69 = 2*s10*s12*s4*s8;
				const Real s70 = 2*s0*s10*s12*s4;
				const Real s71 = -2*s12*s5*s8;
				const Real s72 = 2*s11*s12*s6*s8;
				const Real s73 = 2*s0*s11*s12*s6;
				const Real s74 = -2*s12*s7*s8;
				const Real s75 = -2*s12*s3*s8*s9;
				const Real s76 = 2*s12*s2*s8*s9;
				const Real s77 = 2*s0*s12*s3*s9;
				const Real s78 = -2*s0*s12*s2*s9;
				const Real s79 = s13*s26;
				const Real s80 = s13*s21;
				const Real s81 = s13*s23;
				const Real s82 = -2*s13*s2*s3;
				const Real s83 = s13*s19;
				const Real s84 = -2*s10*s13*s4;
				const Real s85 = s13*s5;
				const Real s86 = -2*s11*s13*s6;
				const Real s87 = s13*s7;
				const Real s88 = 2*s10*s14*s2*s3;
				const Real s89 = -2*s10*s14*s19;
				const Real s90 = 2*s0*s10*s14*s8;
				const Real s91 = -2*s1*s10*s14;
				const Real s92 = -2*s14*s26*s4;
				const Real s93 = -2*s14*s18*s4;
				const Real s94 = -2*s14*s23*s4;
				const Real s95 = 2*s14*s2*s3*s4;
				const Real s96 = 2*s0*s14*s4*s8;
				const Real s97 = 2*s10*s11*s14*s6;
				const Real s98 = 2*s11*s14*s4*s6;
				const Real s99 = -2*s10*s14*s7;
				const Real s100 = -2*s10*s14*s3*s9;
				const Real s101 = 2*s10*s14*s2*s9;
				const Real s102 = 2*s14*s3*s4*s9;
				const Real s103 = -2*s14*s2*s4*s9;
				const Real s104 = -2*s10*s12*s14*s8;
				const Real s105 = 2*s0*s10*s12*s14;
				const Real s106 = 2*s12*s14*s4*s8;
				const Real s107 = -2*s0*s12*s14*s4;
				const Real s108 = s15*s26;
				const Real s109 = s15*s18;
				const Real s110 = s15*s23;
				const Real s111 = -2*s15*s2*s3;
				const Real s112 = s15*s19;
				const Real s113 = -2*s0*s15*s8;
				const Real s114 = s1*s15;
				const Real s115 = -2*s11*s15*s6;
				const Real s116 = s15*s7;
				const Real s117 = 2*s11*s16*s2*s3;
				const Real s118 = -2*s11*s16*s19;
				const Real s119 = 2*s0*s11*s16*s8;
				const Real s120 = -2*s1*s11*s16;
				const Real s121 = 2*s10*s11*s16*s4;
				const Real s122 = -2*s11*s16*s5;
				const Real s123 = -2*s16*s26*s6;
				const Real s124 = -2*s16*s18*s6;
				const Real s125 = -2*s16*s21*s6;
				const Real s126 = 2*s16*s2*s3*s6;
				const Real s127 = 2*s0*s16*s6*s8;
				const Real s128 = 2*s10*s16*s4*s6;
				const Real s129 = -2*s11*s16*s3*s9;
				const Real s130 = 2*s11*s16*s2*s9;
				const Real s131 = 2*s16*s3*s6*s9;
				const Real s132 = -2*s16*s2*s6*s9;
				const Real s133 = -2*s11*s12*s16*s8;
				const Real s134 = 2*s0*s11*s12*s16;
				const Real s135 = 2*s12*s16*s6*s8;
				const Real s136 = -2*s0*s12*s16*s6;
				const Real s137 = -2*s10*s11*s14*s16;
				const Real s138 = 2*s11*s14*s16*s4;
				const Real s139 = 2*s10*s14*s16*s6;
				const Real s140 = -2*s14*s16*s4*s6;
				const Real s141 = s17*s26;
				const Real s142 = s17*s18;
				const Real s143 = s17*s21;
				const Real s144 = -2*s17*s2*s3;
				const Real s145 = s17*s19;
				const Real s146 = -2*s0*s17*s8;
				const Real s147 = s1*s17;
				const Real s148 = -2*s10*s17*s4;
				const Real s149 = s17*s5;
				const Real s150 = s100 + s101 + s102 + s103 + s104 + s105 + s106 + s107 + s108 + s109 + s110 + s111 + s112 + s113 + s114 + s115 + s116 + s117 + s118 + s119 + s120 + s121 + s122 + s123 + s124 + s125 + s126 + s127 + s128 + s129 + s130 + s131 + s132 + s133 + s134 + s135 + s136 + s137 + s138 + s139 + s140 + s141 + s142 + s143 + s144 + s145 + s146 + s147 + s148 + s149 + s20 + s22 + s24 + s25 + s27 + s28 + s29 + s30 + s31 + s32 + s33 + s34 + s35 + s36 + s37 + s38 + s39 + s40 + s41 + s42 + s43 + s44 + s45 + s46 + s47 + s48 + s49 + s50 + s51 + s52 + s54 + s55 + s56 + s57 + s58 + s59 + s60 + s61 + s62 + s63 + s64 + s65 + s66 + s67 + s68 + s69 + s70 + s71 + s72 + s73 + s74 + s75 + s76 + s77 + s78 + s79 + s80 + s81 + s82 + s83 + s84 + s85 + s86 + s87 + s88 + s89 + s90 + s91 + s92 + s93 + s94 + s95 + s96 + s97 + s98 + s99;
				const Real s151 = sqrt(s150);
				const Real s152 = 1/s151;
				const Real s153 = -2*s0*s2*s8;
				const Real s154 = 2*s1*s3;
				const Real s155 = -2*s10*s2*s4;
				const Real s156 = 2*s3*s5;
				const Real s157 = -2*s11*s2*s6;
				const Real s158 = 2*s3*s7;
				const Real s159 = 2*s0*s8*s9;
				const Real s160 = -2*s1*s9;
				const Real s161 = 2*s10*s4*s9;
				const Real s162 = -2*s5*s9;
				const Real s163 = 2*s11*s6*s9;
				const Real s164 = -2*s7*s9;
				const Real s165 = 2*s12*s2*s8;
				const Real s166 = -4*s0*s12*s3;
				const Real s167 = 2*s0*s12*s2;
				const Real s168 = -2*s12*s8*s9;
				const Real s169 = 2*s0*s12*s9;
				const Real s170 = 2*s13*s3;
				const Real s171 = -2*s13*s2;
				const Real s172 = 2*s10*s14*s2;
				const Real s173 = -4*s14*s3*s4;
				const Real s174 = 2*s14*s2*s4;
				const Real s175 = -2*s10*s14*s9;
				const Real s176 = 2*s14*s4*s9;
				const Real s177 = 2*s15*s3;
				const Real s178 = -2*s15*s2;
				const Real s179 = 2*s11*s16*s2;
				const Real s180 = -4*s16*s3*s6;
				const Real s181 = 2*s16*s2*s6;
				const Real s182 = -2*s11*s16*s9;
				const Real s183 = 2*s16*s6*s9;
				const Real s184 = 2*s17*s3;
				const Real s185 = -2*s17*s2;
				const Real s186 = s153 + s154 + s155 + s156 + s157 + s158 + s159 + s160 + s161 + s162 + s163 + s164 + s165 + s166 + s167 + s168 + s169 + s170 + s171 + s172 + s173 + s174 + s175 + s176 + s177 + s178 + s179 + s180 + s181 + s182 + s183 + s184 + s185;
				const Real s187 = -s3;
				const Real s188 = s187 + s2;
				const Real s189 = s188*s188;
				const Real s190 = -s8;
				const Real s191 = s0 + s190;
				const Real s192 = s191*s191;
				const Real s193 = -s10;
				const Real s194 = s193 + s4;
				const Real s195 = s194*s194;
				const Real s196 = -s11;
				const Real s197 = s196 + s6;
				const Real s198 = s197*s197;
				const Real s199 = s187 + s9;
				const Real s200 = s188*s199;
				const Real s201 = s12 + s190;
				const Real s202 = s191*s201;
				const Real s203 = s14 + s193;
				const Real s204 = s194*s203;
				const Real s205 = s16 + s196;
				const Real s206 = s197*s205;
				const Real s207 = s200 + s202 + s204 + s206;
				const Real s208 = s207*s207;
				const Real s209 = -s208;
				const Real s210 = s189 + s192 + s195 + s198;
				const Real s211 = s199*s199;
				const Real s212 = s201*s201;
				const Real s213 = s203*s203;
				const Real s214 = s205*s205;
				const Real s215 = s211 + s212 + s213 + s214;
				const Real s216 = s210*s215;
				const Real s217 = s209 + s216;
				const Real s218 = 1/s217;
				const Real s219 = -(s188*s199);
				const Real s220 = -(s191*s201);
				const Real s221 = -(s194*s203);
				const Real s222 = -(s197*s205);
				const Real s223 = s219 + s220 + s221 + s222;
				const Real s224 = s217*s217;
				const Real s225 = 1/s224;
				const Real s226 = -2*s199*s210;
				const Real s227 = 2*s3;
				const Real s228 = -s2;
				const Real s229 = -s9;
				const Real s230 = s227 + s228 + s229;
				const Real s231 = -2*s207*s230;
				const Real s232 = -2*s188*s215;
				const Real s233 = s226 + s231 + s232;
				const Real s234 = -2*s188*s199*s218;
				const Real s235 = -2*s3;
				const Real s236 = s2 + s235 + s9;
				const Real s237 = -(s218*s223);
				const Real s238 = -(s199*s210*s225*s233);
				const Real s239 = -(s188*s223*s225*s233);
				const Real s240 = -(s210*s218);
				const Real s241 = s188*s218*s236;
				const Real s242 = s234 + s237 + s238 + s239 + s240 + s241;
				const Real s243 = -(s199*s223*s225*s233);
				const Real s244 = -(s188*s215*s225*s233);
				const Real s245 = s199*s218*s236;
				const Real s246 = -(s215*s218);
				const Real s247 = s234 + s237 + s243 + s244 + s245 + s246;
				const Real s248 = s199*s210*s218;
				const Real s249 = s188*s218*s223;
				const Real s250 = s248 + s249;
				const Real s251 = s199*s218*s223;
				const Real s252 = s188*s215*s218;
				const Real s253 = s251 + s252;
				const Real s254 = -(s201*s210*s225*s233);
				const Real s255 = -(s191*s223*s225*s233);
				const Real s256 = s191*s218*s236;
				const Real s257 = -2*s188*s201*s218;
				const Real s258 = s254 + s255 + s256 + s257;
				const Real s259 = -(s201*s223*s225*s233);
				const Real s260 = -(s191*s215*s225*s233);
				const Real s261 = -2*s191*s199*s218;
				const Real s262 = s201*s218*s236;
				const Real s263 = s259 + s260 + s261 + s262;
				const Real s264 = s201*s210*s218;
				const Real s265 = s191*s218*s223;
				const Real s266 = s264 + s265;
				const Real s267 = s201*s218*s223;
				const Real s268 = s191*s215*s218;
				const Real s269 = s267 + s268;
				const Real s270 = -(s203*s210*s225*s233);
				const Real s271 = -(s194*s223*s225*s233);
				const Real s272 = s194*s218*s236;
				const Real s273 = -2*s188*s203*s218;
				const Real s274 = s270 + s271 + s272 + s273;
				const Real s275 = -(s203*s223*s225*s233);
				const Real s276 = -(s194*s215*s225*s233);
				const Real s277 = -2*s194*s199*s218;
				const Real s278 = s203*s218*s236;
				const Real s279 = s275 + s276 + s277 + s278;
				const Real s280 = s203*s210*s218;
				const Real s281 = s194*s218*s223;
				const Real s282 = s280 + s281;
				const Real s283 = s203*s218*s223;
				const Real s284 = s194*s215*s218;
				const Real s285 = s283 + s284;
				const Real s286 = P_D_near__[23*i+0];
				const Real s287 = P_D_near__[23*i+1];
				const Real s288 = 2*s19*s8;
				const Real s289 = -2*s0*s2*s3;
				const Real s290 = -2*s0*s10*s4;
				const Real s291 = 2*s5*s8;
				const Real s292 = -2*s0*s11*s6;
				const Real s293 = 2*s7*s8;
				const Real s294 = -4*s2*s8*s9;
				const Real s295 = 2*s0*s3*s9;
				const Real s296 = 2*s0*s2*s9;
				const Real s297 = 2*s53*s8;
				const Real s298 = -2*s0*s53;
				const Real s299 = 2*s12*s2*s3;
				const Real s300 = -2*s12*s19;
				const Real s301 = 2*s10*s12*s4;
				const Real s302 = -2*s12*s5;
				const Real s303 = 2*s11*s12*s6;
				const Real s304 = -2*s12*s7;
				const Real s305 = -2*s12*s3*s9;
				const Real s306 = 2*s12*s2*s9;
				const Real s307 = 2*s0*s10*s14;
				const Real s308 = -4*s14*s4*s8;
				const Real s309 = 2*s0*s14*s4;
				const Real s310 = -2*s10*s12*s14;
				const Real s311 = 2*s12*s14*s4;
				const Real s312 = 2*s15*s8;
				const Real s313 = -2*s0*s15;
				const Real s314 = 2*s0*s11*s16;
				const Real s315 = -4*s16*s6*s8;
				const Real s316 = 2*s0*s16*s6;
				const Real s317 = -2*s11*s12*s16;
				const Real s318 = 2*s12*s16*s6;
				const Real s319 = 2*s17*s8;
				const Real s320 = -2*s0*s17;
				const Real s321 = s288 + s289 + s290 + s291 + s292 + s293 + s294 + s295 + s296 + s297 + s298 + s299 + s300 + s301 + s302 + s303 + s304 + s305 + s306 + s307 + s308 + s309 + s310 + s311 + s312 + s313 + s314 + s315 + s316 + s317 + s318 + s319 + s320;
				const Real s322 = P_D_near__[23*i+3];
				const Real s323 = P_D_near__[23*i+4];
				const Real s324 = P_D_near__[23*i+5];
				const Real s325 = P_D_near__[23*i+6];
				const Real s326 = P_D_near__[23*i+7];
				const Real s327 = P_D_near__[23*i+8];
				const Real s328 = P_D_near__[23*i+9];
				const Real s329 = P_D_near__[23*i+10];
				const Real s330 = P_D_near__[23*i+11];
				const Real s331 = P_D_near__[23*i+12];
				const Real s332 = P_D_near__[23*i+2];
				const Real s333 = s151/2.;
				const Real s334 = P_D_near__[23*i+13];
				const Real s335 = -2*s201*s210;
				const Real s336 = 2*s8;
				const Real s337 = -s0;
				const Real s338 = -s12;
				const Real s339 = s336 + s337 + s338;
				const Real s340 = -2*s207*s339;
				const Real s341 = -2*s191*s215;
				const Real s342 = s335 + s340 + s341;
				const Real s343 = -2*s8;
				const Real s344 = s0 + s12 + s343;
				const Real s345 = s199*s250;
				const Real s346 = s188*s253;
				const Real s347 = s345 + s346;
				const Real s348 = P_D_near__[23*i+14];
				const Real s349 = -(s199*s210*s218);
				const Real s350 = -(s188*s218*s223);
				const Real s351 = -(s199*s218*s223);
				const Real s352 = -(s188*s215*s218);
				const Real s353 = -(s199*s210*s225*s342);
				const Real s354 = -(s188*s223*s225*s342);
				const Real s355 = s188*s218*s344;
				const Real s356 = s261 + s353 + s354 + s355;
				const Real s357 = -(s199*s223*s225*s342);
				const Real s358 = -(s188*s215*s225*s342);
				const Real s359 = s199*s218*s344;
				const Real s360 = s257 + s357 + s358 + s359;
				const Real s361 = s201*s250;
				const Real s362 = s191*s253;
				const Real s363 = s361 + s362;
				const Real s364 = P_D_near__[23*i+15];
				const Real s365 = s203*s250;
				const Real s366 = s194*s253;
				const Real s367 = s365 + s366;
				const Real s368 = P_D_near__[23*i+16];
				const Real s369 = s205*s250;
				const Real s370 = s197*s253;
				const Real s371 = s369 + s370;
				const Real s372 = P_D_near__[23*i+17];
				const Real s373 = -2*s191*s201*s218;
				const Real s374 = s201*s266;
				const Real s375 = s191*s269;
				const Real s376 = s374 + s375;
				const Real s377 = P_D_near__[23*i+18];
				const Real s378 = -(s201*s210*s225*s342);
				const Real s379 = -(s191*s223*s225*s342);
				const Real s380 = s191*s218*s344;
				const Real s381 = s237 + s240 + s373 + s378 + s379 + s380;
				const Real s382 = -(s201*s223*s225*s342);
				const Real s383 = -(s191*s215*s225*s342);
				const Real s384 = s201*s218*s344;
				const Real s385 = s237 + s246 + s373 + s382 + s383 + s384;
				const Real s386 = s203*s266;
				const Real s387 = s194*s269;
				const Real s388 = s386 + s387;
				const Real s389 = P_D_near__[23*i+19];
				const Real s390 = s205*s266;
				const Real s391 = s197*s269;
				const Real s392 = s390 + s391;
				const Real s393 = P_D_near__[23*i+20];
				const Real s394 = s203*s282;
				const Real s395 = s194*s285;
				const Real s396 = s394 + s395;
				const Real s397 = P_D_near__[23*i+21];
				const Real s398 = -(s203*s210*s225*s342);
				const Real s399 = -(s194*s223*s225*s342);
				const Real s400 = s194*s218*s344;
				const Real s401 = -2*s191*s203*s218;
				const Real s402 = s398 + s399 + s400 + s401;
				const Real s403 = -(s203*s223*s225*s342);
				const Real s404 = -(s194*s215*s225*s342);
				const Real s405 = -2*s194*s201*s218;
				const Real s406 = s203*s218*s344;
				const Real s407 = s403 + s404 + s405 + s406;
				const Real s408 = s205*s282;
				const Real s409 = s197*s285;
				const Real s410 = s408 + s409;
				const Real s411 = P_D_near__[23*i+22];
				const Real s412 = s205*s210*s218;
				const Real s413 = s197*s218*s223;
				const Real s414 = s412 + s413;
				const Real s415 = s205*s414;
				const Real s416 = s205*s218*s223;
				const Real s417 = s197*s215*s218;
				const Real s418 = s416 + s417;
				const Real s419 = s197*s418;
				const Real s420 = s415 + s419;
				const Real s421 = 2*s10*s19;
				const Real s422 = 2*s1*s10;
				const Real s423 = -2*s2*s3*s4;
				const Real s424 = -2*s0*s4*s8;
				const Real s425 = -2*s11*s4*s6;
				const Real s426 = 2*s10*s7;
				const Real s427 = -4*s10*s2*s9;
				const Real s428 = 2*s3*s4*s9;
				const Real s429 = 2*s2*s4*s9;
				const Real s430 = 2*s10*s53;
				const Real s431 = -2*s4*s53;
				const Real s432 = -4*s0*s10*s12;
				const Real s433 = 2*s12*s4*s8;
				const Real s434 = 2*s0*s12*s4;
				const Real s435 = 2*s10*s13;
				const Real s436 = -2*s13*s4;
				const Real s437 = 2*s14*s2*s3;
				const Real s438 = -2*s14*s19;
				const Real s439 = 2*s0*s14*s8;
				const Real s440 = -2*s1*s14;
				const Real s441 = 2*s11*s14*s6;
				const Real s442 = -2*s14*s7;
				const Real s443 = -2*s14*s3*s9;
				const Real s444 = 2*s14*s2*s9;
				const Real s445 = -2*s12*s14*s8;
				const Real s446 = 2*s0*s12*s14;
				const Real s447 = 2*s11*s16*s4;
				const Real s448 = -4*s10*s16*s6;
				const Real s449 = 2*s16*s4*s6;
				const Real s450 = -2*s11*s14*s16;
				const Real s451 = 2*s14*s16*s6;
				const Real s452 = 2*s10*s17;
				const Real s453 = -2*s17*s4;
				const Real s454 = s421 + s422 + s423 + s424 + s425 + s426 + s427 + s428 + s429 + s430 + s431 + s432 + s433 + s434 + s435 + s436 + s437 + s438 + s439 + s440 + s441 + s442 + s443 + s444 + s445 + s446 + s447 + s448 + s449 + s450 + s451 + s452 + s453;
				const Real s455 = -2*s203*s210;
				const Real s456 = 2*s10;
				const Real s457 = -s4;
				const Real s458 = -s14;
				const Real s459 = s456 + s457 + s458;
				const Real s460 = -2*s207*s459;
				const Real s461 = -2*s194*s215;
				const Real s462 = s455 + s460 + s461;
				const Real s463 = -2*s10;
				const Real s464 = s14 + s4 + s463;
				const Real s465 = -(s199*s210*s225*s462);
				const Real s466 = -(s188*s223*s225*s462);
				const Real s467 = s188*s218*s464;
				const Real s468 = s277 + s465 + s466 + s467;
				const Real s469 = -(s199*s223*s225*s462);
				const Real s470 = -(s188*s215*s225*s462);
				const Real s471 = s199*s218*s464;
				const Real s472 = s273 + s469 + s470 + s471;
				const Real s473 = -(s201*s210*s218);
				const Real s474 = -(s191*s218*s223);
				const Real s475 = -(s201*s218*s223);
				const Real s476 = -(s191*s215*s218);
				const Real s477 = -(s201*s210*s225*s462);
				const Real s478 = -(s191*s223*s225*s462);
				const Real s479 = s191*s218*s464;
				const Real s480 = s405 + s477 + s478 + s479;
				const Real s481 = -(s201*s223*s225*s462);
				const Real s482 = -(s191*s215*s225*s462);
				const Real s483 = s201*s218*s464;
				const Real s484 = s401 + s481 + s482 + s483;
				const Real s485 = -2*s194*s203*s218;
				const Real s486 = -(s203*s210*s225*s462);
				const Real s487 = -(s194*s223*s225*s462);
				const Real s488 = s194*s218*s464;
				const Real s489 = s237 + s240 + s485 + s486 + s487 + s488;
				const Real s490 = -(s203*s223*s225*s462);
				const Real s491 = -(s194*s215*s225*s462);
				const Real s492 = s203*s218*s464;
				const Real s493 = s237 + s246 + s485 + s490 + s491 + s492;
				const Real s494 = 2*s11*s19;
				const Real s495 = 2*s1*s11;
				const Real s496 = 2*s11*s5;
				const Real s497 = -2*s2*s3*s6;
				const Real s498 = -2*s0*s6*s8;
				const Real s499 = -2*s10*s4*s6;
				const Real s500 = -4*s11*s2*s9;
				const Real s501 = 2*s3*s6*s9;
				const Real s502 = 2*s2*s6*s9;
				const Real s503 = 2*s11*s53;
				const Real s504 = -2*s53*s6;
				const Real s505 = -4*s0*s11*s12;
				const Real s506 = 2*s12*s6*s8;
				const Real s507 = 2*s0*s12*s6;
				const Real s508 = 2*s11*s13;
				const Real s509 = -2*s13*s6;
				const Real s510 = -4*s11*s14*s4;
				const Real s511 = 2*s10*s14*s6;
				const Real s512 = 2*s14*s4*s6;
				const Real s513 = 2*s11*s15;
				const Real s514 = -2*s15*s6;
				const Real s515 = 2*s16*s2*s3;
				const Real s516 = -2*s16*s19;
				const Real s517 = 2*s0*s16*s8;
				const Real s518 = -2*s1*s16;
				const Real s519 = 2*s10*s16*s4;
				const Real s520 = -2*s16*s5;
				const Real s521 = -2*s16*s3*s9;
				const Real s522 = 2*s16*s2*s9;
				const Real s523 = -2*s12*s16*s8;
				const Real s524 = 2*s0*s12*s16;
				const Real s525 = -2*s10*s14*s16;
				const Real s526 = 2*s14*s16*s4;
				const Real s527 = s494 + s495 + s496 + s497 + s498 + s499 + s500 + s501 + s502 + s503 + s504 + s505 + s506 + s507 + s508 + s509 + s510 + s511 + s512 + s513 + s514 + s515 + s516 + s517 + s518 + s519 + s520 + s521 + s522 + s523 + s524 + s525 + s526;
				const Real s528 = -2*s205*s210;
				const Real s529 = 2*s11;
				const Real s530 = -s6;
				const Real s531 = -s16;
				const Real s532 = s529 + s530 + s531;
				const Real s533 = -2*s207*s532;
				const Real s534 = -2*s197*s215;
				const Real s535 = s528 + s533 + s534;
				const Real s536 = -2*s197*s199*s218;
				const Real s537 = -2*s188*s205*s218;
				const Real s538 = -2*s11;
				const Real s539 = s16 + s538 + s6;
				const Real s540 = -(s199*s210*s225*s535);
				const Real s541 = -(s188*s223*s225*s535);
				const Real s542 = s188*s218*s539;
				const Real s543 = s536 + s540 + s541 + s542;
				const Real s544 = -(s199*s223*s225*s535);
				const Real s545 = -(s188*s215*s225*s535);
				const Real s546 = s199*s218*s539;
				const Real s547 = s537 + s544 + s545 + s546;
				const Real s548 = -2*s197*s201*s218;
				const Real s549 = -2*s191*s205*s218;
				const Real s550 = -(s201*s210*s225*s535);
				const Real s551 = -(s191*s223*s225*s535);
				const Real s552 = s191*s218*s539;
				const Real s553 = s548 + s550 + s551 + s552;
				const Real s554 = -(s201*s223*s225*s535);
				const Real s555 = -(s191*s215*s225*s535);
				const Real s556 = s201*s218*s539;
				const Real s557 = s549 + s554 + s555 + s556;
				const Real s558 = -2*s197*s203*s218;
				const Real s559 = -2*s194*s205*s218;
				const Real s560 = -(s203*s210*s218);
				const Real s561 = -(s194*s218*s223);
				const Real s562 = -(s203*s218*s223);
				const Real s563 = -(s194*s215*s218);
				const Real s564 = -(s203*s210*s225*s535);
				const Real s565 = -(s194*s223*s225*s535);
				const Real s566 = s194*s218*s539;
				const Real s567 = s558 + s564 + s565 + s566;
				const Real s568 = -(s203*s223*s225*s535);
				const Real s569 = -(s194*s215*s225*s535);
				const Real s570 = s203*s218*s539;
				const Real s571 = s559 + s568 + s569 + s570;
				const Real s572 = -2*s197*s205*s218;
				const Real s573 = 2*s18*s2;
				const Real s574 = 2*s2*s21;
				const Real s575 = 2*s2*s23;
				const Real s576 = -2*s0*s3*s8;
				const Real s577 = -2*s10*s3*s4;
				const Real s578 = -2*s11*s3*s6;
				const Real s579 = -2*s18*s9;
				const Real s580 = -2*s21*s9;
				const Real s581 = -2*s23*s9;
				const Real s582 = 2*s12*s3*s8;
				const Real s583 = -4*s12*s2*s8;
				const Real s584 = 2*s0*s12*s3;
				const Real s585 = 2*s12*s8*s9;
				const Real s586 = -2*s0*s12*s9;
				const Real s587 = -2*s13*s3;
				const Real s588 = 2*s13*s2;
				const Real s589 = 2*s10*s14*s3;
				const Real s590 = -4*s10*s14*s2;
				const Real s591 = 2*s14*s3*s4;
				const Real s592 = 2*s10*s14*s9;
				const Real s593 = -2*s14*s4*s9;
				const Real s594 = -2*s15*s3;
				const Real s595 = 2*s15*s2;
				const Real s596 = 2*s11*s16*s3;
				const Real s597 = -4*s11*s16*s2;
				const Real s598 = 2*s16*s3*s6;
				const Real s599 = 2*s11*s16*s9;
				const Real s600 = -2*s16*s6*s9;
				const Real s601 = -2*s17*s3;
				const Real s602 = 2*s17*s2;
				const Real s603 = s159 + s161 + s163 + s573 + s574 + s575 + s576 + s577 + s578 + s579 + s580 + s581 + s582 + s583 + s584 + s585 + s586 + s587 + s588 + s589 + s590 + s591 + s592 + s593 + s594 + s595 + s596 + s597 + s598 + s599 + s600 + s601 + s602;
				const Real s604 = -2*s199*s207;
				const Real s605 = 2*s188*s215;
				const Real s606 = s604 + s605;
				const Real s607 = -(s199*s210*s225*s606);
				const Real s608 = -(s188*s223*s225*s606);
				const Real s609 = s188*s199*s218;
				const Real s610 = s218*s223;
				const Real s611 = s607 + s608 + s609 + s610;
				const Real s612 = -(s199*s223*s225*s606);
				const Real s613 = -(s188*s215*s225*s606);
				const Real s614 = -(s211*s218);
				const Real s615 = s215*s218;
				const Real s616 = s612 + s613 + s614 + s615;
				const Real s617 = -(s201*s210*s225*s606);
				const Real s618 = -(s191*s223*s225*s606);
				const Real s619 = -(s191*s199*s218);
				const Real s620 = 2*s188*s201*s218;
				const Real s621 = s617 + s618 + s619 + s620;
				const Real s622 = -(s201*s223*s225*s606);
				const Real s623 = -(s191*s215*s225*s606);
				const Real s624 = -(s199*s201*s218);
				const Real s625 = s622 + s623 + s624;
				const Real s626 = -(s203*s210*s225*s606);
				const Real s627 = -(s194*s223*s225*s606);
				const Real s628 = -(s194*s199*s218);
				const Real s629 = 2*s188*s203*s218;
				const Real s630 = s626 + s627 + s628 + s629;
				const Real s631 = -(s203*s223*s225*s606);
				const Real s632 = -(s194*s215*s225*s606);
				const Real s633 = -(s199*s203*s218);
				const Real s634 = s631 + s632 + s633;
				const Real s635 = -2*s2*s3*s8;
				const Real s636 = 2*s0*s26;
				const Real s637 = 2*s0*s21;
				const Real s638 = 2*s0*s23;
				const Real s639 = -2*s10*s4*s8;
				const Real s640 = -2*s11*s6*s8;
				const Real s641 = 2*s3*s8*s9;
				const Real s642 = 2*s2*s8*s9;
				const Real s643 = -4*s0*s3*s9;
				const Real s644 = -2*s53*s8;
				const Real s645 = 2*s0*s53;
				const Real s646 = -2*s12*s26;
				const Real s647 = -2*s12*s21;
				const Real s648 = -2*s12*s23;
				const Real s649 = 2*s12*s3*s9;
				const Real s650 = -2*s12*s2*s9;
				const Real s651 = 2*s10*s14*s8;
				const Real s652 = -4*s0*s10*s14;
				const Real s653 = 2*s14*s4*s8;
				const Real s654 = 2*s10*s12*s14;
				const Real s655 = -2*s12*s14*s4;
				const Real s656 = -2*s15*s8;
				const Real s657 = 2*s0*s15;
				const Real s658 = 2*s11*s16*s8;
				const Real s659 = -4*s0*s11*s16;
				const Real s660 = 2*s16*s6*s8;
				const Real s661 = 2*s11*s12*s16;
				const Real s662 = -2*s12*s16*s6;
				const Real s663 = -2*s17*s8;
				const Real s664 = 2*s0*s17;
				const Real s665 = s299 + s301 + s303 + s635 + s636 + s637 + s638 + s639 + s640 + s641 + s642 + s643 + s644 + s645 + s646 + s647 + s648 + s649 + s650 + s651 + s652 + s653 + s654 + s655 + s656 + s657 + s658 + s659 + s660 + s661 + s662 + s663 + s664;
				const Real s666 = -2*s201*s207;
				const Real s667 = 2*s191*s215;
				const Real s668 = s666 + s667;
				const Real s669 = -(s199*s210*s225*s668);
				const Real s670 = -(s188*s223*s225*s668);
				const Real s671 = 2*s191*s199*s218;
				const Real s672 = -(s188*s201*s218);
				const Real s673 = s669 + s670 + s671 + s672;
				const Real s674 = -(s199*s223*s225*s668);
				const Real s675 = -(s188*s215*s225*s668);
				const Real s676 = s624 + s674 + s675;
				const Real s677 = -(s201*s210*s225*s668);
				const Real s678 = -(s191*s223*s225*s668);
				const Real s679 = s191*s201*s218;
				const Real s680 = s610 + s677 + s678 + s679;
				const Real s681 = -(s201*s223*s225*s668);
				const Real s682 = -(s191*s215*s225*s668);
				const Real s683 = -(s212*s218);
				const Real s684 = s615 + s681 + s682 + s683;
				const Real s685 = -(s203*s210*s225*s668);
				const Real s686 = -(s194*s223*s225*s668);
				const Real s687 = -(s194*s201*s218);
				const Real s688 = 2*s191*s203*s218;
				const Real s689 = s685 + s686 + s687 + s688;
				const Real s690 = -(s203*s223*s225*s668);
				const Real s691 = -(s194*s215*s225*s668);
				const Real s692 = -(s201*s203*s218);
				const Real s693 = s690 + s691 + s692;
				const Real s694 = -2*s10*s2*s3;
				const Real s695 = -2*s0*s10*s8;
				const Real s696 = 2*s26*s4;
				const Real s697 = 2*s18*s4;
				const Real s698 = 2*s23*s4;
				const Real s699 = -2*s10*s11*s6;
				const Real s700 = 2*s10*s3*s9;
				const Real s701 = 2*s10*s2*s9;
				const Real s702 = -4*s3*s4*s9;
				const Real s703 = -2*s10*s53;
				const Real s704 = 2*s4*s53;
				const Real s705 = 2*s10*s12*s8;
				const Real s706 = 2*s0*s10*s12;
				const Real s707 = -4*s12*s4*s8;
				const Real s708 = -2*s10*s13;
				const Real s709 = 2*s13*s4;
				const Real s710 = -2*s14*s26;
				const Real s711 = -2*s14*s18;
				const Real s712 = -2*s14*s23;
				const Real s713 = 2*s14*s3*s9;
				const Real s714 = -2*s14*s2*s9;
				const Real s715 = 2*s12*s14*s8;
				const Real s716 = -2*s0*s12*s14;
				const Real s717 = 2*s10*s11*s16;
				const Real s718 = -4*s11*s16*s4;
				const Real s719 = 2*s10*s16*s6;
				const Real s720 = 2*s11*s14*s16;
				const Real s721 = -2*s14*s16*s6;
				const Real s722 = -2*s10*s17;
				const Real s723 = 2*s17*s4;
				const Real s724 = s437 + s439 + s441 + s694 + s695 + s696 + s697 + s698 + s699 + s700 + s701 + s702 + s703 + s704 + s705 + s706 + s707 + s708 + s709 + s710 + s711 + s712 + s713 + s714 + s715 + s716 + s717 + s718 + s719 + s720 + s721 + s722 + s723;
				const Real s725 = -2*s203*s207;
				const Real s726 = 2*s194*s215;
				const Real s727 = s725 + s726;
				const Real s728 = -(s199*s210*s225*s727);
				const Real s729 = -(s188*s223*s225*s727);
				const Real s730 = 2*s194*s199*s218;
				const Real s731 = -(s188*s203*s218);
				const Real s732 = s728 + s729 + s730 + s731;
				const Real s733 = -(s199*s223*s225*s727);
				const Real s734 = -(s188*s215*s225*s727);
				const Real s735 = s633 + s733 + s734;
				const Real s736 = -(s201*s210*s225*s727);
				const Real s737 = -(s191*s223*s225*s727);
				const Real s738 = 2*s194*s201*s218;
				const Real s739 = -(s191*s203*s218);
				const Real s740 = s736 + s737 + s738 + s739;
				const Real s741 = -(s201*s223*s225*s727);
				const Real s742 = -(s191*s215*s225*s727);
				const Real s743 = s692 + s741 + s742;
				const Real s744 = -(s203*s210*s225*s727);
				const Real s745 = -(s194*s223*s225*s727);
				const Real s746 = s194*s203*s218;
				const Real s747 = s610 + s744 + s745 + s746;
				const Real s748 = -(s203*s223*s225*s727);
				const Real s749 = -(s194*s215*s225*s727);
				const Real s750 = -(s213*s218);
				const Real s751 = s615 + s748 + s749 + s750;
				const Real s752 = -2*s11*s2*s3;
				const Real s753 = -2*s0*s11*s8;
				const Real s754 = -2*s10*s11*s4;
				const Real s755 = 2*s26*s6;
				const Real s756 = 2*s18*s6;
				const Real s757 = 2*s21*s6;
				const Real s758 = 2*s11*s3*s9;
				const Real s759 = 2*s11*s2*s9;
				const Real s760 = -4*s3*s6*s9;
				const Real s761 = -2*s11*s53;
				const Real s762 = 2*s53*s6;
				const Real s763 = 2*s11*s12*s8;
				const Real s764 = 2*s0*s11*s12;
				const Real s765 = -4*s12*s6*s8;
				const Real s766 = -2*s11*s13;
				const Real s767 = 2*s13*s6;
				const Real s768 = 2*s10*s11*s14;
				const Real s769 = 2*s11*s14*s4;
				const Real s770 = -4*s10*s14*s6;
				const Real s771 = -2*s11*s15;
				const Real s772 = 2*s15*s6;
				const Real s773 = -2*s16*s26;
				const Real s774 = -2*s16*s18;
				const Real s775 = -2*s16*s21;
				const Real s776 = 2*s16*s3*s9;
				const Real s777 = -2*s16*s2*s9;
				const Real s778 = 2*s12*s16*s8;
				const Real s779 = -2*s0*s12*s16;
				const Real s780 = 2*s10*s14*s16;
				const Real s781 = -2*s14*s16*s4;
				const Real s782 = s515 + s517 + s519 + s752 + s753 + s754 + s755 + s756 + s757 + s758 + s759 + s760 + s761 + s762 + s763 + s764 + s765 + s766 + s767 + s768 + s769 + s770 + s771 + s772 + s773 + s774 + s775 + s776 + s777 + s778 + s779 + s780 + s781;
				const Real s783 = -2*s205*s207;
				const Real s784 = 2*s197*s215;
				const Real s785 = s783 + s784;
				const Real s786 = -(s199*s205*s218);
				const Real s787 = -(s199*s210*s225*s785);
				const Real s788 = -(s188*s223*s225*s785);
				const Real s789 = 2*s197*s199*s218;
				const Real s790 = -(s188*s205*s218);
				const Real s791 = s787 + s788 + s789 + s790;
				const Real s792 = -(s199*s223*s225*s785);
				const Real s793 = -(s188*s215*s225*s785);
				const Real s794 = s786 + s792 + s793;
				const Real s795 = -(s201*s205*s218);
				const Real s796 = -(s201*s210*s225*s785);
				const Real s797 = -(s191*s223*s225*s785);
				const Real s798 = 2*s197*s201*s218;
				const Real s799 = -(s191*s205*s218);
				const Real s800 = s796 + s797 + s798 + s799;
				const Real s801 = -(s201*s223*s225*s785);
				const Real s802 = -(s191*s215*s225*s785);
				const Real s803 = s795 + s801 + s802;
				const Real s804 = -(s203*s205*s218);
				const Real s805 = -(s203*s210*s225*s785);
				const Real s806 = -(s194*s223*s225*s785);
				const Real s807 = 2*s197*s203*s218;
				const Real s808 = -(s194*s205*s218);
				const Real s809 = s805 + s806 + s807 + s808;
				const Real s810 = -(s203*s223*s225*s785);
				const Real s811 = -(s194*s215*s225*s785);
				const Real s812 = s804 + s810 + s811;
				const Real s813 = -2*s18*s2;
				const Real s814 = -2*s2*s21;
				const Real s815 = -2*s2*s23;
				const Real s816 = 2*s0*s3*s8;
				const Real s817 = 2*s0*s2*s8;
				const Real s818 = -2*s1*s3;
				const Real s819 = 2*s10*s3*s4;
				const Real s820 = 2*s10*s2*s4;
				const Real s821 = -2*s3*s5;
				const Real s822 = 2*s11*s3*s6;
				const Real s823 = 2*s11*s2*s6;
				const Real s824 = -2*s3*s7;
				const Real s825 = 2*s18*s9;
				const Real s826 = 2*s21*s9;
				const Real s827 = 2*s23*s9;
				const Real s828 = -4*s0*s8*s9;
				const Real s829 = 2*s1*s9;
				const Real s830 = -4*s10*s4*s9;
				const Real s831 = 2*s5*s9;
				const Real s832 = -4*s11*s6*s9;
				const Real s833 = 2*s7*s9;
				const Real s834 = -2*s12*s3*s8;
				const Real s835 = -2*s0*s12*s2;
				const Real s836 = -2*s10*s14*s3;
				const Real s837 = -2*s14*s2*s4;
				const Real s838 = -2*s11*s16*s3;
				const Real s839 = -2*s16*s2*s6;
				const Real s840 = s165 + s172 + s179 + s584 + s591 + s598 + s813 + s814 + s815 + s816 + s817 + s818 + s819 + s820 + s821 + s822 + s823 + s824 + s825 + s826 + s827 + s828 + s829 + s830 + s831 + s832 + s833 + s834 + s835 + s836 + s837 + s838 + s839;
				const Real s841 = 2*s199*s210;
				const Real s842 = -2*s188*s207;
				const Real s843 = s841 + s842;
				const Real s844 = -(s199*s210*s225*s843);
				const Real s845 = -(s188*s223*s225*s843);
				const Real s846 = -(s189*s218);
				const Real s847 = s210*s218;
				const Real s848 = s844 + s845 + s846 + s847;
				const Real s849 = -(s199*s223*s225*s843);
				const Real s850 = -(s188*s215*s225*s843);
				const Real s851 = s609 + s610 + s849 + s850;
				const Real s852 = -(s201*s210*s225*s843);
				const Real s853 = -(s191*s223*s225*s843);
				const Real s854 = -(s188*s191*s218);
				const Real s855 = s852 + s853 + s854;
				const Real s856 = -(s201*s223*s225*s843);
				const Real s857 = -(s191*s215*s225*s843);
				const Real s858 = s671 + s672 + s856 + s857;
				const Real s859 = -(s203*s210*s225*s843);
				const Real s860 = -(s194*s223*s225*s843);
				const Real s861 = -(s188*s194*s218);
				const Real s862 = s859 + s860 + s861;
				const Real s863 = -(s203*s223*s225*s843);
				const Real s864 = -(s194*s215*s225*s843);
				const Real s865 = s730 + s731 + s863 + s864;
				const Real s866 = 2*s2*s3*s8;
				const Real s867 = -2*s19*s8;
				const Real s868 = -2*s0*s26;
				const Real s869 = -2*s0*s21;
				const Real s870 = -2*s0*s23;
				const Real s871 = 2*s0*s2*s3;
				const Real s872 = 2*s10*s4*s8;
				const Real s873 = 2*s0*s10*s4;
				const Real s874 = -2*s5*s8;
				const Real s875 = 2*s11*s6*s8;
				const Real s876 = 2*s0*s11*s6;
				const Real s877 = -2*s7*s8;
				const Real s878 = -2*s3*s8*s9;
				const Real s879 = -2*s0*s2*s9;
				const Real s880 = 2*s12*s26;
				const Real s881 = 2*s12*s21;
				const Real s882 = 2*s12*s23;
				const Real s883 = -4*s12*s2*s3;
				const Real s884 = 2*s12*s19;
				const Real s885 = -4*s10*s12*s4;
				const Real s886 = 2*s12*s5;
				const Real s887 = -4*s11*s12*s6;
				const Real s888 = 2*s12*s7;
				const Real s889 = -2*s10*s14*s8;
				const Real s890 = -2*s0*s14*s4;
				const Real s891 = -2*s11*s16*s8;
				const Real s892 = -2*s0*s16*s6;
				const Real s893 = s295 + s307 + s314 + s642 + s653 + s660 + s866 + s867 + s868 + s869 + s870 + s871 + s872 + s873 + s874 + s875 + s876 + s877 + s878 + s879 + s880 + s881 + s882 + s883 + s884 + s885 + s886 + s887 + s888 + s889 + s890 + s891 + s892;
				const Real s894 = 2*s201*s210;
				const Real s895 = -2*s191*s207;
				const Real s896 = s894 + s895;
				const Real s897 = -(s199*s210*s225*s896);
				const Real s898 = -(s188*s223*s225*s896);
				const Real s899 = s854 + s897 + s898;
				const Real s900 = -(s199*s223*s225*s896);
				const Real s901 = -(s188*s215*s225*s896);
				const Real s902 = s619 + s620 + s900 + s901;
				const Real s903 = -(s201*s210*s225*s896);
				const Real s904 = -(s191*s223*s225*s896);
				const Real s905 = -(s192*s218);
				const Real s906 = s847 + s903 + s904 + s905;
				const Real s907 = -(s201*s223*s225*s896);
				const Real s908 = -(s191*s215*s225*s896);
				const Real s909 = s610 + s679 + s907 + s908;
				const Real s910 = -(s203*s210*s225*s896);
				const Real s911 = -(s194*s223*s225*s896);
				const Real s912 = -(s191*s194*s218);
				const Real s913 = s910 + s911 + s912;
				const Real s914 = -(s203*s223*s225*s896);
				const Real s915 = -(s194*s215*s225*s896);
				const Real s916 = s738 + s739 + s914 + s915;
				const Real s917 = 2*s10*s2*s3;
				const Real s918 = -2*s10*s19;
				const Real s919 = 2*s0*s10*s8;
				const Real s920 = -2*s1*s10;
				const Real s921 = -2*s26*s4;
				const Real s922 = -2*s18*s4;
				const Real s923 = -2*s23*s4;
				const Real s924 = 2*s2*s3*s4;
				const Real s925 = 2*s0*s4*s8;
				const Real s926 = 2*s10*s11*s6;
				const Real s927 = 2*s11*s4*s6;
				const Real s928 = -2*s10*s7;
				const Real s929 = -2*s10*s3*s9;
				const Real s930 = -2*s2*s4*s9;
				const Real s931 = -2*s10*s12*s8;
				const Real s932 = -2*s0*s12*s4;
				const Real s933 = 2*s14*s26;
				const Real s934 = 2*s14*s18;
				const Real s935 = 2*s14*s23;
				const Real s936 = -4*s14*s2*s3;
				const Real s937 = 2*s14*s19;
				const Real s938 = -4*s0*s14*s8;
				const Real s939 = 2*s1*s14;
				const Real s940 = -4*s11*s14*s6;
				const Real s941 = 2*s14*s7;
				const Real s942 = -2*s10*s11*s16;
				const Real s943 = -2*s16*s4*s6;
				const Real s944 = s428 + s433 + s447 + s701 + s706 + s719 + s917 + s918 + s919 + s920 + s921 + s922 + s923 + s924 + s925 + s926 + s927 + s928 + s929 + s930 + s931 + s932 + s933 + s934 + s935 + s936 + s937 + s938 + s939 + s940 + s941 + s942 + s943;
				const Real s945 = 2*s203*s210;
				const Real s946 = -2*s194*s207;
				const Real s947 = s945 + s946;
				const Real s948 = -(s199*s210*s225*s947);
				const Real s949 = -(s188*s223*s225*s947);
				const Real s950 = s861 + s948 + s949;
				const Real s951 = -(s199*s223*s225*s947);
				const Real s952 = -(s188*s215*s225*s947);
				const Real s953 = s628 + s629 + s951 + s952;
				const Real s954 = -(s201*s210*s225*s947);
				const Real s955 = -(s191*s223*s225*s947);
				const Real s956 = s912 + s954 + s955;
				const Real s957 = -(s201*s223*s225*s947);
				const Real s958 = -(s191*s215*s225*s947);
				const Real s959 = s687 + s688 + s957 + s958;
				const Real s960 = -(s203*s210*s225*s947);
				const Real s961 = -(s194*s223*s225*s947);
				const Real s962 = -(s195*s218);
				const Real s963 = s847 + s960 + s961 + s962;
				const Real s964 = -(s203*s223*s225*s947);
				const Real s965 = -(s194*s215*s225*s947);
				const Real s966 = s610 + s746 + s964 + s965;
				const Real s967 = 2*s11*s2*s3;
				const Real s968 = -2*s11*s19;
				const Real s969 = 2*s0*s11*s8;
				const Real s970 = -2*s1*s11;
				const Real s971 = 2*s10*s11*s4;
				const Real s972 = -2*s11*s5;
				const Real s973 = -2*s26*s6;
				const Real s974 = -2*s18*s6;
				const Real s975 = -2*s21*s6;
				const Real s976 = 2*s2*s3*s6;
				const Real s977 = 2*s0*s6*s8;
				const Real s978 = 2*s10*s4*s6;
				const Real s979 = -2*s11*s3*s9;
				const Real s980 = -2*s2*s6*s9;
				const Real s981 = -2*s11*s12*s8;
				const Real s982 = -2*s0*s12*s6;
				const Real s983 = -2*s10*s11*s14;
				const Real s984 = -2*s14*s4*s6;
				const Real s985 = 2*s16*s26;
				const Real s986 = 2*s16*s18;
				const Real s987 = 2*s16*s21;
				const Real s988 = -4*s16*s2*s3;
				const Real s989 = 2*s16*s19;
				const Real s990 = -4*s0*s16*s8;
				const Real s991 = 2*s1*s16;
				const Real s992 = -4*s10*s16*s4;
				const Real s993 = 2*s16*s5;
				const Real s994 = s501 + s506 + s511 + s759 + s764 + s769 + s967 + s968 + s969 + s970 + s971 + s972 + s973 + s974 + s975 + s976 + s977 + s978 + s979 + s980 + s981 + s982 + s983 + s984 + s985 + s986 + s987 + s988 + s989 + s990 + s991 + s992 + s993;
				const Real s995 = 2*s205*s210;
				const Real s996 = -2*s197*s207;
				const Real s997 = s995 + s996;
				const Real s998 = -(s188*s197*s218);
				const Real s999 = -(s197*s199*s218);
				const Real s1000 = 2*s188*s205*s218;
				const Real s1001 = -(s199*s210*s225*s997);
				const Real s1002 = -(s188*s223*s225*s997);
				const Real s1003 = s1001 + s1002 + s998;
				const Real s1004 = -(s199*s223*s225*s997);
				const Real s1005 = -(s188*s215*s225*s997);
				const Real s1006 = s1000 + s1004 + s1005 + s999;
				const Real s1007 = -(s191*s197*s218);
				const Real s1008 = -(s197*s201*s218);
				const Real s1009 = 2*s191*s205*s218;
				const Real s1010 = -(s201*s210*s225*s997);
				const Real s1011 = -(s191*s223*s225*s997);
				const Real s1012 = s1007 + s1010 + s1011;
				const Real s1013 = -(s201*s223*s225*s997);
				const Real s1014 = -(s191*s215*s225*s997);
				const Real s1015 = s1008 + s1009 + s1013 + s1014;
				const Real s1016 = -(s194*s197*s218);
				const Real s1017 = -(s197*s203*s218);
				const Real s1018 = 2*s194*s205*s218;
				const Real s1019 = -(s203*s210*s225*s997);
				const Real s1020 = -(s194*s223*s225*s997);
				const Real s1021 = s1016 + s1019 + s1020;
				const Real s1022 = -(s203*s223*s225*s997);
				const Real s1023 = -(s194*s215*s225*s997);
				const Real s1024 = s1017 + s1018 + s1022 + s1023;
				const Real s1025 = s197*s205*s218;
				buffer__[12*i+0] += (s152*s186*s286)/4. + (s10*s152*s186*s322)/4. + (s11*s152*s186*s323)/4. + (s152*s186*s2*s324)/4. + (s0*s152*s186*s325)/4. + (s12*s152*s186*s329)/4. + (s14*s152*s186*s330)/4. + (s152*s16*s186*s331)/4. + s287*((s152*s186*s3)/4. + s333) + s334*((s152*s186*s347)/4. + (s151*(s199*s242 + s188*s247 + s349 + s350 + s351 + s352))/2.) + s348*((s151*(s201*s242 + s191*s247))/2. + (s152*s186*s363)/4.) + s364*((s151*(s203*s242 + s194*s247))/2. + (s152*s186*s367)/4.) + s368*((s151*(s205*s242 + s197*s247))/2. + (s152*s186*s371)/4.) + s372*((s151*(s201*s258 + s191*s263))/2. + (s152*s186*s376)/4.) + s377*((s151*(s203*s258 + s194*s263))/2. + (s152*s186*s388)/4.) + s389*((s151*(s205*s258 + s197*s263))/2. + (s152*s186*s392)/4.) + s393*((s151*(s203*s274 + s194*s279))/2. + (s152*s186*s396)/4.) + (s152*s186*s326*s4)/4. + s397*((s151*(s205*s274 + s197*s279))/2. + (s152*s186*s410)/4.) + s411*((s152*s186*s420)/4. + (s151*(s197*(-(s197*s215*s225*s233) - s205*s223*s225*s233 + s205*s218*s236 + s536) + s205*(-(s205*s210*s225*s233) - s197*s223*s225*s233 + s197*s218*s236 + s537)))/2.) + (s152*s186*s327*s6)/4. + (s152*s186*s332*s8)/4. + (s152*s186*s328*s9)/4.;
				buffer__[12*i+1] += (s152*s286*s321)/4. + (s152*s287*s3*s321)/4. + (s10*s152*s321*s322)/4. + (s11*s152*s321*s323)/4. + (s152*s2*s321*s324)/4. + (s0*s152*s321*s325)/4. + (s12*s152*s321*s329)/4. + (s14*s152*s321*s330)/4. + (s152*s16*s321*s331)/4. + s334*((s152*s321*s347)/4. + (s151*(s199*s356 + s188*s360))/2.) + s348*((s151*(s349 + s350 + s351 + s352 + s201*s356 + s191*s360))/2. + (s152*s321*s363)/4.) + s364*((s151*(s203*s356 + s194*s360))/2. + (s152*s321*s367)/4.) + s368*((s151*(s205*s356 + s197*s360))/2. + (s152*s321*s371)/4.) + s377*((s151*(s203*s381 + s194*s385))/2. + (s152*s321*s388)/4.) + s389*((s151*(s205*s381 + s197*s385))/2. + (s152*s321*s392)/4.) + (s152*s321*s326*s4)/4. + s393*((s152*s321*s396)/4. + (s151*(s203*s402 + s194*s407))/2.) + s397*((s151*(s205*s402 + s197*s407))/2. + (s152*s321*s410)/4.) + s372*((s152*s321*s376)/4. + (s151*(s201*s381 + s191*s385 + s473 + s474 + s475 + s476))/2.) + s411*((s152*s321*s420)/4. + (s151*(s197*(-(s197*s215*s225*s342) - s205*s223*s225*s342 + s205*s218*s344 + s548) + s205*(-(s205*s210*s225*s342) - s197*s223*s225*s342 + s197*s218*s344 + s549)))/2.) + (s152*s321*s327*s6)/4. + s332*(s333 + (s152*s321*s8)/4.) + (s152*s321*s328*s9)/4.;
				buffer__[12*i+2] += (s152*s286*s454)/4. + (s152*s287*s3*s454)/4. + (s11*s152*s323*s454)/4. + (s152*s2*s324*s454)/4. + (s0*s152*s325*s454)/4. + (s12*s152*s329*s454)/4. + (s14*s152*s330*s454)/4. + (s152*s16*s331*s454)/4. + (s152*s326*s4*s454)/4. + s322*(s333 + (s10*s152*s454)/4.) + s334*((s152*s347*s454)/4. + (s151*(s199*s468 + s188*s472))/2.) + s348*((s152*s363*s454)/4. + (s151*(s201*s468 + s191*s472))/2.) + s364*((s152*s367*s454)/4. + (s151*(s349 + s350 + s351 + s352 + s203*s468 + s194*s472))/2.) + s368*((s152*s371*s454)/4. + (s151*(s205*s468 + s197*s472))/2.) + s372*((s152*s376*s454)/4. + (s151*(s201*s480 + s191*s484))/2.) + s377*((s152*s388*s454)/4. + (s151*(s473 + s474 + s475 + s476 + s203*s480 + s194*s484))/2.) + s389*((s152*s392*s454)/4. + (s151*(s205*s480 + s197*s484))/2.) + s397*((s152*s410*s454)/4. + (s151*(s205*s489 + s197*s493))/2.) + s411*((s152*s420*s454)/4. + (s151*(s197*(-(s197*s215*s225*s462) - s205*s223*s225*s462 + s205*s218*s464 + s558) + s205*(-(s205*s210*s225*s462) - s197*s223*s225*s462 + s197*s218*s464 + s559)))/2.) + s393*((s152*s396*s454)/4. + (s151*(s203*s489 + s194*s493 + s560 + s561 + s562 + s563))/2.) + (s152*s327*s454*s6)/4. + (s152*s332*s454*s8)/4. + (s152*s328*s454*s9)/4.;
				buffer__[12*i+3] += (s152*s286*s527)/4. + (s152*s287*s3*s527)/4. + (s10*s152*s322*s527)/4. + (s152*s2*s324*s527)/4. + (s0*s152*s325*s527)/4. + (s12*s152*s329*s527)/4. + (s14*s152*s330*s527)/4. + (s152*s16*s331*s527)/4. + (s152*s326*s4*s527)/4. + s323*(s333 + (s11*s152*s527)/4.) + s334*((s152*s347*s527)/4. + (s151*(s199*s543 + s188*s547))/2.) + s348*((s152*s363*s527)/4. + (s151*(s201*s543 + s191*s547))/2.) + s364*((s152*s367*s527)/4. + (s151*(s203*s543 + s194*s547))/2.) + s368*((s152*s371*s527)/4. + (s151*(s349 + s350 + s351 + s352 + s205*s543 + s197*s547))/2.) + s372*((s152*s376*s527)/4. + (s151*(s201*s553 + s191*s557))/2.) + s377*((s152*s388*s527)/4. + (s151*(s203*s553 + s194*s557))/2.) + s389*((s152*s392*s527)/4. + (s151*(s473 + s474 + s475 + s476 + s205*s553 + s197*s557))/2.) + s393*((s152*s396*s527)/4. + (s151*(s203*s567 + s194*s571))/2.) + s397*((s152*s410*s527)/4. + (s151*(s560 + s561 + s562 + s563 + s205*s567 + s197*s571))/2.) + s411*((s152*s420*s527)/4. + (s151*(-(s205*s210*s218) - s197*s215*s218 - s197*s218*s223 - s205*s218*s223 + s205*(s237 + s240 - s205*s210*s225*s535 - s197*s223*s225*s535 + s197*s218*s539 + s572) + s197*(s237 + s246 - s197*s215*s225*s535 - s205*s223*s225*s535 + s205*s218*s539 + s572)))/2.) + (s152*s327*s527*s6)/4. + (s152*s332*s527*s8)/4. + (s152*s328*s527*s9)/4.;
				buffer__[12*i+4] += (s152*s286*s603)/4. + (s152*s287*s3*s603)/4. + (s10*s152*s322*s603)/4. + (s11*s152*s323*s603)/4. + (s0*s152*s325*s603)/4. + (s12*s152*s329*s603)/4. + (s14*s152*s330*s603)/4. + (s152*s16*s331*s603)/4. + (s152*s326*s4*s603)/4. + (s152*s327*s6*s603)/4. + s324*(s333 + (s152*s2*s603)/4.) + s334*((s152*s347*s603)/4. + (s151*(s251 + s252 + s199*s611 + s188*s616))/2.) + s348*((s152*s363*s603)/4. + (s151*(s201*s611 + s191*s616))/2.) + s364*((s152*s367*s603)/4. + (s151*(s203*s611 + s194*s616))/2.) + s368*((s152*s371*s603)/4. + (s151*(s205*s611 + s197*s616))/2.) + s372*((s152*s376*s603)/4. + (s151*(s201*s621 + s191*s625))/2.) + s377*((s152*s388*s603)/4. + (s151*(s203*s621 + s194*s625))/2.) + s389*((s152*s392*s603)/4. + (s151*(s205*s621 + s197*s625))/2.) + s393*((s152*s396*s603)/4. + (s151*(s203*s630 + s194*s634))/2.) + s397*((s152*s410*s603)/4. + (s151*(s205*s630 + s197*s634))/2.) + (s152*s332*s603*s8)/4. + (s152*s328*s603*s9)/4. + s411*((s152*s420*s603)/4. + (s151*(s197*(-(s197*s215*s225*s606) - s205*s223*s225*s606 + s786) + s205*(s1000 - s205*s210*s225*s606 - s197*s223*s225*s606 + s999)))/2.);
				buffer__[12*i+5] += (s152*s286*s665)/4. + (s152*s287*s3*s665)/4. + (s10*s152*s322*s665)/4. + (s11*s152*s323*s665)/4. + (s152*s2*s324*s665)/4. + (s12*s152*s329*s665)/4. + (s14*s152*s330*s665)/4. + (s152*s16*s331*s665)/4. + (s152*s326*s4*s665)/4. + (s152*s327*s6*s665)/4. + s325*(s333 + (s0*s152*s665)/4.) + s334*((s152*s347*s665)/4. + (s151*(s199*s673 + s188*s676))/2.) + s348*((s152*s363*s665)/4. + (s151*(s251 + s252 + s201*s673 + s191*s676))/2.) + s364*((s152*s367*s665)/4. + (s151*(s203*s673 + s194*s676))/2.) + s368*((s152*s371*s665)/4. + (s151*(s205*s673 + s197*s676))/2.) + s372*((s152*s376*s665)/4. + (s151*(s267 + s268 + s201*s680 + s191*s684))/2.) + s377*((s152*s388*s665)/4. + (s151*(s203*s680 + s194*s684))/2.) + s389*((s152*s392*s665)/4. + (s151*(s205*s680 + s197*s684))/2.) + s393*((s152*s396*s665)/4. + (s151*(s203*s689 + s194*s693))/2.) + s397*((s152*s410*s665)/4. + (s151*(s205*s689 + s197*s693))/2.) + s411*((s152*s420*s665)/4. + (s151*(s205*(s1008 + s1009 - s205*s210*s225*s668 - s197*s223*s225*s668) + s197*(-(s197*s215*s225*s668) - s205*s223*s225*s668 + s795)))/2.) + (s152*s332*s665*s8)/4. + (s152*s328*s665*s9)/4.;
				buffer__[12*i+6] += (s152*s286*s724)/4. + (s152*s287*s3*s724)/4. + (s10*s152*s322*s724)/4. + (s11*s152*s323*s724)/4. + (s152*s2*s324*s724)/4. + (s0*s152*s325*s724)/4. + (s12*s152*s329*s724)/4. + (s14*s152*s330*s724)/4. + (s152*s16*s331*s724)/4. + (s152*s327*s6*s724)/4. + s326*(s333 + (s152*s4*s724)/4.) + s334*((s152*s347*s724)/4. + (s151*(s199*s732 + s188*s735))/2.) + s348*((s152*s363*s724)/4. + (s151*(s201*s732 + s191*s735))/2.) + s364*((s152*s367*s724)/4. + (s151*(s251 + s252 + s203*s732 + s194*s735))/2.) + s368*((s152*s371*s724)/4. + (s151*(s205*s732 + s197*s735))/2.) + s372*((s152*s376*s724)/4. + (s151*(s201*s740 + s191*s743))/2.) + s377*((s152*s388*s724)/4. + (s151*(s267 + s268 + s203*s740 + s194*s743))/2.) + s389*((s152*s392*s724)/4. + (s151*(s205*s740 + s197*s743))/2.) + s393*((s152*s396*s724)/4. + (s151*(s283 + s284 + s203*s747 + s194*s751))/2.) + s397*((s152*s410*s724)/4. + (s151*(s205*s747 + s197*s751))/2.) + (s152*s332*s724*s8)/4. + s411*((s152*s420*s724)/4. + (s151*(s205*(s1017 + s1018 - s205*s210*s225*s727 - s197*s223*s225*s727) + s197*(-(s197*s215*s225*s727) - s205*s223*s225*s727 + s804)))/2.) + (s152*s328*s724*s9)/4.;
				buffer__[12*i+7] += (s152*s286*s782)/4. + (s152*s287*s3*s782)/4. + (s10*s152*s322*s782)/4. + (s11*s152*s323*s782)/4. + (s152*s2*s324*s782)/4. + (s0*s152*s325*s782)/4. + (s12*s152*s329*s782)/4. + (s14*s152*s330*s782)/4. + (s152*s16*s331*s782)/4. + (s152*s326*s4*s782)/4. + s327*(s333 + (s152*s6*s782)/4.) + s411*((s152*s420*s782)/4. + (s151*(s416 + s417 + s205*(s1025 + s610 - s205*s210*s225*s785 - s197*s223*s225*s785) + s197*(-(s214*s218) + s615 - s197*s215*s225*s785 - s205*s223*s225*s785)))/2.) + s334*((s152*s347*s782)/4. + (s151*(s199*s791 + s188*s794))/2.) + s348*((s152*s363*s782)/4. + (s151*(s201*s791 + s191*s794))/2.) + s364*((s152*s367*s782)/4. + (s151*(s203*s791 + s194*s794))/2.) + s368*((s152*s371*s782)/4. + (s151*(s251 + s252 + s205*s791 + s197*s794))/2.) + (s152*s332*s782*s8)/4. + s372*((s152*s376*s782)/4. + (s151*(s201*s800 + s191*s803))/2.) + s377*((s152*s388*s782)/4. + (s151*(s203*s800 + s194*s803))/2.) + s389*((s152*s392*s782)/4. + (s151*(s267 + s268 + s205*s800 + s197*s803))/2.) + s393*((s152*s396*s782)/4. + (s151*(s203*s809 + s194*s812))/2.) + s397*((s152*s410*s782)/4. + (s151*(s283 + s284 + s205*s809 + s197*s812))/2.) + (s152*s328*s782*s9)/4.;
				buffer__[12*i+8] += (s152*s286*s840)/4. + (s152*s287*s3*s840)/4. + (s10*s152*s322*s840)/4. + (s11*s152*s323*s840)/4. + (s152*s2*s324*s840)/4. + (s0*s152*s325*s840)/4. + (s12*s152*s329*s840)/4. + (s14*s152*s330*s840)/4. + (s152*s16*s331*s840)/4. + (s152*s326*s4*s840)/4. + (s152*s327*s6*s840)/4. + (s152*s332*s8*s840)/4. + s334*((s152*s347*s840)/4. + (s151*(s248 + s249 + s199*s848 + s188*s851))/2.) + s348*((s152*s363*s840)/4. + (s151*(s201*s848 + s191*s851))/2.) + s364*((s152*s367*s840)/4. + (s151*(s203*s848 + s194*s851))/2.) + s368*((s152*s371*s840)/4. + (s151*(s205*s848 + s197*s851))/2.) + s372*((s152*s376*s840)/4. + (s151*(s201*s855 + s191*s858))/2.) + s377*((s152*s388*s840)/4. + (s151*(s203*s855 + s194*s858))/2.) + s389*((s152*s392*s840)/4. + (s151*(s205*s855 + s197*s858))/2.) + s393*((s152*s396*s840)/4. + (s151*(s203*s862 + s194*s865))/2.) + s397*((s152*s410*s840)/4. + (s151*(s205*s862 + s197*s865))/2.) + s328*(s333 + (s152*s840*s9)/4.) + s411*((s152*s420*s840)/4. + (s151*(s197*(s789 + s790 - s197*s215*s225*s843 - s205*s223*s225*s843) + s205*(-(s205*s210*s225*s843) - s197*s223*s225*s843 + s998)))/2.);
				buffer__[12*i+9] += (s152*s286*s893)/4. + (s152*s287*s3*s893)/4. + (s10*s152*s322*s893)/4. + (s11*s152*s323*s893)/4. + (s152*s2*s324*s893)/4. + (s0*s152*s325*s893)/4. + (s14*s152*s330*s893)/4. + (s152*s16*s331*s893)/4. + (s152*s326*s4*s893)/4. + (s152*s327*s6*s893)/4. + (s152*s332*s8*s893)/4. + s329*(s333 + (s12*s152*s893)/4.) + s411*((s152*s420*s893)/4. + (s151*(s205*(s1007 - s205*s210*s225*s896 - s197*s223*s225*s896) + s197*(s798 + s799 - s197*s215*s225*s896 - s205*s223*s225*s896)))/2.) + (s152*s328*s893*s9)/4. + s334*((s152*s347*s893)/4. + (s151*(s199*s899 + s188*s902))/2.) + s348*((s152*s363*s893)/4. + (s151*(s248 + s249 + s201*s899 + s191*s902))/2.) + s364*((s152*s367*s893)/4. + (s151*(s203*s899 + s194*s902))/2.) + s368*((s152*s371*s893)/4. + (s151*(s205*s899 + s197*s902))/2.) + s372*((s152*s376*s893)/4. + (s151*(s264 + s265 + s201*s906 + s191*s909))/2.) + s377*((s152*s388*s893)/4. + (s151*(s203*s906 + s194*s909))/2.) + s389*((s152*s392*s893)/4. + (s151*(s205*s906 + s197*s909))/2.) + s393*((s152*s396*s893)/4. + (s151*(s203*s913 + s194*s916))/2.) + s397*((s152*s410*s893)/4. + (s151*(s205*s913 + s197*s916))/2.);
				buffer__[12*i+10] += (s152*s286*s944)/4. + (s152*s287*s3*s944)/4. + (s10*s152*s322*s944)/4. + (s11*s152*s323*s944)/4. + (s152*s2*s324*s944)/4. + (s0*s152*s325*s944)/4. + (s12*s152*s329*s944)/4. + (s152*s16*s331*s944)/4. + (s152*s326*s4*s944)/4. + (s152*s327*s6*s944)/4. + (s152*s332*s8*s944)/4. + (s152*s328*s9*s944)/4. + s330*(s333 + (s14*s152*s944)/4.) + s411*((s152*s420*s944)/4. + (s151*(s205*(s1016 - s205*s210*s225*s947 - s197*s223*s225*s947) + s197*(s807 + s808 - s197*s215*s225*s947 - s205*s223*s225*s947)))/2.) + s334*((s152*s347*s944)/4. + (s151*(s199*s950 + s188*s953))/2.) + s348*((s152*s363*s944)/4. + (s151*(s201*s950 + s191*s953))/2.) + s364*((s152*s367*s944)/4. + (s151*(s248 + s249 + s203*s950 + s194*s953))/2.) + s368*((s152*s371*s944)/4. + (s151*(s205*s950 + s197*s953))/2.) + s372*((s152*s376*s944)/4. + (s151*(s201*s956 + s191*s959))/2.) + s377*((s152*s388*s944)/4. + (s151*(s264 + s265 + s203*s956 + s194*s959))/2.) + s389*((s152*s392*s944)/4. + (s151*(s205*s956 + s197*s959))/2.) + s393*((s152*s396*s944)/4. + (s151*(s280 + s281 + s203*s963 + s194*s966))/2.) + s397*((s152*s410*s944)/4. + (s151*(s205*s963 + s197*s966))/2.);
				buffer__[12*i+11] += (s152*s286*s994)/4. + (s152*s287*s3*s994)/4. + (s10*s152*s322*s994)/4. + (s11*s152*s323*s994)/4. + (s152*s2*s324*s994)/4. + (s0*s152*s325*s994)/4. + (s12*s152*s329*s994)/4. + (s14*s152*s330*s994)/4. + (s152*s326*s4*s994)/4. + (s152*s327*s6*s994)/4. + (s152*s332*s8*s994)/4. + (s152*s328*s9*s994)/4. + s331*(s333 + (s152*s16*s994)/4.) + s334*((s151*(s1006*s188 + s1003*s199))/2. + (s152*s347*s994)/4.) + s348*((s151*(s1006*s191 + s1003*s201))/2. + (s152*s363*s994)/4.) + s364*((s151*(s1006*s194 + s1003*s203))/2. + (s152*s367*s994)/4.) + s368*((s151*(s1006*s197 + s1003*s205 + s248 + s249))/2. + (s152*s371*s994)/4.) + s372*((s151*(s1015*s191 + s1012*s201))/2. + (s152*s376*s994)/4.) + s377*((s151*(s1015*s194 + s1012*s203))/2. + (s152*s388*s994)/4.) + s389*((s151*(s1015*s197 + s1012*s205 + s264 + s265))/2. + (s152*s392*s994)/4.) + s393*((s151*(s1024*s194 + s1021*s203))/2. + (s152*s396*s994)/4.) + s397*((s151*(s1024*s197 + s1021*s205 + s280 + s281))/2. + (s152*s410*s994)/4.) + s411*((s152*s420*s994)/4. + (s151*(s412 + s413 + s205*(-(s198*s218) + s847 - s205*s210*s225*s997 - s197*s223*s225*s997) + s197*(s1025 + s610 - s197*s215*s225*s997 - s205*s223*s225*s997)))/2.);
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[4*simplices__[3*i+1]+1];
				const Real s1 = s0*s0;
				const Real s2 = V_coords__[4*simplices__[3*i+1]+0];
				const Real s3 = V_coords__[4*simplices__[3*i+0]+0];
				const Real s4 = V_coords__[4*simplices__[3*i+1]+2];
				const Real s5 = s4*s4;
				const Real s6 = V_coords__[4*simplices__[3*i+1]+3];
				const Real s7 = s6*s6;
				const Real s8 = V_coords__[4*simplices__[3*i+0]+1];
				const Real s9 = V_coords__[4*simplices__[3*i+2]+0];
				const Real s10 = V_coords__[4*simplices__[3*i+0]+2];
				const Real s11 = V_coords__[4*simplices__[3*i+0]+3];
				const Real s12 = V_coords__[4*simplices__[3*i+2]+1];
				const Real s13 = s12*s12;
				const Real s14 = V_coords__[4*simplices__[3*i+2]+2];
				const Real s15 = s14*s14;
				const Real s16 = V_coords__[4*simplices__[3*i+2]+3];
				const Real s17 = s16*s16;
				const Real s18 = s8*s8;
				const Real s19 = s2*s2;
				const Real s20 = s18*s19;
				const Real s21 = s10*s10;
				const Real s22 = s19*s21;
				const Real s23 = s11*s11;
				const Real s24 = s19*s23;
				const Real s25 = -2*s0*s2*s3*s8;
				const Real s26 = s3*s3;
				const Real s27 = s1*s26;
				const Real s28 = s1*s21;
				const Real s29 = s1*s23;
				const Real s30 = -2*s10*s2*s3*s4;
				const Real s31 = -2*s0*s10*s4*s8;
				const Real s32 = s26*s5;
				const Real s33 = s18*s5;
				const Real s34 = s23*s5;
				const Real s35 = -2*s11*s2*s3*s6;
				const Real s36 = -2*s0*s11*s6*s8;
				const Real s37 = -2*s10*s11*s4*s6;
				const Real s38 = s26*s7;
				const Real s39 = s18*s7;
				const Real s40 = s21*s7;
				const Real s41 = -2*s18*s2*s9;
				const Real s42 = -2*s2*s21*s9;
				const Real s43 = -2*s2*s23*s9;
				const Real s44 = 2*s0*s3*s8*s9;
				const Real s45 = 2*s0*s2*s8*s9;
				const Real s46 = -2*s1*s3*s9;
				const Real s47 = 2*s10*s3*s4*s9;
				const Real s48 = 2*s10*s2*s4*s9;
				const Real s49 = -2*s3*s5*s9;
				const Real s50 = 2*s11*s3*s6*s9;
				const Real s51 = 2*s11*s2*s6*s9;
				const Real s52 = -2*s3*s7*s9;
				const Real s53 = s9*s9;
				const Real s54 = s18*s53;
				const Real s55 = s21*s53;
				const Real s56 = s23*s53;
				const Real s57 = -2*s0*s53*s8;
				const Real s58 = s1*s53;
				const Real s59 = -2*s10*s4*s53;
				const Real s60 = s5*s53;
				const Real s61 = -2*s11*s53*s6;
				const Real s62 = s53*s7;
				const Real s63 = 2*s12*s2*s3*s8;
				const Real s64 = -2*s12*s19*s8;
				const Real s65 = -2*s0*s12*s26;
				const Real s66 = -2*s0*s12*s21;
				const Real s67 = -2*s0*s12*s23;
				const Real s68 = 2*s0*s12*s2*s3;
				const Real s69 = 2*s10*s12*s4*s8;
				const Real s70 = 2*s0*s10*s12*s4;
				const Real s71 = -2*s12*s5*s8;
				const Real s72 = 2*s11*s12*s6*s8;
				const Real s73 = 2*s0*s11*s12*s6;
				const Real s74 = -2*s12*s7*s8;
				const Real s75 = -2*s12*s3*s8*s9;
				const Real s76 = 2*s12*s2*s8*s9;
				const Real s77 = 2*s0*s12*s3*s9;
				const Real s78 = -2*s0*s12*s2*s9;
				const Real s79 = s13*s26;
				const Real s80 = s13*s21;
				const Real s81 = s13*s23;
				const Real s82 = -2*s13*s2*s3;
				const Real s83 = s13*s19;
				const Real s84 = -2*s10*s13*s4;
				const Real s85 = s13*s5;
				const Real s86 = -2*s11*s13*s6;
				const Real s87 = s13*s7;
				const Real s88 = 2*s10*s14*s2*s3;
				const Real s89 = -2*s10*s14*s19;
				const Real s90 = 2*s0*s10*s14*s8;
				const Real s91 = -2*s1*s10*s14;
				const Real s92 = -2*s14*s26*s4;
				const Real s93 = -2*s14*s18*s4;
				const Real s94 = -2*s14*s23*s4;
				const Real s95 = 2*s14*s2*s3*s4;
				const Real s96 = 2*s0*s14*s4*s8;
				const Real s97 = 2*s10*s11*s14*s6;
				const Real s98 = 2*s11*s14*s4*s6;
				const Real s99 = -2*s10*s14*s7;
				const Real s100 = -2*s10*s14*s3*s9;
				const Real s101 = 2*s10*s14*s2*s9;
				const Real s102 = 2*s14*s3*s4*s9;
				const Real s103 = -2*s14*s2*s4*s9;
				const Real s104 = -2*s10*s12*s14*s8;
				const Real s105 = 2*s0*s10*s12*s14;
				const Real s106 = 2*s12*s14*s4*s8;
				const Real s107 = -2*s0*s12*s14*s4;
				const Real s108 = s15*s26;
				const Real s109 = s15*s18;
				const Real s110 = s15*s23;
				const Real s111 = -2*s15*s2*s3;
				const Real s112 = s15*s19;
				const Real s113 = -2*s0*s15*s8;
				const Real s114 = s1*s15;
				const Real s115 = -2*s11*s15*s6;
				const Real s116 = s15*s7;
				const Real s117 = 2*s11*s16*s2*s3;
				const Real s118 = -2*s11*s16*s19;
				const Real s119 = 2*s0*s11*s16*s8;
				const Real s120 = -2*s1*s11*s16;
				const Real s121 = 2*s10*s11*s16*s4;
				const Real s122 = -2*s11*s16*s5;
				const Real s123 = -2*s16*s26*s6;
				const Real s124 = -2*s16*s18*s6;
				const Real s125 = -2*s16*s21*s6;
				const Real s126 = 2*s16*s2*s3*s6;
				const Real s127 = 2*s0*s16*s6*s8;
				const Real s128 = 2*s10*s16*s4*s6;
				const Real s129 = -2*s11*s16*s3*s9;
				const Real s130 = 2*s11*s16*s2*s9;
				const Real s131 = 2*s16*s3*s6*s9;
				const Real s132 = -2*s16*s2*s6*s9;
				const Real s133 = -2*s11*s12*s16*s8;
				const Real s134 = 2*s0*s11*s12*s16;
				const Real s135 = 2*s12*s16*s6*s8;
				const Real s136 = -2*s0*s12*s16*s6;
				const Real s137 = -2*s10*s11*s14*s16;
				const Real s138 = 2*s11*s14*s16*s4;
				const Real s139 = 2*s10*s14*s16*s6;
				const Real s140 = -2*s14*s16*s4*s6;
				const Real s141 = s17*s26;
				const Real s142 = s17*s18;
				const Real s143 = s17*s21;
				const Real s144 = -2*s17*s2*s3;
				const Real s145 = s17*s19;
				const Real s146 = -2*s0*s17*s8;
				const Real s147 = s1*s17;
				const Real s148 = -2*s10*s17*s4;
				const Real s149 = s17*s5;
				const Real s150 = s100 + s101 + s102 + s103 + s104 + s105 + s106 + s107 + s108 + s109 + s110 + s111 + s112 + s113 + s114 + s115 + s116 + s117 + s118 + s119 + s120 + s121 + s122 + s123 + s124 + s125 + s126 + s127 + s128 + s129 + s130 + s131 + s132 + s133 + s134 + s135 + s136 + s137 + s138 + s139 + s140 + s141 + s142 + s143 + s144 + s145 + s146 + s147 + s148 + s149 + s20 + s22 + s24 + s25 + s27 + s28 + s29 + s30 + s31 + s32 + s33 + s34 + s35 + s36 + s37 + s38 + s39 + s40 + s41 + s42 + s43 + s44 + s45 + s46 + s47 + s48 + s49 + s50 + s51 + s52 + s54 + s55 + s56 + s57 + s58 + s59 + s60 + s61 + s62 + s63 + s64 + s65 + s66 + s67 + s68 + s69 + s70 + s71 + s72 + s73 + s74 + s75 + s76 + s77 + s78 + s79 + s80 + s81 + s82 + s83 + s84 + s85 + s86 + s87 + s88 + s89 + s90 + s91 + s92 + s93 + s94 + s95 + s96 + s97 + s98 + s99;
				const Real s151 = sqrt(s150);
				const Real s152 = 1/s151;
				const Real s153 = -2*s0*s2*s8;
				const Real s154 = 2*s1*s3;
				const Real s155 = -2*s10*s2*s4;
				const Real s156 = 2*s3*s5;
				const Real s157 = -2*s11*s2*s6;
				const Real s158 = 2*s3*s7;
				const Real s159 = 2*s0*s8*s9;
				const Real s160 = -2*s1*s9;
				const Real s161 = 2*s10*s4*s9;
				const Real s162 = -2*s5*s9;
				const Real s163 = 2*s11*s6*s9;
				const Real s164 = -2*s7*s9;
				const Real s165 = 2*s12*s2*s8;
				const Real s166 = -4*s0*s12*s3;
				const Real s167 = 2*s0*s12*s2;
				const Real s168 = -2*s12*s8*s9;
				const Real s169 = 2*s0*s12*s9;
				const Real s170 = 2*s13*s3;
				const Real s171 = -2*s13*s2;
				const Real s172 = 2*s10*s14*s2;
				const Real s173 = -4*s14*s3*s4;
				const Real s174 = 2*s14*s2*s4;
				const Real s175 = -2*s10*s14*s9;
				const Real s176 = 2*s14*s4*s9;
				const Real s177 = 2*s15*s3;
				const Real s178 = -2*s15*s2;
				const Real s179 = 2*s11*s16*s2;
				const Real s180 = -4*s16*s3*s6;
				const Real s181 = 2*s16*s2*s6;
				const Real s182 = -2*s11*s16*s9;
				const Real s183 = 2*s16*s6*s9;
				const Real s184 = 2*s17*s3;
				const Real s185 = -2*s17*s2;
				const Real s186 = s153 + s154 + s155 + s156 + s157 + s158 + s159 + s160 + s161 + s162 + s163 + s164 + s165 + s166 + s167 + s168 + s169 + s170 + s171 + s172 + s173 + s174 + s175 + s176 + s177 + s178 + s179 + s180 + s181 + s182 + s183 + s184 + s185;
				const Real s187 = -s3;
				const Real s188 = s187 + s2;
				const Real s189 = s188*s188;
				const Real s190 = -s8;
				const Real s191 = s0 + s190;
				const Real s192 = s191*s191;
				const Real s193 = -s10;
				const Real s194 = s193 + s4;
				const Real s195 = s194*s194;
				const Real s196 = -s11;
				const Real s197 = s196 + s6;
				const Real s198 = s197*s197;
				const Real s199 = s187 + s9;
				const Real s200 = s188*s199;
				const Real s201 = s12 + s190;
				const Real s202 = s191*s201;
				const Real s203 = s14 + s193;
				const Real s204 = s194*s203;
				const Real s205 = s16 + s196;
				const Real s206 = s197*s205;
				const Real s207 = s200 + s202 + s204 + s206;
				const Real s208 = s207*s207;
				const Real s209 = -s208;
				const Real s210 = s189 + s192 + s195 + s198;
				const Real s211 = s199*s199;
				const Real s212 = s201*s201;
				const Real s213 = s203*s203;
				const Real s214 = s205*s205;
				const Real s215 = s211 + s212 + s213 + s214;
				const Real s216 = s210*s215;
				const Real s217 = s209 + s216;
				const Real s218 = 1/s217;
				const Real s219 = -(s188*s199);
				const Real s220 = -(s191*s201);
				const Real s221 = -(s194*s203);
				const Real s222 = -(s197*s205);
				const Real s223 = s219 + s220 + s221 + s222;
				const Real s224 = s217*s217;
				const Real s225 = 1/s224;
				const Real s226 = -2*s199*s210;
				const Real s227 = 2*s3;
				const Real s228 = -s2;
				const Real s229 = -s9;
				const Real s230 = s227 + s228 + s229;
				const Real s231 = -2*s207*s230;
				const Real s232 = -2*s188*s215;
				const Real s233 = s226 + s231 + s232;
				const Real s234 = -2*s188*s199*s218;
				const Real s235 = -2*s3;
				const Real s236 = s2 + s235 + s9;
				const Real s237 = -(s218*s223);
				const Real s238 = -(s199*s210*s225*s233);
				const Real s239 = -(s188*s223*s225*s233);
				const Real s240 = -(s210*s218);
				const Real s241 = s188*s218*s236;
				const Real s242 = s234 + s237 + s238 + s239 + s240 + s241;
				const Real s243 = -(s199*s223*s225*s233);
				const Real s244 = -(s188*s215*s225*s233);
				const Real s245 = s199*s218*s236;
				const Real s246 = -(s215*s218);
				const Real s247 = s234 + s237 + s243 + s244 + s245 + s246;
				const Real s248 = s199*s210*s218;
				const Real s249 = s188*s218*s223;
				const Real s250 = s248 + s249;
				const Real s251 = s199*s218*s223;
				const Real s252 = s188*s215*s218;
				const Real s253 = s251 + s252;
				const Real s254 = -(s201*s210*s225*s233);
				const Real s255 = -(s191*s223*s225*s233);
				const Real s256 = s191*s218*s236;
				const Real s257 = -2*s188*s201*s218;
				const Real s258 = s254 + s255 + s256 + s257;
				const Real s259 = -(s201*s223*s225*s233);
				const Real s260 = -(s191*s215*s225*s233);
				const Real s261 = -2*s191*s199*s218;
				const Real s262 = s201*s218*s236;
				const Real s263 = s259 + s260 + s261 + s262;
				const Real s264 = s201*s210*s218;
				const Real s265 = s191*s218*s223;
				const Real s266 = s264 + s265;
				const Real s267 = s201*s218*s223;
				const Real s268 = s191*s215*s218;
				const Real s269 = s267 + s268;
				const Real s270 = -(s203*s210*s225*s233);
				const Real s271 = -(s194*s223*s225*s233);
				const Real s272 = s194*s218*s236;
				const Real s273 = -2*s188*s203*s218;
				const Real s274 = s270 + s271 + s272 + s273;
				const Real s275 = -(s203*s223*s225*s233);
				const Real s276 = -(s194*s215*s225*s233);
				const Real s277 = -2*s194*s199*s218;
				const Real s278 = s203*s218*s236;
				const Real s279 = s275 + s276 + s277 + s278;
				const Real s280 = s203*s210*s218;
				const Real s281 = s194*s218*s223;
				const Real s282 = s280 + s281;
				const Real s283 = s203*s218*s223;
				const Real s284 = s194*s215*s218;
				const Real s285 = s283 + s284;
				const Real s286 = P_D_near__[23*i+0];
				const Real s287 = P_D_near__[23*i+1];
				const Real s288 = 2*s19*s8;
				const Real s289 = -2*s0*s2*s3;
				const Real s290 = -2*s0*s10*s4;
				const Real s291 = 2*s5*s8;
				const Real s292 = -2*s0*s11*s6;
				const Real s293 = 2*s7*s8;
				const Real s294 = -4*s2*s8*s9;
				const Real s295 = 2*s0*s3*s9;
				const Real s296 = 2*s0*s2*s9;
				const Real s297 = 2*s53*s8;
				const Real s298 = -2*s0*s53;
				const Real s299 = 2*s12*s2*s3;
				const Real s300 = -2*s12*s19;
				const Real s301 = 2*s10*s12*s4;
				const Real s302 = -2*s12*s5;
				const Real s303 = 2*s11*s12*s6;
				const Real s304 = -2*s12*s7;
				const Real s305 = -2*s12*s3*s9;
				const Real s306 = 2*s12*s2*s9;
				const Real s307 = 2*s0*s10*s14;
				const Real s308 = -4*s14*s4*s8;
				const Real s309 = 2*s0*s14*s4;
				const Real s310 = -2*s10*s12*s14;
				const Real s311 = 2*s12*s14*s4;
				const Real s312 = 2*s15*s8;
				const Real s313 = -2*s0*s15;
				const Real s314 = 2*s0*s11*s16;
				const Real s315 = -4*s16*s6*s8;
				const Real s316 = 2*s0*s16*s6;
				const Real s317 = -2*s11*s12*s16;
				const Real s318 = 2*s12*s16*s6;
				const Real s319 = 2*s17*s8;
				const Real s320 = -2*s0*s17;
				const Real s321 = s288 + s289 + s290 + s291 + s292 + s293 + s294 + s295 + s296 + s297 + s298 + s299 + s300 + s301 + s302 + s303 + s304 + s305 + s306 + s307 + s308 + s309 + s310 + s311 + s312 + s313 + s314 + s315 + s316 + s317 + s318 + s319 + s320;
				const Real s322 = P_D_near__[23*i+3];
				const Real s323 = P_D_near__[23*i+4];
				const Real s324 = P_D_near__[23*i+5];
				const Real s325 = P_D_near__[23*i+6];
				const Real s326 = P_D_near__[23*i+7];
				const Real s327 = P_D_near__[23*i+8];
				const Real s328 = P_D_near__[23*i+9];
				const Real s329 = P_D_near__[23*i+10];
				const Real s330 = P_D_near__[23*i+11];
				const Real s331 = P_D_near__[23*i+12];
				const Real s332 = P_D_near__[23*i+2];
				const Real s333 = s151/2.;
				const Real s334 = P_D_near__[23*i+13];
				const Real s335 = -2*s201*s210;
				const Real s336 = 2*s8;
				const Real s337 = -s0;
				const Real s338 = -s12;
				const Real s339 = s336 + s337 + s338;
				const Real s340 = -2*s207*s339;
				const Real s341 = -2*s191*s215;
				const Real s342 = s335 + s340 + s341;
				const Real s343 = -2*s8;
				const Real s344 = s0 + s12 + s343;
				const Real s345 = s199*s250;
				const Real s346 = s188*s253;
				const Real s347 = s345 + s346;
				const Real s348 = P_D_near__[23*i+14];
				const Real s349 = -(s199*s210*s218);
				const Real s350 = -(s188*s218*s223);
				const Real s351 = -(s199*s218*s223);
				const Real s352 = -(s188*s215*s218);
				const Real s353 = -(s199*s210*s225*s342);
				const Real s354 = -(s188*s223*s225*s342);
				const Real s355 = s188*s218*s344;
				const Real s356 = s261 + s353 + s354 + s355;
				const Real s357 = -(s199*s223*s225*s342);
				const Real s358 = -(s188*s215*s225*s342);
				const Real s359 = s199*s218*s344;
				const Real s360 = s257 + s357 + s358 + s359;
				const Real s361 = s201*s250;
				const Real s362 = s191*s253;
				const Real s363 = s361 + s362;
				const Real s364 = P_D_near__[23*i+15];
				const Real s365 = s203*s250;
				const Real s366 = s194*s253;
				const Real s367 = s365 + s366;
				const Real s368 = P_D_near__[23*i+16];
				const Real s369 = s205*s250;
				const Real s370 = s197*s253;
				const Real s371 = s369 + s370;
				const Real s372 = P_D_near__[23*i+17];
				const Real s373 = -2*s191*s201*s218;
				const Real s374 = s201*s266;
				const Real s375 = s191*s269;
				const Real s376 = s374 + s375;
				const Real s377 = P_D_near__[23*i+18];
				const Real s378 = -(s201*s210*s225*s342);
				const Real s379 = -(s191*s223*s225*s342);
				const Real s380 = s191*s218*s344;
				const Real s381 = s237 + s240 + s373 + s378 + s379 + s380;
				const Real s382 = -(s201*s223*s225*s342);
				const Real s383 = -(s191*s215*s225*s342);
				const Real s384 = s201*s218*s344;
				const Real s385 = s237 + s246 + s373 + s382 + s383 + s384;
				const Real s386 = s203*s266;
				const Real s387 = s194*s269;
				const Real s388 = s386 + s387;
				const Real s389 = P_D_near__[23*i+19];
				const Real s390 = s205*s266;
				const Real s391 = s197*s269;
				const Real s392 = s390 + s391;
				const Real s393 = P_D_near__[23*i+20];
				const Real s394 = s203*s282;
				const Real s395 = s194*s285;
				const Real s396 = s394 + s395;
				const Real s397 = P_D_near__[23*i+21];
				const Real s398 = -(s203*s210*s225*s342);
				const Real s399 = -(s194*s223*s225*s342);
				const Real s400 = s194*s218*s344;
				const Real s401 = -2*s191*s203*s218;
				const Real s402 = s398 + s399 + s400 + s401;
				const Real s403 = -(s203*s223*s225*s342);
				const Real s404 = -(s194*s215*s225*s342);
				const Real s405 = -2*s194*s201*s218;
				const Real s406 = s203*s218*s344;
				const Real s407 = s403 + s404 + s405 + s406;
				const Real s408 = s205*s282;
				const Real s409 = s197*s285;
				const Real s410 = s408 + s409;
				const Real s411 = P_D_near__[23*i+22];
				const Real s412 = s205*s210*s218;
				const Real s413 = s197*s218*s223;
				const Real s414 = s412 + s413;
				const Real s415 = s205*s414;
				const Real s416 = s205*s218*s223;
				const Real s417 = s197*s215*s218;
				const Real s418 = s416 + s417;
				const Real s419 = s197*s418;
				const Real s420 = s415 + s419;
				const Real s421 = 2*s10*s19;
				const Real s422 = 2*s1*s10;
				const Real s423 = -2*s2*s3*s4;
				const Real s424 = -2*s0*s4*s8;
				const Real s425 = -2*s11*s4*s6;
				const Real s426 = 2*s10*s7;
				const Real s427 = -4*s10*s2*s9;
				const Real s428 = 2*s3*s4*s9;
				const Real s429 = 2*s2*s4*s9;
				const Real s430 = 2*s10*s53;
				const Real s431 = -2*s4*s53;
				const Real s432 = -4*s0*s10*s12;
				const Real s433 = 2*s12*s4*s8;
				const Real s434 = 2*s0*s12*s4;
				const Real s435 = 2*s10*s13;
				const Real s436 = -2*s13*s4;
				const Real s437 = 2*s14*s2*s3;
				const Real s438 = -2*s14*s19;
				const Real s439 = 2*s0*s14*s8;
				const Real s440 = -2*s1*s14;
				const Real s441 = 2*s11*s14*s6;
				const Real s442 = -2*s14*s7;
				const Real s443 = -2*s14*s3*s9;
				const Real s444 = 2*s14*s2*s9;
				const Real s445 = -2*s12*s14*s8;
				const Real s446 = 2*s0*s12*s14;
				const Real s447 = 2*s11*s16*s4;
				const Real s448 = -4*s10*s16*s6;
				const Real s449 = 2*s16*s4*s6;
				const Real s450 = -2*s11*s14*s16;
				const Real s451 = 2*s14*s16*s6;
				const Real s452 = 2*s10*s17;
				const Real s453 = -2*s17*s4;
				const Real s454 = s421 + s422 + s423 + s424 + s425 + s426 + s427 + s428 + s429 + s430 + s431 + s432 + s433 + s434 + s435 + s436 + s437 + s438 + s439 + s440 + s441 + s442 + s443 + s444 + s445 + s446 + s447 + s448 + s449 + s450 + s451 + s452 + s453;
				const Real s455 = -2*s203*s210;
				const Real s456 = 2*s10;
				const Real s457 = -s4;
				const Real s458 = -s14;
				const Real s459 = s456 + s457 + s458;
				const Real s460 = -2*s207*s459;
				const Real s461 = -2*s194*s215;
				const Real s462 = s455 + s460 + s461;
				const Real s463 = -2*s10;
				const Real s464 = s14 + s4 + s463;
				const Real s465 = -(s199*s210*s225*s462);
				const Real s466 = -(s188*s223*s225*s462);
				const Real s467 = s188*s218*s464;
				const Real s468 = s277 + s465 + s466 + s467;
				const Real s469 = -(s199*s223*s225*s462);
				const Real s470 = -(s188*s215*s225*s462);
				const Real s471 = s199*s218*s464;
				const Real s472 = s273 + s469 + s470 + s471;
				const Real s473 = -(s201*s210*s218);
				const Real s474 = -(s191*s218*s223);
				const Real s475 = -(s201*s218*s223);
				const Real s476 = -(s191*s215*s218);
				const Real s477 = -(s201*s210*s225*s462);
				const Real s478 = -(s191*s223*s225*s462);
				const Real s479 = s191*s218*s464;
				const Real s480 = s405 + s477 + s478 + s479;
				const Real s481 = -(s201*s223*s225*s462);
				const Real s482 = -(s191*s215*s225*s462);
				const Real s483 = s201*s218*s464;
				const Real s484 = s401 + s481 + s482 + s483;
				const Real s485 = -2*s194*s203*s218;
				const Real s486 = -(s203*s210*s225*s462);
				const Real s487 = -(s194*s223*s225*s462);
				const Real s488 = s194*s218*s464;
				const Real s489 = s237 + s240 + s485 + s486 + s487 + s488;
				const Real s490 = -(s203*s223*s225*s462);
				const Real s491 = -(s194*s215*s225*s462);
				const Real s492 = s203*s218*s464;
				const Real s493 = s237 + s246 + s485 + s490 + s491 + s492;
				const Real s494 = 2*s11*s19;
				const Real s495 = 2*s1*s11;
				const Real s496 = 2*s11*s5;
				const Real s497 = -2*s2*s3*s6;
				const Real s498 = -2*s0*s6*s8;
				const Real s499 = -2*s10*s4*s6;
				const Real s500 = -4*s11*s2*s9;
				const Real s501 = 2*s3*s6*s9;
				const Real s502 = 2*s2*s6*s9;
				const Real s503 = 2*s11*s53;
				const Real s504 = -2*s53*s6;
				const Real s505 = -4*s0*s11*s12;
				const Real s506 = 2*s12*s6*s8;
				const Real s507 = 2*s0*s12*s6;
				const Real s508 = 2*s11*s13;
				const Real s509 = -2*s13*s6;
				const Real s510 = -4*s11*s14*s4;
				const Real s511 = 2*s10*s14*s6;
				const Real s512 = 2*s14*s4*s6;
				const Real s513 = 2*s11*s15;
				const Real s514 = -2*s15*s6;
				const Real s515 = 2*s16*s2*s3;
				const Real s516 = -2*s16*s19;
				const Real s517 = 2*s0*s16*s8;
				const Real s518 = -2*s1*s16;
				const Real s519 = 2*s10*s16*s4;
				const Real s520 = -2*s16*s5;
				const Real s521 = -2*s16*s3*s9;
				const Real s522 = 2*s16*s2*s9;
				const Real s523 = -2*s12*s16*s8;
				const Real s524 = 2*s0*s12*s16;
				const Real s525 = -2*s10*s14*s16;
				const Real s526 = 2*s14*s16*s4;
				const Real s527 = s494 + s495 + s496 + s497 + s498 + s499 + s500 + s501 + s502 + s503 + s504 + s505 + s506 + s507 + s508 + s509 + s510 + s511 + s512 + s513 + s514 + s515 + s516 + s517 + s518 + s519 + s520 + s521 + s522 + s523 + s524 + s525 + s526;
				const Real s528 = -2*s205*s210;
				const Real s529 = 2*s11;
				const Real s530 = -s6;
				const Real s531 = -s16;
				const Real s532 = s529 + s530 + s531;
				const Real s533 = -2*s207*s532;
				const Real s534 = -2*s197*s215;
				const Real s535 = s528 + s533 + s534;
				const Real s536 = -2*s197*s199*s218;
				const Real s537 = -2*s188*s205*s218;
				const Real s538 = -2*s11;
				const Real s539 = s16 + s538 + s6;
				const Real s540 = -(s199*s210*s225*s535);
				const Real s541 = -(s188*s223*s225*s535);
				const Real s542 = s188*s218*s539;
				const Real s543 = s536 + s540 + s541 + s542;
				const Real s544 = -(s199*s223*s225*s535);
				const Real s545 = -(s188*s215*s225*s535);
				const Real s546 = s199*s218*s539;
				const Real s547 = s537 + s544 + s545 + s546;
				const Real s548 = -2*s197*s201*s218;
				const Real s549 = -2*s191*s205*s218;
				const Real s550 = -(s201*s210*s225*s535);
				const Real s551 = -(s191*s223*s225*s535);
				const Real s552 = s191*s218*s539;
				const Real s553 = s548 + s550 + s551 + s552;
				const Real s554 = -(s201*s223*s225*s535);
				const Real s555 = -(s191*s215*s225*s535);
				const Real s556 = s201*s218*s539;
				const Real s557 = s549 + s554 + s555 + s556;
				const Real s558 = -2*s197*s203*s218;
				const Real s559 = -2*s194*s205*s218;
				const Real s560 = -(s203*s210*s218);
				const Real s561 = -(s194*s218*s223);
				const Real s562 = -(s203*s218*s223);
				const Real s563 = -(s194*s215*s218);
				const Real s564 = -(s203*s210*s225*s535);
				const Real s565 = -(s194*s223*s225*s535);
				const Real s566 = s194*s218*s539;
				const Real s567 = s558 + s564 + s565 + s566;
				const Real s568 = -(s203*s223*s225*s535);
				const Real s569 = -(s194*s215*s225*s535);
				const Real s570 = s203*s218*s539;
				const Real s571 = s559 + s568 + s569 + s570;
				const Real s572 = -2*s197*s205*s218;
				const Real s573 = 2*s18*s2;
				const Real s574 = 2*s2*s21;
				const Real s575 = 2*s2*s23;
				const Real s576 = -2*s0*s3*s8;
				const Real s577 = -2*s10*s3*s4;
				const Real s578 = -2*s11*s3*s6;
				const Real s579 = -2*s18*s9;
				const Real s580 = -2*s21*s9;
				const Real s581 = -2*s23*s9;
				const Real s582 = 2*s12*s3*s8;
				const Real s583 = -4*s12*s2*s8;
				const Real s584 = 2*s0*s12*s3;
				const Real s585 = 2*s12*s8*s9;
				const Real s586 = -2*s0*s12*s9;
				const Real s587 = -2*s13*s3;
				const Real s588 = 2*s13*s2;
				const Real s589 = 2*s10*s14*s3;
				const Real s590 = -4*s10*s14*s2;
				const Real s591 = 2*s14*s3*s4;
				const Real s592 = 2*s10*s14*s9;
				const Real s593 = -2*s14*s4*s9;
				const Real s594 = -2*s15*s3;
				const Real s595 = 2*s15*s2;
				const Real s596 = 2*s11*s16*s3;
				const Real s597 = -4*s11*s16*s2;
				const Real s598 = 2*s16*s3*s6;
				const Real s599 = 2*s11*s16*s9;
				const Real s600 = -2*s16*s6*s9;
				const Real s601 = -2*s17*s3;
				const Real s602 = 2*s17*s2;
				const Real s603 = s159 + s161 + s163 + s573 + s574 + s575 + s576 + s577 + s578 + s579 + s580 + s581 + s582 + s583 + s584 + s585 + s586 + s587 + s588 + s589 + s590 + s591 + s592 + s593 + s594 + s595 + s596 + s597 + s598 + s599 + s600 + s601 + s602;
				const Real s604 = -2*s199*s207;
				const Real s605 = 2*s188*s215;
				const Real s606 = s604 + s605;
				const Real s607 = -(s199*s210*s225*s606);
				const Real s608 = -(s188*s223*s225*s606);
				const Real s609 = s188*s199*s218;
				const Real s610 = s218*s223;
				const Real s611 = s607 + s608 + s609 + s610;
				const Real s612 = -(s199*s223*s225*s606);
				const Real s613 = -(s188*s215*s225*s606);
				const Real s614 = -(s211*s218);
				const Real s615 = s215*s218;
				const Real s616 = s612 + s613 + s614 + s615;
				const Real s617 = -(s201*s210*s225*s606);
				const Real s618 = -(s191*s223*s225*s606);
				const Real s619 = -(s191*s199*s218);
				const Real s620 = 2*s188*s201*s218;
				const Real s621 = s617 + s618 + s619 + s620;
				const Real s622 = -(s201*s223*s225*s606);
				const Real s623 = -(s191*s215*s225*s606);
				const Real s624 = -(s199*s201*s218);
				const Real s625 = s622 + s623 + s624;
				const Real s626 = -(s203*s210*s225*s606);
				const Real s627 = -(s194*s223*s225*s606);
				const Real s628 = -(s194*s199*s218);
				const Real s629 = 2*s188*s203*s218;
				const Real s630 = s626 + s627 + s628 + s629;
				const Real s631 = -(s203*s223*s225*s606);
				const Real s632 = -(s194*s215*s225*s606);
				const Real s633 = -(s199*s203*s218);
				const Real s634 = s631 + s632 + s633;
				const Real s635 = -2*s2*s3*s8;
				const Real s636 = 2*s0*s26;
				const Real s637 = 2*s0*s21;
				const Real s638 = 2*s0*s23;
				const Real s639 = -2*s10*s4*s8;
				const Real s640 = -2*s11*s6*s8;
				const Real s641 = 2*s3*s8*s9;
				const Real s642 = 2*s2*s8*s9;
				const Real s643 = -4*s0*s3*s9;
				const Real s644 = -2*s53*s8;
				const Real s645 = 2*s0*s53;
				const Real s646 = -2*s12*s26;
				const Real s647 = -2*s12*s21;
				const Real s648 = -2*s12*s23;
				const Real s649 = 2*s12*s3*s9;
				const Real s650 = -2*s12*s2*s9;
				const Real s651 = 2*s10*s14*s8;
				const Real s652 = -4*s0*s10*s14;
				const Real s653 = 2*s14*s4*s8;
				const Real s654 = 2*s10*s12*s14;
				const Real s655 = -2*s12*s14*s4;
				const Real s656 = -2*s15*s8;
				const Real s657 = 2*s0*s15;
				const Real s658 = 2*s11*s16*s8;
				const Real s659 = -4*s0*s11*s16;
				const Real s660 = 2*s16*s6*s8;
				const Real s661 = 2*s11*s12*s16;
				const Real s662 = -2*s12*s16*s6;
				const Real s663 = -2*s17*s8;
				const Real s664 = 2*s0*s17;
				const Real s665 = s299 + s301 + s303 + s635 + s636 + s637 + s638 + s639 + s640 + s641 + s642 + s643 + s644 + s645 + s646 + s647 + s648 + s649 + s650 + s651 + s652 + s653 + s654 + s655 + s656 + s657 + s658 + s659 + s660 + s661 + s662 + s663 + s664;
				const Real s666 = -2*s201*s207;
				const Real s667 = 2*s191*s215;
				const Real s668 = s666 + s667;
				const Real s669 = -(s199*s210*s225*s668);
				const Real s670 = -(s188*s223*s225*s668);
				const Real s671 = 2*s191*s199*s218;
				const Real s672 = -(s188*s201*s218);
				const Real s673 = s669 + s670 + s671 + s672;
				const Real s674 = -(s199*s223*s225*s668);
				const Real s675 = -(s188*s215*s225*s668);
				const Real s676 = s624 + s674 + s675;
				const Real s677 = -(s201*s210*s225*s668);
				const Real s678 = -(s191*s223*s225*s668);
				const Real s679 = s191*s201*s218;
				const Real s680 = s610 + s677 + s678 + s679;
				const Real s681 = -(s201*s223*s225*s668);
				const Real s682 = -(s191*s215*s225*s668);
				const Real s683 = -(s212*s218);
				const Real s684 = s615 + s681 + s682 + s683;
				const Real s685 = -(s203*s210*s225*s668);
				const Real s686 = -(s194*s223*s225*s668);
				const Real s687 = -(s194*s201*s218);
				const Real s688 = 2*s191*s203*s218;
				const Real s689 = s685 + s686 + s687 + s688;
				const Real s690 = -(s203*s223*s225*s668);
				const Real s691 = -(s194*s215*s225*s668);
				const Real s692 = -(s201*s203*s218);
				const Real s693 = s690 + s691 + s692;
				const Real s694 = -2*s10*s2*s3;
				const Real s695 = -2*s0*s10*s8;
				const Real s696 = 2*s26*s4;
				const Real s697 = 2*s18*s4;
				const Real s698 = 2*s23*s4;
				const Real s699 = -2*s10*s11*s6;
				const Real s700 = 2*s10*s3*s9;
				const Real s701 = 2*s10*s2*s9;
				const Real s702 = -4*s3*s4*s9;
				const Real s703 = -2*s10*s53;
				const Real s704 = 2*s4*s53;
				const Real s705 = 2*s10*s12*s8;
				const Real s706 = 2*s0*s10*s12;
				const Real s707 = -4*s12*s4*s8;
				const Real s708 = -2*s10*s13;
				const Real s709 = 2*s13*s4;
				const Real s710 = -2*s14*s26;
				const Real s711 = -2*s14*s18;
				const Real s712 = -2*s14*s23;
				const Real s713 = 2*s14*s3*s9;
				const Real s714 = -2*s14*s2*s9;
				const Real s715 = 2*s12*s14*s8;
				const Real s716 = -2*s0*s12*s14;
				const Real s717 = 2*s10*s11*s16;
				const Real s718 = -4*s11*s16*s4;
				const Real s719 = 2*s10*s16*s6;
				const Real s720 = 2*s11*s14*s16;
				const Real s721 = -2*s14*s16*s6;
				const Real s722 = -2*s10*s17;
				const Real s723 = 2*s17*s4;
				const Real s724 = s437 + s439 + s441 + s694 + s695 + s696 + s697 + s698 + s699 + s700 + s701 + s702 + s703 + s704 + s705 + s706 + s707 + s708 + s709 + s710 + s711 + s712 + s713 + s714 + s715 + s716 + s717 + s718 + s719 + s720 + s721 + s722 + s723;
				const Real s725 = -2*s203*s207;
				const Real s726 = 2*s194*s215;
				const Real s727 = s725 + s726;
				const Real s728 = -(s199*s210*s225*s727);
				const Real s729 = -(s188*s223*s225*s727);
				const Real s730 = 2*s194*s199*s218;
				const Real s731 = -(s188*s203*s218);
				const Real s732 = s728 + s729 + s730 + s731;
				const Real s733 = -(s199*s223*s225*s727);
				const Real s734 = -(s188*s215*s225*s727);
				const Real s735 = s633 + s733 + s734;
				const Real s736 = -(s201*s210*s225*s727);
				const Real s737 = -(s191*s223*s225*s727);
				const Real s738 = 2*s194*s201*s218;
				const Real s739 = -(s191*s203*s218);
				const Real s740 = s736 + s737 + s738 + s739;
				const Real s741 = -(s201*s223*s225*s727);
				const Real s742 = -(s191*s215*s225*s727);
				const Real s743 = s692 + s741 + s742;
				const Real s744 = -(s203*s210*s225*s727);
				const Real s745 = -(s194*s223*s225*s727);
				const Real s746 = s194*s203*s218;
				const Real s747 = s610 + s744 + s745 + s746;
				const Real s748 = -(s203*s223*s225*s727);
				const Real s749 = -(s194*s215*s225*s727);
				const Real s750 = -(s213*s218);
				const Real s751 = s615 + s748 + s749 + s750;
				const Real s752 = -2*s11*s2*s3;
				const Real s753 = -2*s0*s11*s8;
				const Real s754 = -2*s10*s11*s4;
				const Real s755 = 2*s26*s6;
				const Real s756 = 2*s18*s6;
				const Real s757 = 2*s21*s6;
				const Real s758 = 2*s11*s3*s9;
				const Real s759 = 2*s11*s2*s9;
				const Real s760 = -4*s3*s6*s9;
				const Real s761 = -2*s11*s53;
				const Real s762 = 2*s53*s6;
				const Real s763 = 2*s11*s12*s8;
				const Real s764 = 2*s0*s11*s12;
				const Real s765 = -4*s12*s6*s8;
				const Real s766 = -2*s11*s13;
				const Real s767 = 2*s13*s6;
				const Real s768 = 2*s10*s11*s14;
				const Real s769 = 2*s11*s14*s4;
				const Real s770 = -4*s10*s14*s6;
				const Real s771 = -2*s11*s15;
				const Real s772 = 2*s15*s6;
				const Real s773 = -2*s16*s26;
				const Real s774 = -2*s16*s18;
				const Real s775 = -2*s16*s21;
				const Real s776 = 2*s16*s3*s9;
				const Real s777 = -2*s16*s2*s9;
				const Real s778 = 2*s12*s16*s8;
				const Real s779 = -2*s0*s12*s16;
				const Real s780 = 2*s10*s14*s16;
				const Real s781 = -2*s14*s16*s4;
				const Real s782 = s515 + s517 + s519 + s752 + s753 + s754 + s755 + s756 + s757 + s758 + s759 + s760 + s761 + s762 + s763 + s764 + s765 + s766 + s767 + s768 + s769 + s770 + s771 + s772 + s773 + s774 + s775 + s776 + s777 + s778 + s779 + s780 + s781;
				const Real s783 = -2*s205*s207;
				const Real s784 = 2*s197*s215;
				const Real s785 = s783 + s784;
				const Real s786 = -(s199*s205*s218);
				const Real s787 = -(s199*s210*s225*s785);
				const Real s788 = -(s188*s223*s225*s785);
				const Real s789 = 2*s197*s199*s218;
				const Real s790 = -(s188*s205*s218);
				const Real s791 = s787 + s788 + s789 + s790;
				const Real s792 = -(s199*s223*s225*s785);
				const Real s793 = -(s188*s215*s225*s785);
				const Real s794 = s786 + s792 + s793;
				const Real s795 = -(s201*s205*s218);
				const Real s796 = -(s201*s210*s225*s785);
				const Real s797 = -(s191*s223*s225*s785);
				const Real s798 = 2*s197*s201*s218;
				const Real s799 = -(s191*s205*s218);
				const Real s800 = s796 + s797 + s798 + s799;
				const Real s801 = -(s201*s223*s225*s785);
				const Real s802 = -(s191*s215*s225*s785);
				const Real s803 = s795 + s801 + s802;
				const Real s804 = -(s203*s205*s218);
				const Real s805 = -(s203*s210*s225*s785);
				const Real s806 = -(s194*s223*s225*s785);
				const Real s807 = 2*s197*s203*s218;
				const Real s808 = -(s194*s205*s218);
				const Real s809 = s805 + s806 + s807 + s808;
				const Real s810 = -(s203*s223*s225*s785);
				const Real s811 = -(s194*s215*s225*s785);
				const Real s812 = s804 + s810 + s811;
				const Real s813 = -2*s18*s2;
				const Real s814 = -2*s2*s21;
				const Real s815 = -2*s2*s23;
				const Real s816 = 2*s0*s3*s8;
				const Real s817 = 2*s0*s2*s8;
				const Real s818 = -2*s1*s3;
				const Real s819 = 2*s10*s3*s4;
				const Real s820 = 2*s10*s2*s4;
				const Real s821 = -2*s3*s5;
				const Real s822 = 2*s11*s3*s6;
				const Real s823 = 2*s11*s2*s6;
				const Real s824 = -2*s3*s7;
				const Real s825 = 2*s18*s9;
				const Real s826 = 2*s21*s9;
				const Real s827 = 2*s23*s9;
				const Real s828 = -4*s0*s8*s9;
				const Real s829 = 2*s1*s9;
				const Real s830 = -4*s10*s4*s9;
				const Real s831 = 2*s5*s9;
				const Real s832 = -4*s11*s6*s9;
				const Real s833 = 2*s7*s9;
				const Real s834 = -2*s12*s3*s8;
				const Real s835 = -2*s0*s12*s2;
				const Real s836 = -2*s10*s14*s3;
				const Real s837 = -2*s14*s2*s4;
				const Real s838 = -2*s11*s16*s3;
				const Real s839 = -2*s16*s2*s6;
				const Real s840 = s165 + s172 + s179 + s584 + s591 + s598 + s813 + s814 + s815 + s816 + s817 + s818 + s819 + s820 + s821 + s822 + s823 + s824 + s825 + s826 + s827 + s828 + s829 + s830 + s831 + s832 + s833 + s834 + s835 + s836 + s837 + s838 + s839;
				const Real s841 = 2*s199*s210;
				const Real s842 = -2*s188*s207;
				const Real s843 = s841 + s842;
				const Real s844 = -(s199*s210*s225*s843);
				const Real s845 = -(s188*s223*s225*s843);
				const Real s846 = -(s189*s218);
				const Real s847 = s210*s218;
				const Real s848 = s844 + s845 + s846 + s847;
				const Real s849 = -(s199*s223*s225*s843);
				const Real s850 = -(s188*s215*s225*s843);
				const Real s851 = s609 + s610 + s849 + s850;
				const Real s852 = -(s201*s210*s225*s843);
				const Real s853 = -(s191*s223*s225*s843);
				const Real s854 = -(s188*s191*s218);
				const Real s855 = s852 + s853 + s854;
				const Real s856 = -(s201*s223*s225*s843);
				const Real s857 = -(s191*s215*s225*s843);
				const Real s858 = s671 + s672 + s856 + s857;
				const Real s859 = -(s203*s210*s225*s843);
				const Real s860 = -(s194*s223*s225*s843);
				const Real s861 = -(s188*s194*s218);
				const Real s862 = s859 + s860 + s861;
				const Real s863 = -(s203*s223*s225*s843);
				const Real s864 = -(s194*s215*s225*s843);
				const Real s865 = s730 + s731 + s863 + s864;
				const Real s866 = 2*s2*s3*s8;
				const Real s867 = -2*s19*s8;
				const Real s868 = -2*s0*s26;
				const Real s869 = -2*s0*s21;
				const Real s870 = -2*s0*s23;
				const Real s871 = 2*s0*s2*s3;
				const Real s872 = 2*s10*s4*s8;
				const Real s873 = 2*s0*s10*s4;
				const Real s874 = -2*s5*s8;
				const Real s875 = 2*s11*s6*s8;
				const Real s876 = 2*s0*s11*s6;
				const Real s877 = -2*s7*s8;
				const Real s878 = -2*s3*s8*s9;
				const Real s879 = -2*s0*s2*s9;
				const Real s880 = 2*s12*s26;
				const Real s881 = 2*s12*s21;
				const Real s882 = 2*s12*s23;
				const Real s883 = -4*s12*s2*s3;
				const Real s884 = 2*s12*s19;
				const Real s885 = -4*s10*s12*s4;
				const Real s886 = 2*s12*s5;
				const Real s887 = -4*s11*s12*s6;
				const Real s888 = 2*s12*s7;
				const Real s889 = -2*s10*s14*s8;
				const Real s890 = -2*s0*s14*s4;
				const Real s891 = -2*s11*s16*s8;
				const Real s892 = -2*s0*s16*s6;
				const Real s893 = s295 + s307 + s314 + s642 + s653 + s660 + s866 + s867 + s868 + s869 + s870 + s871 + s872 + s873 + s874 + s875 + s876 + s877 + s878 + s879 + s880 + s881 + s882 + s883 + s884 + s885 + s886 + s887 + s888 + s889 + s890 + s891 + s892;
				const Real s894 = 2*s201*s210;
				const Real s895 = -2*s191*s207;
				const Real s896 = s894 + s895;
				const Real s897 = -(s199*s210*s225*s896);
				const Real s898 = -(s188*s223*s225*s896);
				const Real s899 = s854 + s897 + s898;
				const Real s900 = -(s199*s223*s225*s896);
				const Real s901 = -(s188*s215*s225*s896);
				const Real s902 = s619 + s620 + s900 + s901;
				const Real s903 = -(s201*s210*s225*s896);
				const Real s904 = -(s191*s223*s225*s896);
				const Real s905 = -(s192*s218);
				const Real s906 = s847 + s903 + s904 + s905;
				const Real s907 = -(s201*s223*s225*s896);
				const Real s908 = -(s191*s215*s225*s896);
				const Real s909 = s610 + s679 + s907 + s908;
				const Real s910 = -(s203*s210*s225*s896);
				const Real s911 = -(s194*s223*s225*s896);
				const Real s912 = -(s191*s194*s218);
				const Real s913 = s910 + s911 + s912;
				const Real s914 = -(s203*s223*s225*s896);
				const Real s915 = -(s194*s215*s225*s896);
				const Real s916 = s738 + s739 + s914 + s915;
				const Real s917 = 2*s10*s2*s3;
				const Real s918 = -2*s10*s19;
				const Real s919 = 2*s0*s10*s8;
				const Real s920 = -2*s1*s10;
				const Real s921 = -2*s26*s4;
				const Real s922 = -2*s18*s4;
				const Real s923 = -2*s23*s4;
				const Real s924 = 2*s2*s3*s4;
				const Real s925 = 2*s0*s4*s8;
				const Real s926 = 2*s10*s11*s6;
				const Real s927 = 2*s11*s4*s6;
				const Real s928 = -2*s10*s7;
				const Real s929 = -2*s10*s3*s9;
				const Real s930 = -2*s2*s4*s9;
				const Real s931 = -2*s10*s12*s8;
				const Real s932 = -2*s0*s12*s4;
				const Real s933 = 2*s14*s26;
				const Real s934 = 2*s14*s18;
				const Real s935 = 2*s14*s23;
				const Real s936 = -4*s14*s2*s3;
				const Real s937 = 2*s14*s19;
				const Real s938 = -4*s0*s14*s8;
				const Real s939 = 2*s1*s14;
				const Real s940 = -4*s11*s14*s6;
				const Real s941 = 2*s14*s7;
				const Real s942 = -2*s10*s11*s16;
				const Real s943 = -2*s16*s4*s6;
				const Real s944 = s428 + s433 + s447 + s701 + s706 + s719 + s917 + s918 + s919 + s920 + s921 + s922 + s923 + s924 + s925 + s926 + s927 + s928 + s929 + s930 + s931 + s932 + s933 + s934 + s935 + s936 + s937 + s938 + s939 + s940 + s941 + s942 + s943;
				const Real s945 = 2*s203*s210;
				const Real s946 = -2*s194*s207;
				const Real s947 = s945 + s946;
				const Real s948 = -(s199*s210*s225*s947);
				const Real s949 = -(s188*s223*s225*s947);
				const Real s950 = s861 + s948 + s949;
				const Real s951 = -(s199*s223*s225*s947);
				const Real s952 = -(s188*s215*s225*s947);
				const Real s953 = s628 + s629 + s951 + s952;
				const Real s954 = -(s201*s210*s225*s947);
				const Real s955 = -(s191*s223*s225*s947);
				const Real s956 = s912 + s954 + s955;
				const Real s957 = -(s201*s223*s225*s947);
				const Real s958 = -(s191*s215*s225*s947);
				const Real s959 = s687 + s688 + s957 + s958;
				const Real s960 = -(s203*s210*s225*s947);
				const Real s961 = -(s194*s223*s225*s947);
				const Real s962 = -(s195*s218);
				const Real s963 = s847 + s960 + s961 + s962;
				const Real s964 = -(s203*s223*s225*s947);
				const Real s965 = -(s194*s215*s225*s947);
				const Real s966 = s610 + s746 + s964 + s965;
				const Real s967 = 2*s11*s2*s3;
				const Real s968 = -2*s11*s19;
				const Real s969 = 2*s0*s11*s8;
				const Real s970 = -2*s1*s11;
				const Real s971 = 2*s10*s11*s4;
				const Real s972 = -2*s11*s5;
				const Real s973 = -2*s26*s6;
				const Real s974 = -2*s18*s6;
				const Real s975 = -2*s21*s6;
				const Real s976 = 2*s2*s3*s6;
				const Real s977 = 2*s0*s6*s8;
				const Real s978 = 2*s10*s4*s6;
				const Real s979 = -2*s11*s3*s9;
				const Real s980 = -2*s2*s6*s9;
				const Real s981 = -2*s11*s12*s8;
				const Real s982 = -2*s0*s12*s6;
				const Real s983 = -2*s10*s11*s14;
				const Real s984 = -2*s14*s4*s6;
				const Real s985 = 2*s16*s26;
				const Real s986 = 2*s16*s18;
				const Real s987 = 2*s16*s21;
				const Real s988 = -4*s16*s2*s3;
				const Real s989 = 2*s16*s19;
				const Real s990 = -4*s0*s16*s8;
				const Real s991 = 2*s1*s16;
				const Real s992 = -4*s10*s16*s4;
				const Real s993 = 2*s16*s5;
				const Real s994 = s501 + s506 + s511 + s759 + s764 + s769 + s967 + s968 + s969 + s970 + s971 + s972 + s973 + s974 + s975 + s976 + s977 + s978 + s979 + s980 + s981 + s982 + s983 + s984 + s985 + s986 + s987 + s988 + s989 + s990 + s991 + s992 + s993;
				const Real s995 = 2*s205*s210;
				const Real s996 = -2*s197*s207;
				const Real s997 = s995 + s996;
				const Real s998 = -(s188*s197*s218);
				const Real s999 = -(s197*s199*s218);
				const Real s1000 = 2*s188*s205*s218;
				const Real s1001 = -(s199*s210*s225*s997);
				const Real s1002 = -(s188*s223*s225*s997);
				const Real s1003 = s1001 + s1002 + s998;
				const Real s1004 = -(s199*s223*s225*s997);
				const Real s1005 = -(s188*s215*s225*s997);
				const Real s1006 = s1000 + s1004 + s1005 + s999;
				const Real s1007 = -(s191*s197*s218);
				const Real s1008 = -(s197*s201*s218);
				const Real s1009 = 2*s191*s205*s218;
				const Real s1010 = -(s201*s210*s225*s997);
				const Real s1011 = -(s191*s223*s225*s997);
				const Real s1012 = s1007 + s1010 + s1011;
				const Real s1013 = -(s201*s223*s225*s997);
				const Real s1014 = -(s191*s215*s225*s997);
				const Real s1015 = s1008 + s1009 + s1013 + s1014;
				const Real s1016 = -(s194*s197*s218);
				const Real s1017 = -(s197*s203*s218);
				const Real s1018 = 2*s194*s205*s218;
				const Real s1019 = -(s203*s210*s225*s997);
				const Real s1020 = -(s194*s223*s225*s997);
				const Real s1021 = s1016 + s1019 + s1020;
				const Real s1022 = -(s203*s223*s225*s997);
				const Real s1023 = -(s194*s215*s225*s997);
				const Real s1024 = s1017 + s1018 + s1022 + s1023;
				const Real s1025 = s197*s205*s218;
				buffer__[12*i+0] = (s152*s186*s286)/4. + (s10*s152*s186*s322)/4. + (s11*s152*s186*s323)/4. + (s152*s186*s2*s324)/4. + (s0*s152*s186*s325)/4. + (s12*s152*s186*s329)/4. + (s14*s152*s186*s330)/4. + (s152*s16*s186*s331)/4. + s287*((s152*s186*s3)/4. + s333) + s334*((s152*s186*s347)/4. + (s151*(s199*s242 + s188*s247 + s349 + s350 + s351 + s352))/2.) + s348*((s151*(s201*s242 + s191*s247))/2. + (s152*s186*s363)/4.) + s364*((s151*(s203*s242 + s194*s247))/2. + (s152*s186*s367)/4.) + s368*((s151*(s205*s242 + s197*s247))/2. + (s152*s186*s371)/4.) + s372*((s151*(s201*s258 + s191*s263))/2. + (s152*s186*s376)/4.) + s377*((s151*(s203*s258 + s194*s263))/2. + (s152*s186*s388)/4.) + s389*((s151*(s205*s258 + s197*s263))/2. + (s152*s186*s392)/4.) + s393*((s151*(s203*s274 + s194*s279))/2. + (s152*s186*s396)/4.) + (s152*s186*s326*s4)/4. + s397*((s151*(s205*s274 + s197*s279))/2. + (s152*s186*s410)/4.) + s411*((s152*s186*s420)/4. + (s151*(s197*(-(s197*s215*s225*s233) - s205*s223*s225*s233 + s205*s218*s236 + s536) + s205*(-(s205*s210*s225*s233) - s197*s223*s225*s233 + s197*s218*s236 + s537)))/2.) + (s152*s186*s327*s6)/4. + (s152*s186*s332*s8)/4. + (s152*s186*s328*s9)/4.;
				buffer__[12*i+1] = (s152*s286*s321)/4. + (s152*s287*s3*s321)/4. + (s10*s152*s321*s322)/4. + (s11*s152*s321*s323)/4. + (s152*s2*s321*s324)/4. + (s0*s152*s321*s325)/4. + (s12*s152*s321*s329)/4. + (s14*s152*s321*s330)/4. + (s152*s16*s321*s331)/4. + s334*((s152*s321*s347)/4. + (s151*(s199*s356 + s188*s360))/2.) + s348*((s151*(s349 + s350 + s351 + s352 + s201*s356 + s191*s360))/2. + (s152*s321*s363)/4.) + s364*((s151*(s203*s356 + s194*s360))/2. + (s152*s321*s367)/4.) + s368*((s151*(s205*s356 + s197*s360))/2. + (s152*s321*s371)/4.) + s377*((s151*(s203*s381 + s194*s385))/2. + (s152*s321*s388)/4.) + s389*((s151*(s205*s381 + s197*s385))/2. + (s152*s321*s392)/4.) + (s152*s321*s326*s4)/4. + s393*((s152*s321*s396)/4. + (s151*(s203*s402 + s194*s407))/2.) + s397*((s151*(s205*s402 + s197*s407))/2. + (s152*s321*s410)/4.) + s372*((s152*s321*s376)/4. + (s151*(s201*s381 + s191*s385 + s473 + s474 + s475 + s476))/2.) + s411*((s152*s321*s420)/4. + (s151*(s197*(-(s197*s215*s225*s342) - s205*s223*s225*s342 + s205*s218*s344 + s548) + s205*(-(s205*s210*s225*s342) - s197*s223*s225*s342 + s197*s218*s344 + s549)))/2.) + (s152*s321*s327*s6)/4. + s332*(s333 + (s152*s321*s8)/4.) + (s152*s321*s328*s9)/4.;
				buffer__[12*i+2] = (s152*s286*s454)/4. + (s152*s287*s3*s454)/4. + (s11*s152*s323*s454)/4. + (s152*s2*s324*s454)/4. + (s0*s152*s325*s454)/4. + (s12*s152*s329*s454)/4. + (s14*s152*s330*s454)/4. + (s152*s16*s331*s454)/4. + (s152*s326*s4*s454)/4. + s322*(s333 + (s10*s152*s454)/4.) + s334*((s152*s347*s454)/4. + (s151*(s199*s468 + s188*s472))/2.) + s348*((s152*s363*s454)/4. + (s151*(s201*s468 + s191*s472))/2.) + s364*((s152*s367*s454)/4. + (s151*(s349 + s350 + s351 + s352 + s203*s468 + s194*s472))/2.) + s368*((s152*s371*s454)/4. + (s151*(s205*s468 + s197*s472))/2.) + s372*((s152*s376*s454)/4. + (s151*(s201*s480 + s191*s484))/2.) + s377*((s152*s388*s454)/4. + (s151*(s473 + s474 + s475 + s476 + s203*s480 + s194*s484))/2.) + s389*((s152*s392*s454)/4. + (s151*(s205*s480 + s197*s484))/2.) + s397*((s152*s410*s454)/4. + (s151*(s205*s489 + s197*s493))/2.) + s411*((s152*s420*s454)/4. + (s151*(s197*(-(s197*s215*s225*s462) - s205*s223*s225*s462 + s205*s218*s464 + s558) + s205*(-(s205*s210*s225*s462) - s197*s223*s225*s462 + s197*s218*s464 + s559)))/2.) + s393*((s152*s396*s454)/4. + (s151*(s203*s489 + s194*s493 + s560 + s561 + s562 + s563))/2.) + (s152*s327*s454*s6)/4. + (s152*s332*s454*s8)/4. + (s152*s328*s454*s9)/4.;
				buffer__[12*i+3] = (s152*s286*s527)/4. + (s152*s287*s3*s527)/4. + (s10*s152*s322*s527)/4. + (s152*s2*s324*s527)/4. + (s0*s152*s325*s527)/4. + (s12*s152*s329*s527)/4. + (s14*s152*s330*s527)/4. + (s152*s16*s331*s527)/4. + (s152*s326*s4*s527)/4. + s323*(s333 + (s11*s152*s527)/4.) + s334*((s152*s347*s527)/4. + (s151*(s199*s543 + s188*s547))/2.) + s348*((s152*s363*s527)/4. + (s151*(s201*s543 + s191*s547))/2.) + s364*((s152*s367*s527)/4. + (s151*(s203*s543 + s194*s547))/2.) + s368*((s152*s371*s527)/4. + (s151*(s349 + s350 + s351 + s352 + s205*s543 + s197*s547))/2.) + s372*((s152*s376*s527)/4. + (s151*(s201*s553 + s191*s557))/2.) + s377*((s152*s388*s527)/4. + (s151*(s203*s553 + s194*s557))/2.) + s389*((s152*s392*s527)/4. + (s151*(s473 + s474 + s475 + s476 + s205*s553 + s197*s557))/2.) + s393*((s152*s396*s527)/4. + (s151*(s203*s567 + s194*s571))/2.) + s397*((s152*s410*s527)/4. + (s151*(s560 + s561 + s562 + s563 + s205*s567 + s197*s571))/2.) + s411*((s152*s420*s527)/4. + (s151*(-(s205*s210*s218) - s197*s215*s218 - s197*s218*s223 - s205*s218*s223 + s205*(s237 + s240 - s205*s210*s225*s535 - s197*s223*s225*s535 + s197*s218*s539 + s572) + s197*(s237 + s246 - s197*s215*s225*s535 - s205*s223*s225*s535 + s205*s218*s539 + s572)))/2.) + (s152*s327*s527*s6)/4. + (s152*s332*s527*s8)/4. + (s152*s328*s527*s9)/4.;
				buffer__[12*i+4] = (s152*s286*s603)/4. + (s152*s287*s3*s603)/4. + (s10*s152*s322*s603)/4. + (s11*s152*s323*s603)/4. + (s0*s152*s325*s603)/4. + (s12*s152*s329*s603)/4. + (s14*s152*s330*s603)/4. + (s152*s16*s331*s603)/4. + (s152*s326*s4*s603)/4. + (s152*s327*s6*s603)/4. + s324*(s333 + (s152*s2*s603)/4.) + s334*((s152*s347*s603)/4. + (s151*(s251 + s252 + s199*s611 + s188*s616))/2.) + s348*((s152*s363*s603)/4. + (s151*(s201*s611 + s191*s616))/2.) + s364*((s152*s367*s603)/4. + (s151*(s203*s611 + s194*s616))/2.) + s368*((s152*s371*s603)/4. + (s151*(s205*s611 + s197*s616))/2.) + s372*((s152*s376*s603)/4. + (s151*(s201*s621 + s191*s625))/2.) + s377*((s152*s388*s603)/4. + (s151*(s203*s621 + s194*s625))/2.) + s389*((s152*s392*s603)/4. + (s151*(s205*s621 + s197*s625))/2.) + s393*((s152*s396*s603)/4. + (s151*(s203*s630 + s194*s634))/2.) + s397*((s152*s410*s603)/4. + (s151*(s205*s630 + s197*s634))/2.) + (s152*s332*s603*s8)/4. + (s152*s328*s603*s9)/4. + s411*((s152*s420*s603)/4. + (s151*(s197*(-(s197*s215*s225*s606) - s205*s223*s225*s606 + s786) + s205*(s1000 - s205*s210*s225*s606 - s197*s223*s225*s606 + s999)))/2.);
				buffer__[12*i+5] = (s152*s286*s665)/4. + (s152*s287*s3*s665)/4. + (s10*s152*s322*s665)/4. + (s11*s152*s323*s665)/4. + (s152*s2*s324*s665)/4. + (s12*s152*s329*s665)/4. + (s14*s152*s330*s665)/4. + (s152*s16*s331*s665)/4. + (s152*s326*s4*s665)/4. + (s152*s327*s6*s665)/4. + s325*(s333 + (s0*s152*s665)/4.) + s334*((s152*s347*s665)/4. + (s151*(s199*s673 + s188*s676))/2.) + s348*((s152*s363*s665)/4. + (s151*(s251 + s252 + s201*s673 + s191*s676))/2.) + s364*((s152*s367*s665)/4. + (s151*(s203*s673 + s194*s676))/2.) + s368*((s152*s371*s665)/4. + (s151*(s205*s673 + s197*s676))/2.) + s372*((s152*s376*s665)/4. + (s151*(s267 + s268 + s201*s680 + s191*s684))/2.) + s377*((s152*s388*s665)/4. + (s151*(s203*s680 + s194*s684))/2.) + s389*((s152*s392*s665)/4. + (s151*(s205*s680 + s197*s684))/2.) + s393*((s152*s396*s665)/4. + (s151*(s203*s689 + s194*s693))/2.) + s397*((s152*s410*s665)/4. + (s151*(s205*s689 + s197*s693))/2.) + s411*((s152*s420*s665)/4. + (s151*(s205*(s1008 + s1009 - s205*s210*s225*s668 - s197*s223*s225*s668) + s197*(-(s197*s215*s225*s668) - s205*s223*s225*s668 + s795)))/2.) + (s152*s332*s665*s8)/4. + (s152*s328*s665*s9)/4.;
				buffer__[12*i+6] = (s152*s286*s724)/4. + (s152*s287*s3*s724)/4. + (s10*s152*s322*s724)/4. + (s11*s152*s323*s724)/4. + (s152*s2*s324*s724)/4. + (s0*s152*s325*s724)/4. + (s12*s152*s329*s724)/4. + (s14*s152*s330*s724)/4. + (s152*s16*s331*s724)/4. + (s152*s327*s6*s724)/4. + s326*(s333 + (s152*s4*s724)/4.) + s334*((s152*s347*s724)/4. + (s151*(s199*s732 + s188*s735))/2.) + s348*((s152*s363*s724)/4. + (s151*(s201*s732 + s191*s735))/2.) + s364*((s152*s367*s724)/4. + (s151*(s251 + s252 + s203*s732 + s194*s735))/2.) + s368*((s152*s371*s724)/4. + (s151*(s205*s732 + s197*s735))/2.) + s372*((s152*s376*s724)/4. + (s151*(s201*s740 + s191*s743))/2.) + s377*((s152*s388*s724)/4. + (s151*(s267 + s268 + s203*s740 + s194*s743))/2.) + s389*((s152*s392*s724)/4. + (s151*(s205*s740 + s197*s743))/2.) + s393*((s152*s396*s724)/4. + (s151*(s283 + s284 + s203*s747 + s194*s751))/2.) + s397*((s152*s410*s724)/4. + (s151*(s205*s747 + s197*s751))/2.) + (s152*s332*s724*s8)/4. + s411*((s152*s420*s724)/4. + (s151*(s205*(s1017 + s1018 - s205*s210*s225*s727 - s197*s223*s225*s727) + s197*(-(s197*s215*s225*s727) - s205*s223*s225*s727 + s804)))/2.) + (s152*s328*s724*s9)/4.;
				buffer__[12*i+7] = (s152*s286*s782)/4. + (s152*s287*s3*s782)/4. + (s10*s152*s322*s782)/4. + (s11*s152*s323*s782)/4. + (s152*s2*s324*s782)/4. + (s0*s152*s325*s782)/4. + (s12*s152*s329*s782)/4. + (s14*s152*s330*s782)/4. + (s152*s16*s331*s782)/4. + (s152*s326*s4*s782)/4. + s327*(s333 + (s152*s6*s782)/4.) + s411*((s152*s420*s782)/4. + (s151*(s416 + s417 + s205*(s1025 + s610 - s205*s210*s225*s785 - s197*s223*s225*s785) + s197*(-(s214*s218) + s615 - s197*s215*s225*s785 - s205*s223*s225*s785)))/2.) + s334*((s152*s347*s782)/4. + (s151*(s199*s791 + s188*s794))/2.) + s348*((s152*s363*s782)/4. + (s151*(s201*s791 + s191*s794))/2.) + s364*((s152*s367*s782)/4. + (s151*(s203*s791 + s194*s794))/2.) + s368*((s152*s371*s782)/4. + (s151*(s251 + s252 + s205*s791 + s197*s794))/2.) + (s152*s332*s782*s8)/4. + s372*((s152*s376*s782)/4. + (s151*(s201*s800 + s191*s803))/2.) + s377*((s152*s388*s782)/4. + (s151*(s203*s800 + s194*s803))/2.) + s389*((s152*s392*s782)/4. + (s151*(s267 + s268 + s205*s800 + s197*s803))/2.) + s393*((s152*s396*s782)/4. + (s151*(s203*s809 + s194*s812))/2.) + s397*((s152*s410*s782)/4. + (s151*(s283 + s284 + s205*s809 + s197*s812))/2.) + (s152*s328*s782*s9)/4.;
				buffer__[12*i+8] = (s152*s286*s840)/4. + (s152*s287*s3*s840)/4. + (s10*s152*s322*s840)/4. + (s11*s152*s323*s840)/4. + (s152*s2*s324*s840)/4. + (s0*s152*s325*s840)/4. + (s12*s152*s329*s840)/4. + (s14*s152*s330*s840)/4. + (s152*s16*s331*s840)/4. + (s152*s326*s4*s840)/4. + (s152*s327*s6*s840)/4. + (s152*s332*s8*s840)/4. + s334*((s152*s347*s840)/4. + (s151*(s248 + s249 + s199*s848 + s188*s851))/2.) + s348*((s152*s363*s840)/4. + (s151*(s201*s848 + s191*s851))/2.) + s364*((s152*s367*s840)/4. + (s151*(s203*s848 + s194*s851))/2.) + s368*((s152*s371*s840)/4. + (s151*(s205*s848 + s197*s851))/2.) + s372*((s152*s376*s840)/4. + (s151*(s201*s855 + s191*s858))/2.) + s377*((s152*s388*s840)/4. + (s151*(s203*s855 + s194*s858))/2.) + s389*((s152*s392*s840)/4. + (s151*(s205*s855 + s197*s858))/2.) + s393*((s152*s396*s840)/4. + (s151*(s203*s862 + s194*s865))/2.) + s397*((s152*s410*s840)/4. + (s151*(s205*s862 + s197*s865))/2.) + s328*(s333 + (s152*s840*s9)/4.) + s411*((s152*s420*s840)/4. + (s151*(s197*(s789 + s790 - s197*s215*s225*s843 - s205*s223*s225*s843) + s205*(-(s205*s210*s225*s843) - s197*s223*s225*s843 + s998)))/2.);
				buffer__[12*i+9] = (s152*s286*s893)/4. + (s152*s287*s3*s893)/4. + (s10*s152*s322*s893)/4. + (s11*s152*s323*s893)/4. + (s152*s2*s324*s893)/4. + (s0*s152*s325*s893)/4. + (s14*s152*s330*s893)/4. + (s152*s16*s331*s893)/4. + (s152*s326*s4*s893)/4. + (s152*s327*s6*s893)/4. + (s152*s332*s8*s893)/4. + s329*(s333 + (s12*s152*s893)/4.) + s411*((s152*s420*s893)/4. + (s151*(s205*(s1007 - s205*s210*s225*s896 - s197*s223*s225*s896) + s197*(s798 + s799 - s197*s215*s225*s896 - s205*s223*s225*s896)))/2.) + (s152*s328*s893*s9)/4. + s334*((s152*s347*s893)/4. + (s151*(s199*s899 + s188*s902))/2.) + s348*((s152*s363*s893)/4. + (s151*(s248 + s249 + s201*s899 + s191*s902))/2.) + s364*((s152*s367*s893)/4. + (s151*(s203*s899 + s194*s902))/2.) + s368*((s152*s371*s893)/4. + (s151*(s205*s899 + s197*s902))/2.) + s372*((s152*s376*s893)/4. + (s151*(s264 + s265 + s201*s906 + s191*s909))/2.) + s377*((s152*s388*s893)/4. + (s151*(s203*s906 + s194*s909))/2.) + s389*((s152*s392*s893)/4. + (s151*(s205*s906 + s197*s909))/2.) + s393*((s152*s396*s893)/4. + (s151*(s203*s913 + s194*s916))/2.) + s397*((s152*s410*s893)/4. + (s151*(s205*s913 + s197*s916))/2.);
				buffer__[12*i+10] = (s152*s286*s944)/4. + (s152*s287*s3*s944)/4. + (s10*s152*s322*s944)/4. + (s11*s152*s323*s944)/4. + (s152*s2*s324*s944)/4. + (s0*s152*s325*s944)/4. + (s12*s152*s329*s944)/4. + (s152*s16*s331*s944)/4. + (s152*s326*s4*s944)/4. + (s152*s327*s6*s944)/4. + (s152*s332*s8*s944)/4. + (s152*s328*s9*s944)/4. + s330*(s333 + (s14*s152*s944)/4.) + s411*((s152*s420*s944)/4. + (s151*(s205*(s1016 - s205*s210*s225*s947 - s197*s223*s225*s947) + s197*(s807 + s808 - s197*s215*s225*s947 - s205*s223*s225*s947)))/2.) + s334*((s152*s347*s944)/4. + (s151*(s199*s950 + s188*s953))/2.) + s348*((s152*s363*s944)/4. + (s151*(s201*s950 + s191*s953))/2.) + s364*((s152*s367*s944)/4. + (s151*(s248 + s249 + s203*s950 + s194*s953))/2.) + s368*((s152*s371*s944)/4. + (s151*(s205*s950 + s197*s953))/2.) + s372*((s152*s376*s944)/4. + (s151*(s201*s956 + s191*s959))/2.) + s377*((s152*s388*s944)/4. + (s151*(s264 + s265 + s203*s956 + s194*s959))/2.) + s389*((s152*s392*s944)/4. + (s151*(s205*s956 + s197*s959))/2.) + s393*((s152*s396*s944)/4. + (s151*(s280 + s281 + s203*s963 + s194*s966))/2.) + s397*((s152*s410*s944)/4. + (s151*(s205*s963 + s197*s966))/2.);
				buffer__[12*i+11] = (s152*s286*s994)/4. + (s152*s287*s3*s994)/4. + (s10*s152*s322*s994)/4. + (s11*s152*s323*s994)/4. + (s152*s2*s324*s994)/4. + (s0*s152*s325*s994)/4. + (s12*s152*s329*s994)/4. + (s14*s152*s330*s994)/4. + (s152*s326*s4*s994)/4. + (s152*s327*s6*s994)/4. + (s152*s332*s8*s994)/4. + (s152*s328*s9*s994)/4. + s331*(s333 + (s152*s16*s994)/4.) + s334*((s151*(s1006*s188 + s1003*s199))/2. + (s152*s347*s994)/4.) + s348*((s151*(s1006*s191 + s1003*s201))/2. + (s152*s363*s994)/4.) + s364*((s151*(s1006*s194 + s1003*s203))/2. + (s152*s367*s994)/4.) + s368*((s151*(s1006*s197 + s1003*s205 + s248 + s249))/2. + (s152*s371*s994)/4.) + s372*((s151*(s1015*s191 + s1012*s201))/2. + (s152*s376*s994)/4.) + s377*((s151*(s1015*s194 + s1012*s203))/2. + (s152*s388*s994)/4.) + s389*((s151*(s1015*s197 + s1012*s205 + s264 + s265))/2. + (s152*s392*s994)/4.) + s393*((s151*(s1024*s194 + s1021*s203))/2. + (s152*s396*s994)/4.) + s397*((s151*(s1024*s197 + s1021*s205 + s280 + s281))/2. + (s152*s410*s994)/4.) + s411*((s152*s420*s994)/4. + (s151*(s412 + s413 + s205*(-(s198*s218) + s847 - s205*s210*s225*s997 - s197*s223*s225*s997) + s197*(s1025 + s610 - s197*s215*s225*s997 - s205*s223*s225*s997)))/2.);
			}
		}

        ptoc(ClassName()+"::DNearToHulls");
        
    }

    void DFarToHulls( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_far, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+"::DFarToHulls");

        if( P_D_far.Dimension(1) != 15 )
        {
            eprint("in DFarToHulls: P_D_far.Dimension(1) != 15. Aborting");
        }
        
		ptr<Real> V_coords__  = V_coords.data();
        ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_far__   = P_D_far.data();
              Real * restrict const buffer__    = buffer.data();

		if( addTo )
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[4*simplices__[3*i+1]+1];
				const Real s1 = s0*s0;
				const Real s2 = V_coords__[4*simplices__[3*i+1]+0];
				const Real s3 = V_coords__[4*simplices__[3*i+0]+0];
				const Real s4 = V_coords__[4*simplices__[3*i+1]+2];
				const Real s5 = s4*s4;
				const Real s6 = V_coords__[4*simplices__[3*i+1]+3];
				const Real s7 = s6*s6;
				const Real s8 = V_coords__[4*simplices__[3*i+0]+1];
				const Real s9 = V_coords__[4*simplices__[3*i+2]+0];
				const Real s10 = V_coords__[4*simplices__[3*i+0]+2];
				const Real s11 = V_coords__[4*simplices__[3*i+0]+3];
				const Real s12 = V_coords__[4*simplices__[3*i+2]+1];
				const Real s13 = s12*s12;
				const Real s14 = V_coords__[4*simplices__[3*i+2]+2];
				const Real s15 = s14*s14;
				const Real s16 = V_coords__[4*simplices__[3*i+2]+3];
				const Real s17 = s16*s16;
				const Real s18 = s8*s8;
				const Real s19 = s2*s2;
				const Real s20 = s18*s19;
				const Real s21 = s10*s10;
				const Real s22 = s19*s21;
				const Real s23 = s11*s11;
				const Real s24 = s19*s23;
				const Real s25 = -2*s0*s2*s3*s8;
				const Real s26 = s3*s3;
				const Real s27 = s1*s26;
				const Real s28 = s1*s21;
				const Real s29 = s1*s23;
				const Real s30 = -2*s10*s2*s3*s4;
				const Real s31 = -2*s0*s10*s4*s8;
				const Real s32 = s26*s5;
				const Real s33 = s18*s5;
				const Real s34 = s23*s5;
				const Real s35 = -2*s11*s2*s3*s6;
				const Real s36 = -2*s0*s11*s6*s8;
				const Real s37 = -2*s10*s11*s4*s6;
				const Real s38 = s26*s7;
				const Real s39 = s18*s7;
				const Real s40 = s21*s7;
				const Real s41 = -2*s18*s2*s9;
				const Real s42 = -2*s2*s21*s9;
				const Real s43 = -2*s2*s23*s9;
				const Real s44 = 2*s0*s3*s8*s9;
				const Real s45 = 2*s0*s2*s8*s9;
				const Real s46 = -2*s1*s3*s9;
				const Real s47 = 2*s10*s3*s4*s9;
				const Real s48 = 2*s10*s2*s4*s9;
				const Real s49 = -2*s3*s5*s9;
				const Real s50 = 2*s11*s3*s6*s9;
				const Real s51 = 2*s11*s2*s6*s9;
				const Real s52 = -2*s3*s7*s9;
				const Real s53 = s9*s9;
				const Real s54 = s18*s53;
				const Real s55 = s21*s53;
				const Real s56 = s23*s53;
				const Real s57 = -2*s0*s53*s8;
				const Real s58 = s1*s53;
				const Real s59 = -2*s10*s4*s53;
				const Real s60 = s5*s53;
				const Real s61 = -2*s11*s53*s6;
				const Real s62 = s53*s7;
				const Real s63 = 2*s12*s2*s3*s8;
				const Real s64 = -2*s12*s19*s8;
				const Real s65 = -2*s0*s12*s26;
				const Real s66 = -2*s0*s12*s21;
				const Real s67 = -2*s0*s12*s23;
				const Real s68 = 2*s0*s12*s2*s3;
				const Real s69 = 2*s10*s12*s4*s8;
				const Real s70 = 2*s0*s10*s12*s4;
				const Real s71 = -2*s12*s5*s8;
				const Real s72 = 2*s11*s12*s6*s8;
				const Real s73 = 2*s0*s11*s12*s6;
				const Real s74 = -2*s12*s7*s8;
				const Real s75 = -2*s12*s3*s8*s9;
				const Real s76 = 2*s12*s2*s8*s9;
				const Real s77 = 2*s0*s12*s3*s9;
				const Real s78 = -2*s0*s12*s2*s9;
				const Real s79 = s13*s26;
				const Real s80 = s13*s21;
				const Real s81 = s13*s23;
				const Real s82 = -2*s13*s2*s3;
				const Real s83 = s13*s19;
				const Real s84 = -2*s10*s13*s4;
				const Real s85 = s13*s5;
				const Real s86 = -2*s11*s13*s6;
				const Real s87 = s13*s7;
				const Real s88 = 2*s10*s14*s2*s3;
				const Real s89 = -2*s10*s14*s19;
				const Real s90 = 2*s0*s10*s14*s8;
				const Real s91 = -2*s1*s10*s14;
				const Real s92 = -2*s14*s26*s4;
				const Real s93 = -2*s14*s18*s4;
				const Real s94 = -2*s14*s23*s4;
				const Real s95 = 2*s14*s2*s3*s4;
				const Real s96 = 2*s0*s14*s4*s8;
				const Real s97 = 2*s10*s11*s14*s6;
				const Real s98 = 2*s11*s14*s4*s6;
				const Real s99 = -2*s10*s14*s7;
				const Real s100 = -2*s10*s14*s3*s9;
				const Real s101 = 2*s10*s14*s2*s9;
				const Real s102 = 2*s14*s3*s4*s9;
				const Real s103 = -2*s14*s2*s4*s9;
				const Real s104 = -2*s10*s12*s14*s8;
				const Real s105 = 2*s0*s10*s12*s14;
				const Real s106 = 2*s12*s14*s4*s8;
				const Real s107 = -2*s0*s12*s14*s4;
				const Real s108 = s15*s26;
				const Real s109 = s15*s18;
				const Real s110 = s15*s23;
				const Real s111 = -2*s15*s2*s3;
				const Real s112 = s15*s19;
				const Real s113 = -2*s0*s15*s8;
				const Real s114 = s1*s15;
				const Real s115 = -2*s11*s15*s6;
				const Real s116 = s15*s7;
				const Real s117 = 2*s11*s16*s2*s3;
				const Real s118 = -2*s11*s16*s19;
				const Real s119 = 2*s0*s11*s16*s8;
				const Real s120 = -2*s1*s11*s16;
				const Real s121 = 2*s10*s11*s16*s4;
				const Real s122 = -2*s11*s16*s5;
				const Real s123 = -2*s16*s26*s6;
				const Real s124 = -2*s16*s18*s6;
				const Real s125 = -2*s16*s21*s6;
				const Real s126 = 2*s16*s2*s3*s6;
				const Real s127 = 2*s0*s16*s6*s8;
				const Real s128 = 2*s10*s16*s4*s6;
				const Real s129 = -2*s11*s16*s3*s9;
				const Real s130 = 2*s11*s16*s2*s9;
				const Real s131 = 2*s16*s3*s6*s9;
				const Real s132 = -2*s16*s2*s6*s9;
				const Real s133 = -2*s11*s12*s16*s8;
				const Real s134 = 2*s0*s11*s12*s16;
				const Real s135 = 2*s12*s16*s6*s8;
				const Real s136 = -2*s0*s12*s16*s6;
				const Real s137 = -2*s10*s11*s14*s16;
				const Real s138 = 2*s11*s14*s16*s4;
				const Real s139 = 2*s10*s14*s16*s6;
				const Real s140 = -2*s14*s16*s4*s6;
				const Real s141 = s17*s26;
				const Real s142 = s17*s18;
				const Real s143 = s17*s21;
				const Real s144 = -2*s17*s2*s3;
				const Real s145 = s17*s19;
				const Real s146 = -2*s0*s17*s8;
				const Real s147 = s1*s17;
				const Real s148 = -2*s10*s17*s4;
				const Real s149 = s17*s5;
				const Real s150 = s100 + s101 + s102 + s103 + s104 + s105 + s106 + s107 + s108 + s109 + s110 + s111 + s112 + s113 + s114 + s115 + s116 + s117 + s118 + s119 + s120 + s121 + s122 + s123 + s124 + s125 + s126 + s127 + s128 + s129 + s130 + s131 + s132 + s133 + s134 + s135 + s136 + s137 + s138 + s139 + s140 + s141 + s142 + s143 + s144 + s145 + s146 + s147 + s148 + s149 + s20 + s22 + s24 + s25 + s27 + s28 + s29 + s30 + s31 + s32 + s33 + s34 + s35 + s36 + s37 + s38 + s39 + s40 + s41 + s42 + s43 + s44 + s45 + s46 + s47 + s48 + s49 + s50 + s51 + s52 + s54 + s55 + s56 + s57 + s58 + s59 + s60 + s61 + s62 + s63 + s64 + s65 + s66 + s67 + s68 + s69 + s70 + s71 + s72 + s73 + s74 + s75 + s76 + s77 + s78 + s79 + s80 + s81 + s82 + s83 + s84 + s85 + s86 + s87 + s88 + s89 + s90 + s91 + s92 + s93 + s94 + s95 + s96 + s97 + s98 + s99;
				const Real s151 = sqrt(s150);
				const Real s152 = 1/s151;
				const Real s153 = -2*s0*s2*s8;
				const Real s154 = 2*s1*s3;
				const Real s155 = -2*s10*s2*s4;
				const Real s156 = 2*s3*s5;
				const Real s157 = -2*s11*s2*s6;
				const Real s158 = 2*s3*s7;
				const Real s159 = 2*s0*s8*s9;
				const Real s160 = -2*s1*s9;
				const Real s161 = 2*s10*s4*s9;
				const Real s162 = -2*s5*s9;
				const Real s163 = 2*s11*s6*s9;
				const Real s164 = -2*s7*s9;
				const Real s165 = 2*s12*s2*s8;
				const Real s166 = -4*s0*s12*s3;
				const Real s167 = 2*s0*s12*s2;
				const Real s168 = -2*s12*s8*s9;
				const Real s169 = 2*s0*s12*s9;
				const Real s170 = 2*s13*s3;
				const Real s171 = -2*s13*s2;
				const Real s172 = 2*s10*s14*s2;
				const Real s173 = -4*s14*s3*s4;
				const Real s174 = 2*s14*s2*s4;
				const Real s175 = -2*s10*s14*s9;
				const Real s176 = 2*s14*s4*s9;
				const Real s177 = 2*s15*s3;
				const Real s178 = -2*s15*s2;
				const Real s179 = 2*s11*s16*s2;
				const Real s180 = -4*s16*s3*s6;
				const Real s181 = 2*s16*s2*s6;
				const Real s182 = -2*s11*s16*s9;
				const Real s183 = 2*s16*s6*s9;
				const Real s184 = 2*s17*s3;
				const Real s185 = -2*s17*s2;
				const Real s186 = s153 + s154 + s155 + s156 + s157 + s158 + s159 + s160 + s161 + s162 + s163 + s164 + s165 + s166 + s167 + s168 + s169 + s170 + s171 + s172 + s173 + s174 + s175 + s176 + s177 + s178 + s179 + s180 + s181 + s182 + s183 + s184 + s185;
				const Real s187 = -s3;
				const Real s188 = s187 + s2;
				const Real s189 = s188*s188;
				const Real s190 = -s8;
				const Real s191 = s0 + s190;
				const Real s192 = s191*s191;
				const Real s193 = -s10;
				const Real s194 = s193 + s4;
				const Real s195 = s194*s194;
				const Real s196 = -s11;
				const Real s197 = s196 + s6;
				const Real s198 = s197*s197;
				const Real s199 = s187 + s9;
				const Real s200 = s188*s199;
				const Real s201 = s12 + s190;
				const Real s202 = s191*s201;
				const Real s203 = s14 + s193;
				const Real s204 = s194*s203;
				const Real s205 = s16 + s196;
				const Real s206 = s197*s205;
				const Real s207 = s200 + s202 + s204 + s206;
				const Real s208 = s207*s207;
				const Real s209 = -s208;
				const Real s210 = s189 + s192 + s195 + s198;
				const Real s211 = s199*s199;
				const Real s212 = s201*s201;
				const Real s213 = s203*s203;
				const Real s214 = s205*s205;
				const Real s215 = s211 + s212 + s213 + s214;
				const Real s216 = s210*s215;
				const Real s217 = s209 + s216;
				const Real s218 = 1/s217;
				const Real s219 = -(s188*s199);
				const Real s220 = -(s191*s201);
				const Real s221 = -(s194*s203);
				const Real s222 = -(s197*s205);
				const Real s223 = s219 + s220 + s221 + s222;
				const Real s224 = s217*s217;
				const Real s225 = 1/s224;
				const Real s226 = -2*s199*s210;
				const Real s227 = 2*s3;
				const Real s228 = -s2;
				const Real s229 = -s9;
				const Real s230 = s227 + s228 + s229;
				const Real s231 = -2*s207*s230;
				const Real s232 = -2*s188*s215;
				const Real s233 = s226 + s231 + s232;
				const Real s234 = -2*s188*s199*s218;
				const Real s235 = -2*s3;
				const Real s236 = s2 + s235 + s9;
				const Real s237 = -(s218*s223);
				const Real s238 = -(s199*s210*s225*s233);
				const Real s239 = -(s188*s223*s225*s233);
				const Real s240 = -(s210*s218);
				const Real s241 = s188*s218*s236;
				const Real s242 = s234 + s237 + s238 + s239 + s240 + s241;
				const Real s243 = -(s199*s223*s225*s233);
				const Real s244 = -(s188*s215*s225*s233);
				const Real s245 = s199*s218*s236;
				const Real s246 = -(s215*s218);
				const Real s247 = s234 + s237 + s243 + s244 + s245 + s246;
				const Real s248 = s199*s210*s218;
				const Real s249 = s188*s218*s223;
				const Real s250 = s248 + s249;
				const Real s251 = s199*s218*s223;
				const Real s252 = s188*s215*s218;
				const Real s253 = s251 + s252;
				const Real s254 = -(s201*s210*s225*s233);
				const Real s255 = -(s191*s223*s225*s233);
				const Real s256 = s191*s218*s236;
				const Real s257 = -2*s188*s201*s218;
				const Real s258 = s254 + s255 + s256 + s257;
				const Real s259 = -(s201*s223*s225*s233);
				const Real s260 = -(s191*s215*s225*s233);
				const Real s261 = -2*s191*s199*s218;
				const Real s262 = s201*s218*s236;
				const Real s263 = s259 + s260 + s261 + s262;
				const Real s264 = s201*s210*s218;
				const Real s265 = s191*s218*s223;
				const Real s266 = s264 + s265;
				const Real s267 = s201*s218*s223;
				const Real s268 = s191*s215*s218;
				const Real s269 = s267 + s268;
				const Real s270 = -(s203*s210*s225*s233);
				const Real s271 = -(s194*s223*s225*s233);
				const Real s272 = s194*s218*s236;
				const Real s273 = -2*s188*s203*s218;
				const Real s274 = s270 + s271 + s272 + s273;
				const Real s275 = -(s203*s223*s225*s233);
				const Real s276 = -(s194*s215*s225*s233);
				const Real s277 = -2*s194*s199*s218;
				const Real s278 = s203*s218*s236;
				const Real s279 = s275 + s276 + s277 + s278;
				const Real s280 = s203*s210*s218;
				const Real s281 = s194*s218*s223;
				const Real s282 = s280 + s281;
				const Real s283 = s203*s218*s223;
				const Real s284 = s194*s215*s218;
				const Real s285 = s283 + s284;
				const Real s286 = P_D_far__[15*i+0];
				const Real s287 = P_D_far__[15*i+1];
				const Real s288 = s2 + s3 + s9;
				const Real s289 = 2*s19*s8;
				const Real s290 = -2*s0*s2*s3;
				const Real s291 = -2*s0*s10*s4;
				const Real s292 = 2*s5*s8;
				const Real s293 = -2*s0*s11*s6;
				const Real s294 = 2*s7*s8;
				const Real s295 = -4*s2*s8*s9;
				const Real s296 = 2*s0*s3*s9;
				const Real s297 = 2*s0*s2*s9;
				const Real s298 = 2*s53*s8;
				const Real s299 = -2*s0*s53;
				const Real s300 = 2*s12*s2*s3;
				const Real s301 = -2*s12*s19;
				const Real s302 = 2*s10*s12*s4;
				const Real s303 = -2*s12*s5;
				const Real s304 = 2*s11*s12*s6;
				const Real s305 = -2*s12*s7;
				const Real s306 = -2*s12*s3*s9;
				const Real s307 = 2*s12*s2*s9;
				const Real s308 = 2*s0*s10*s14;
				const Real s309 = -4*s14*s4*s8;
				const Real s310 = 2*s0*s14*s4;
				const Real s311 = -2*s10*s12*s14;
				const Real s312 = 2*s12*s14*s4;
				const Real s313 = 2*s15*s8;
				const Real s314 = -2*s0*s15;
				const Real s315 = 2*s0*s11*s16;
				const Real s316 = -4*s16*s6*s8;
				const Real s317 = 2*s0*s16*s6;
				const Real s318 = -2*s11*s12*s16;
				const Real s319 = 2*s12*s16*s6;
				const Real s320 = 2*s17*s8;
				const Real s321 = -2*s0*s17;
				const Real s322 = s289 + s290 + s291 + s292 + s293 + s294 + s295 + s296 + s297 + s298 + s299 + s300 + s301 + s302 + s303 + s304 + s305 + s306 + s307 + s308 + s309 + s310 + s311 + s312 + s313 + s314 + s315 + s316 + s317 + s318 + s319 + s320 + s321;
				const Real s323 = P_D_far__[15*i+3];
				const Real s324 = s10 + s14 + s4;
				const Real s325 = P_D_far__[15*i+4];
				const Real s326 = s11 + s16 + s6;
				const Real s327 = P_D_far__[15*i+2];
				const Real s328 = s0 + s12 + s8;
				const Real s329 = s151/6.;
				const Real s330 = P_D_far__[15*i+5];
				const Real s331 = -2*s201*s210;
				const Real s332 = 2*s8;
				const Real s333 = -s0;
				const Real s334 = -s12;
				const Real s335 = s332 + s333 + s334;
				const Real s336 = -2*s207*s335;
				const Real s337 = -2*s191*s215;
				const Real s338 = s331 + s336 + s337;
				const Real s339 = -2*s8;
				const Real s340 = s0 + s12 + s339;
				const Real s341 = s199*s250;
				const Real s342 = s188*s253;
				const Real s343 = s341 + s342;
				const Real s344 = P_D_far__[15*i+6];
				const Real s345 = -(s199*s210*s218);
				const Real s346 = -(s188*s218*s223);
				const Real s347 = -(s199*s218*s223);
				const Real s348 = -(s188*s215*s218);
				const Real s349 = -(s199*s210*s225*s338);
				const Real s350 = -(s188*s223*s225*s338);
				const Real s351 = s188*s218*s340;
				const Real s352 = s261 + s349 + s350 + s351;
				const Real s353 = -(s199*s223*s225*s338);
				const Real s354 = -(s188*s215*s225*s338);
				const Real s355 = s199*s218*s340;
				const Real s356 = s257 + s353 + s354 + s355;
				const Real s357 = s201*s250;
				const Real s358 = s191*s253;
				const Real s359 = s357 + s358;
				const Real s360 = P_D_far__[15*i+7];
				const Real s361 = s203*s250;
				const Real s362 = s194*s253;
				const Real s363 = s361 + s362;
				const Real s364 = P_D_far__[15*i+8];
				const Real s365 = s205*s250;
				const Real s366 = s197*s253;
				const Real s367 = s365 + s366;
				const Real s368 = P_D_far__[15*i+9];
				const Real s369 = -2*s191*s201*s218;
				const Real s370 = s201*s266;
				const Real s371 = s191*s269;
				const Real s372 = s370 + s371;
				const Real s373 = P_D_far__[15*i+10];
				const Real s374 = -(s201*s210*s225*s338);
				const Real s375 = -(s191*s223*s225*s338);
				const Real s376 = s191*s218*s340;
				const Real s377 = s237 + s240 + s369 + s374 + s375 + s376;
				const Real s378 = -(s201*s223*s225*s338);
				const Real s379 = -(s191*s215*s225*s338);
				const Real s380 = s201*s218*s340;
				const Real s381 = s237 + s246 + s369 + s378 + s379 + s380;
				const Real s382 = s203*s266;
				const Real s383 = s194*s269;
				const Real s384 = s382 + s383;
				const Real s385 = P_D_far__[15*i+11];
				const Real s386 = s205*s266;
				const Real s387 = s197*s269;
				const Real s388 = s386 + s387;
				const Real s389 = P_D_far__[15*i+12];
				const Real s390 = s203*s282;
				const Real s391 = s194*s285;
				const Real s392 = s390 + s391;
				const Real s393 = P_D_far__[15*i+13];
				const Real s394 = -(s203*s210*s225*s338);
				const Real s395 = -(s194*s223*s225*s338);
				const Real s396 = s194*s218*s340;
				const Real s397 = -2*s191*s203*s218;
				const Real s398 = s394 + s395 + s396 + s397;
				const Real s399 = -(s203*s223*s225*s338);
				const Real s400 = -(s194*s215*s225*s338);
				const Real s401 = -2*s194*s201*s218;
				const Real s402 = s203*s218*s340;
				const Real s403 = s399 + s400 + s401 + s402;
				const Real s404 = s205*s282;
				const Real s405 = s197*s285;
				const Real s406 = s404 + s405;
				const Real s407 = P_D_far__[15*i+14];
				const Real s408 = s205*s210*s218;
				const Real s409 = s197*s218*s223;
				const Real s410 = s408 + s409;
				const Real s411 = s205*s410;
				const Real s412 = s205*s218*s223;
				const Real s413 = s197*s215*s218;
				const Real s414 = s412 + s413;
				const Real s415 = s197*s414;
				const Real s416 = s411 + s415;
				const Real s417 = 2*s10*s19;
				const Real s418 = 2*s1*s10;
				const Real s419 = -2*s2*s3*s4;
				const Real s420 = -2*s0*s4*s8;
				const Real s421 = -2*s11*s4*s6;
				const Real s422 = 2*s10*s7;
				const Real s423 = -4*s10*s2*s9;
				const Real s424 = 2*s3*s4*s9;
				const Real s425 = 2*s2*s4*s9;
				const Real s426 = 2*s10*s53;
				const Real s427 = -2*s4*s53;
				const Real s428 = -4*s0*s10*s12;
				const Real s429 = 2*s12*s4*s8;
				const Real s430 = 2*s0*s12*s4;
				const Real s431 = 2*s10*s13;
				const Real s432 = -2*s13*s4;
				const Real s433 = 2*s14*s2*s3;
				const Real s434 = -2*s14*s19;
				const Real s435 = 2*s0*s14*s8;
				const Real s436 = -2*s1*s14;
				const Real s437 = 2*s11*s14*s6;
				const Real s438 = -2*s14*s7;
				const Real s439 = -2*s14*s3*s9;
				const Real s440 = 2*s14*s2*s9;
				const Real s441 = -2*s12*s14*s8;
				const Real s442 = 2*s0*s12*s14;
				const Real s443 = 2*s11*s16*s4;
				const Real s444 = -4*s10*s16*s6;
				const Real s445 = 2*s16*s4*s6;
				const Real s446 = -2*s11*s14*s16;
				const Real s447 = 2*s14*s16*s6;
				const Real s448 = 2*s10*s17;
				const Real s449 = -2*s17*s4;
				const Real s450 = s417 + s418 + s419 + s420 + s421 + s422 + s423 + s424 + s425 + s426 + s427 + s428 + s429 + s430 + s431 + s432 + s433 + s434 + s435 + s436 + s437 + s438 + s439 + s440 + s441 + s442 + s443 + s444 + s445 + s446 + s447 + s448 + s449;
				const Real s451 = -2*s203*s210;
				const Real s452 = 2*s10;
				const Real s453 = -s4;
				const Real s454 = -s14;
				const Real s455 = s452 + s453 + s454;
				const Real s456 = -2*s207*s455;
				const Real s457 = -2*s194*s215;
				const Real s458 = s451 + s456 + s457;
				const Real s459 = -2*s10;
				const Real s460 = s14 + s4 + s459;
				const Real s461 = -(s199*s210*s225*s458);
				const Real s462 = -(s188*s223*s225*s458);
				const Real s463 = s188*s218*s460;
				const Real s464 = s277 + s461 + s462 + s463;
				const Real s465 = -(s199*s223*s225*s458);
				const Real s466 = -(s188*s215*s225*s458);
				const Real s467 = s199*s218*s460;
				const Real s468 = s273 + s465 + s466 + s467;
				const Real s469 = -(s201*s210*s218);
				const Real s470 = -(s191*s218*s223);
				const Real s471 = -(s201*s218*s223);
				const Real s472 = -(s191*s215*s218);
				const Real s473 = -(s201*s210*s225*s458);
				const Real s474 = -(s191*s223*s225*s458);
				const Real s475 = s191*s218*s460;
				const Real s476 = s401 + s473 + s474 + s475;
				const Real s477 = -(s201*s223*s225*s458);
				const Real s478 = -(s191*s215*s225*s458);
				const Real s479 = s201*s218*s460;
				const Real s480 = s397 + s477 + s478 + s479;
				const Real s481 = -2*s194*s203*s218;
				const Real s482 = -(s203*s210*s225*s458);
				const Real s483 = -(s194*s223*s225*s458);
				const Real s484 = s194*s218*s460;
				const Real s485 = s237 + s240 + s481 + s482 + s483 + s484;
				const Real s486 = -(s203*s223*s225*s458);
				const Real s487 = -(s194*s215*s225*s458);
				const Real s488 = s203*s218*s460;
				const Real s489 = s237 + s246 + s481 + s486 + s487 + s488;
				const Real s490 = 2*s11*s19;
				const Real s491 = 2*s1*s11;
				const Real s492 = 2*s11*s5;
				const Real s493 = -2*s2*s3*s6;
				const Real s494 = -2*s0*s6*s8;
				const Real s495 = -2*s10*s4*s6;
				const Real s496 = -4*s11*s2*s9;
				const Real s497 = 2*s3*s6*s9;
				const Real s498 = 2*s2*s6*s9;
				const Real s499 = 2*s11*s53;
				const Real s500 = -2*s53*s6;
				const Real s501 = -4*s0*s11*s12;
				const Real s502 = 2*s12*s6*s8;
				const Real s503 = 2*s0*s12*s6;
				const Real s504 = 2*s11*s13;
				const Real s505 = -2*s13*s6;
				const Real s506 = -4*s11*s14*s4;
				const Real s507 = 2*s10*s14*s6;
				const Real s508 = 2*s14*s4*s6;
				const Real s509 = 2*s11*s15;
				const Real s510 = -2*s15*s6;
				const Real s511 = 2*s16*s2*s3;
				const Real s512 = -2*s16*s19;
				const Real s513 = 2*s0*s16*s8;
				const Real s514 = -2*s1*s16;
				const Real s515 = 2*s10*s16*s4;
				const Real s516 = -2*s16*s5;
				const Real s517 = -2*s16*s3*s9;
				const Real s518 = 2*s16*s2*s9;
				const Real s519 = -2*s12*s16*s8;
				const Real s520 = 2*s0*s12*s16;
				const Real s521 = -2*s10*s14*s16;
				const Real s522 = 2*s14*s16*s4;
				const Real s523 = s490 + s491 + s492 + s493 + s494 + s495 + s496 + s497 + s498 + s499 + s500 + s501 + s502 + s503 + s504 + s505 + s506 + s507 + s508 + s509 + s510 + s511 + s512 + s513 + s514 + s515 + s516 + s517 + s518 + s519 + s520 + s521 + s522;
				const Real s524 = -2*s205*s210;
				const Real s525 = 2*s11;
				const Real s526 = -s6;
				const Real s527 = -s16;
				const Real s528 = s525 + s526 + s527;
				const Real s529 = -2*s207*s528;
				const Real s530 = -2*s197*s215;
				const Real s531 = s524 + s529 + s530;
				const Real s532 = -2*s197*s199*s218;
				const Real s533 = -2*s188*s205*s218;
				const Real s534 = -2*s11;
				const Real s535 = s16 + s534 + s6;
				const Real s536 = -(s199*s210*s225*s531);
				const Real s537 = -(s188*s223*s225*s531);
				const Real s538 = s188*s218*s535;
				const Real s539 = s532 + s536 + s537 + s538;
				const Real s540 = -(s199*s223*s225*s531);
				const Real s541 = -(s188*s215*s225*s531);
				const Real s542 = s199*s218*s535;
				const Real s543 = s533 + s540 + s541 + s542;
				const Real s544 = -2*s197*s201*s218;
				const Real s545 = -2*s191*s205*s218;
				const Real s546 = -(s201*s210*s225*s531);
				const Real s547 = -(s191*s223*s225*s531);
				const Real s548 = s191*s218*s535;
				const Real s549 = s544 + s546 + s547 + s548;
				const Real s550 = -(s201*s223*s225*s531);
				const Real s551 = -(s191*s215*s225*s531);
				const Real s552 = s201*s218*s535;
				const Real s553 = s545 + s550 + s551 + s552;
				const Real s554 = -2*s197*s203*s218;
				const Real s555 = -2*s194*s205*s218;
				const Real s556 = -(s203*s210*s218);
				const Real s557 = -(s194*s218*s223);
				const Real s558 = -(s203*s218*s223);
				const Real s559 = -(s194*s215*s218);
				const Real s560 = -(s203*s210*s225*s531);
				const Real s561 = -(s194*s223*s225*s531);
				const Real s562 = s194*s218*s535;
				const Real s563 = s554 + s560 + s561 + s562;
				const Real s564 = -(s203*s223*s225*s531);
				const Real s565 = -(s194*s215*s225*s531);
				const Real s566 = s203*s218*s535;
				const Real s567 = s555 + s564 + s565 + s566;
				const Real s568 = -2*s197*s205*s218;
				const Real s569 = 2*s18*s2;
				const Real s570 = 2*s2*s21;
				const Real s571 = 2*s2*s23;
				const Real s572 = -2*s0*s3*s8;
				const Real s573 = -2*s10*s3*s4;
				const Real s574 = -2*s11*s3*s6;
				const Real s575 = -2*s18*s9;
				const Real s576 = -2*s21*s9;
				const Real s577 = -2*s23*s9;
				const Real s578 = 2*s12*s3*s8;
				const Real s579 = -4*s12*s2*s8;
				const Real s580 = 2*s0*s12*s3;
				const Real s581 = 2*s12*s8*s9;
				const Real s582 = -2*s0*s12*s9;
				const Real s583 = -2*s13*s3;
				const Real s584 = 2*s13*s2;
				const Real s585 = 2*s10*s14*s3;
				const Real s586 = -4*s10*s14*s2;
				const Real s587 = 2*s14*s3*s4;
				const Real s588 = 2*s10*s14*s9;
				const Real s589 = -2*s14*s4*s9;
				const Real s590 = -2*s15*s3;
				const Real s591 = 2*s15*s2;
				const Real s592 = 2*s11*s16*s3;
				const Real s593 = -4*s11*s16*s2;
				const Real s594 = 2*s16*s3*s6;
				const Real s595 = 2*s11*s16*s9;
				const Real s596 = -2*s16*s6*s9;
				const Real s597 = -2*s17*s3;
				const Real s598 = 2*s17*s2;
				const Real s599 = s159 + s161 + s163 + s569 + s570 + s571 + s572 + s573 + s574 + s575 + s576 + s577 + s578 + s579 + s580 + s581 + s582 + s583 + s584 + s585 + s586 + s587 + s588 + s589 + s590 + s591 + s592 + s593 + s594 + s595 + s596 + s597 + s598;
				const Real s600 = -2*s199*s207;
				const Real s601 = 2*s188*s215;
				const Real s602 = s600 + s601;
				const Real s603 = -(s199*s210*s225*s602);
				const Real s604 = -(s188*s223*s225*s602);
				const Real s605 = s188*s199*s218;
				const Real s606 = s218*s223;
				const Real s607 = s603 + s604 + s605 + s606;
				const Real s608 = -(s199*s223*s225*s602);
				const Real s609 = -(s188*s215*s225*s602);
				const Real s610 = -(s211*s218);
				const Real s611 = s215*s218;
				const Real s612 = s608 + s609 + s610 + s611;
				const Real s613 = -(s201*s210*s225*s602);
				const Real s614 = -(s191*s223*s225*s602);
				const Real s615 = -(s191*s199*s218);
				const Real s616 = 2*s188*s201*s218;
				const Real s617 = s613 + s614 + s615 + s616;
				const Real s618 = -(s201*s223*s225*s602);
				const Real s619 = -(s191*s215*s225*s602);
				const Real s620 = -(s199*s201*s218);
				const Real s621 = s618 + s619 + s620;
				const Real s622 = -(s203*s210*s225*s602);
				const Real s623 = -(s194*s223*s225*s602);
				const Real s624 = -(s194*s199*s218);
				const Real s625 = 2*s188*s203*s218;
				const Real s626 = s622 + s623 + s624 + s625;
				const Real s627 = -(s203*s223*s225*s602);
				const Real s628 = -(s194*s215*s225*s602);
				const Real s629 = -(s199*s203*s218);
				const Real s630 = s627 + s628 + s629;
				const Real s631 = -2*s2*s3*s8;
				const Real s632 = 2*s0*s26;
				const Real s633 = 2*s0*s21;
				const Real s634 = 2*s0*s23;
				const Real s635 = -2*s10*s4*s8;
				const Real s636 = -2*s11*s6*s8;
				const Real s637 = 2*s3*s8*s9;
				const Real s638 = 2*s2*s8*s9;
				const Real s639 = -4*s0*s3*s9;
				const Real s640 = -2*s53*s8;
				const Real s641 = 2*s0*s53;
				const Real s642 = -2*s12*s26;
				const Real s643 = -2*s12*s21;
				const Real s644 = -2*s12*s23;
				const Real s645 = 2*s12*s3*s9;
				const Real s646 = -2*s12*s2*s9;
				const Real s647 = 2*s10*s14*s8;
				const Real s648 = -4*s0*s10*s14;
				const Real s649 = 2*s14*s4*s8;
				const Real s650 = 2*s10*s12*s14;
				const Real s651 = -2*s12*s14*s4;
				const Real s652 = -2*s15*s8;
				const Real s653 = 2*s0*s15;
				const Real s654 = 2*s11*s16*s8;
				const Real s655 = -4*s0*s11*s16;
				const Real s656 = 2*s16*s6*s8;
				const Real s657 = 2*s11*s12*s16;
				const Real s658 = -2*s12*s16*s6;
				const Real s659 = -2*s17*s8;
				const Real s660 = 2*s0*s17;
				const Real s661 = s300 + s302 + s304 + s631 + s632 + s633 + s634 + s635 + s636 + s637 + s638 + s639 + s640 + s641 + s642 + s643 + s644 + s645 + s646 + s647 + s648 + s649 + s650 + s651 + s652 + s653 + s654 + s655 + s656 + s657 + s658 + s659 + s660;
				const Real s662 = -2*s201*s207;
				const Real s663 = 2*s191*s215;
				const Real s664 = s662 + s663;
				const Real s665 = -(s199*s210*s225*s664);
				const Real s666 = -(s188*s223*s225*s664);
				const Real s667 = 2*s191*s199*s218;
				const Real s668 = -(s188*s201*s218);
				const Real s669 = s665 + s666 + s667 + s668;
				const Real s670 = -(s199*s223*s225*s664);
				const Real s671 = -(s188*s215*s225*s664);
				const Real s672 = s620 + s670 + s671;
				const Real s673 = -(s201*s210*s225*s664);
				const Real s674 = -(s191*s223*s225*s664);
				const Real s675 = s191*s201*s218;
				const Real s676 = s606 + s673 + s674 + s675;
				const Real s677 = -(s201*s223*s225*s664);
				const Real s678 = -(s191*s215*s225*s664);
				const Real s679 = -(s212*s218);
				const Real s680 = s611 + s677 + s678 + s679;
				const Real s681 = -(s203*s210*s225*s664);
				const Real s682 = -(s194*s223*s225*s664);
				const Real s683 = -(s194*s201*s218);
				const Real s684 = 2*s191*s203*s218;
				const Real s685 = s681 + s682 + s683 + s684;
				const Real s686 = -(s203*s223*s225*s664);
				const Real s687 = -(s194*s215*s225*s664);
				const Real s688 = -(s201*s203*s218);
				const Real s689 = s686 + s687 + s688;
				const Real s690 = -2*s10*s2*s3;
				const Real s691 = -2*s0*s10*s8;
				const Real s692 = 2*s26*s4;
				const Real s693 = 2*s18*s4;
				const Real s694 = 2*s23*s4;
				const Real s695 = -2*s10*s11*s6;
				const Real s696 = 2*s10*s3*s9;
				const Real s697 = 2*s10*s2*s9;
				const Real s698 = -4*s3*s4*s9;
				const Real s699 = -2*s10*s53;
				const Real s700 = 2*s4*s53;
				const Real s701 = 2*s10*s12*s8;
				const Real s702 = 2*s0*s10*s12;
				const Real s703 = -4*s12*s4*s8;
				const Real s704 = -2*s10*s13;
				const Real s705 = 2*s13*s4;
				const Real s706 = -2*s14*s26;
				const Real s707 = -2*s14*s18;
				const Real s708 = -2*s14*s23;
				const Real s709 = 2*s14*s3*s9;
				const Real s710 = -2*s14*s2*s9;
				const Real s711 = 2*s12*s14*s8;
				const Real s712 = -2*s0*s12*s14;
				const Real s713 = 2*s10*s11*s16;
				const Real s714 = -4*s11*s16*s4;
				const Real s715 = 2*s10*s16*s6;
				const Real s716 = 2*s11*s14*s16;
				const Real s717 = -2*s14*s16*s6;
				const Real s718 = -2*s10*s17;
				const Real s719 = 2*s17*s4;
				const Real s720 = s433 + s435 + s437 + s690 + s691 + s692 + s693 + s694 + s695 + s696 + s697 + s698 + s699 + s700 + s701 + s702 + s703 + s704 + s705 + s706 + s707 + s708 + s709 + s710 + s711 + s712 + s713 + s714 + s715 + s716 + s717 + s718 + s719;
				const Real s721 = -2*s203*s207;
				const Real s722 = 2*s194*s215;
				const Real s723 = s721 + s722;
				const Real s724 = -(s199*s210*s225*s723);
				const Real s725 = -(s188*s223*s225*s723);
				const Real s726 = 2*s194*s199*s218;
				const Real s727 = -(s188*s203*s218);
				const Real s728 = s724 + s725 + s726 + s727;
				const Real s729 = -(s199*s223*s225*s723);
				const Real s730 = -(s188*s215*s225*s723);
				const Real s731 = s629 + s729 + s730;
				const Real s732 = -(s201*s210*s225*s723);
				const Real s733 = -(s191*s223*s225*s723);
				const Real s734 = 2*s194*s201*s218;
				const Real s735 = -(s191*s203*s218);
				const Real s736 = s732 + s733 + s734 + s735;
				const Real s737 = -(s201*s223*s225*s723);
				const Real s738 = -(s191*s215*s225*s723);
				const Real s739 = s688 + s737 + s738;
				const Real s740 = -(s203*s210*s225*s723);
				const Real s741 = -(s194*s223*s225*s723);
				const Real s742 = s194*s203*s218;
				const Real s743 = s606 + s740 + s741 + s742;
				const Real s744 = -(s203*s223*s225*s723);
				const Real s745 = -(s194*s215*s225*s723);
				const Real s746 = -(s213*s218);
				const Real s747 = s611 + s744 + s745 + s746;
				const Real s748 = -2*s11*s2*s3;
				const Real s749 = -2*s0*s11*s8;
				const Real s750 = -2*s10*s11*s4;
				const Real s751 = 2*s26*s6;
				const Real s752 = 2*s18*s6;
				const Real s753 = 2*s21*s6;
				const Real s754 = 2*s11*s3*s9;
				const Real s755 = 2*s11*s2*s9;
				const Real s756 = -4*s3*s6*s9;
				const Real s757 = -2*s11*s53;
				const Real s758 = 2*s53*s6;
				const Real s759 = 2*s11*s12*s8;
				const Real s760 = 2*s0*s11*s12;
				const Real s761 = -4*s12*s6*s8;
				const Real s762 = -2*s11*s13;
				const Real s763 = 2*s13*s6;
				const Real s764 = 2*s10*s11*s14;
				const Real s765 = 2*s11*s14*s4;
				const Real s766 = -4*s10*s14*s6;
				const Real s767 = -2*s11*s15;
				const Real s768 = 2*s15*s6;
				const Real s769 = -2*s16*s26;
				const Real s770 = -2*s16*s18;
				const Real s771 = -2*s16*s21;
				const Real s772 = 2*s16*s3*s9;
				const Real s773 = -2*s16*s2*s9;
				const Real s774 = 2*s12*s16*s8;
				const Real s775 = -2*s0*s12*s16;
				const Real s776 = 2*s10*s14*s16;
				const Real s777 = -2*s14*s16*s4;
				const Real s778 = s511 + s513 + s515 + s748 + s749 + s750 + s751 + s752 + s753 + s754 + s755 + s756 + s757 + s758 + s759 + s760 + s761 + s762 + s763 + s764 + s765 + s766 + s767 + s768 + s769 + s770 + s771 + s772 + s773 + s774 + s775 + s776 + s777;
				const Real s779 = -2*s205*s207;
				const Real s780 = 2*s197*s215;
				const Real s781 = s779 + s780;
				const Real s782 = -(s199*s205*s218);
				const Real s783 = -(s199*s210*s225*s781);
				const Real s784 = -(s188*s223*s225*s781);
				const Real s785 = 2*s197*s199*s218;
				const Real s786 = -(s188*s205*s218);
				const Real s787 = s783 + s784 + s785 + s786;
				const Real s788 = -(s199*s223*s225*s781);
				const Real s789 = -(s188*s215*s225*s781);
				const Real s790 = s782 + s788 + s789;
				const Real s791 = -(s201*s205*s218);
				const Real s792 = -(s201*s210*s225*s781);
				const Real s793 = -(s191*s223*s225*s781);
				const Real s794 = 2*s197*s201*s218;
				const Real s795 = -(s191*s205*s218);
				const Real s796 = s792 + s793 + s794 + s795;
				const Real s797 = -(s201*s223*s225*s781);
				const Real s798 = -(s191*s215*s225*s781);
				const Real s799 = s791 + s797 + s798;
				const Real s800 = -(s203*s205*s218);
				const Real s801 = -(s203*s210*s225*s781);
				const Real s802 = -(s194*s223*s225*s781);
				const Real s803 = 2*s197*s203*s218;
				const Real s804 = -(s194*s205*s218);
				const Real s805 = s801 + s802 + s803 + s804;
				const Real s806 = -(s203*s223*s225*s781);
				const Real s807 = -(s194*s215*s225*s781);
				const Real s808 = s800 + s806 + s807;
				const Real s809 = -2*s18*s2;
				const Real s810 = -2*s2*s21;
				const Real s811 = -2*s2*s23;
				const Real s812 = 2*s0*s3*s8;
				const Real s813 = 2*s0*s2*s8;
				const Real s814 = -2*s1*s3;
				const Real s815 = 2*s10*s3*s4;
				const Real s816 = 2*s10*s2*s4;
				const Real s817 = -2*s3*s5;
				const Real s818 = 2*s11*s3*s6;
				const Real s819 = 2*s11*s2*s6;
				const Real s820 = -2*s3*s7;
				const Real s821 = 2*s18*s9;
				const Real s822 = 2*s21*s9;
				const Real s823 = 2*s23*s9;
				const Real s824 = -4*s0*s8*s9;
				const Real s825 = 2*s1*s9;
				const Real s826 = -4*s10*s4*s9;
				const Real s827 = 2*s5*s9;
				const Real s828 = -4*s11*s6*s9;
				const Real s829 = 2*s7*s9;
				const Real s830 = -2*s12*s3*s8;
				const Real s831 = -2*s0*s12*s2;
				const Real s832 = -2*s10*s14*s3;
				const Real s833 = -2*s14*s2*s4;
				const Real s834 = -2*s11*s16*s3;
				const Real s835 = -2*s16*s2*s6;
				const Real s836 = s165 + s172 + s179 + s580 + s587 + s594 + s809 + s810 + s811 + s812 + s813 + s814 + s815 + s816 + s817 + s818 + s819 + s820 + s821 + s822 + s823 + s824 + s825 + s826 + s827 + s828 + s829 + s830 + s831 + s832 + s833 + s834 + s835;
				const Real s837 = 2*s199*s210;
				const Real s838 = -2*s188*s207;
				const Real s839 = s837 + s838;
				const Real s840 = -(s199*s210*s225*s839);
				const Real s841 = -(s188*s223*s225*s839);
				const Real s842 = -(s189*s218);
				const Real s843 = s210*s218;
				const Real s844 = s840 + s841 + s842 + s843;
				const Real s845 = -(s199*s223*s225*s839);
				const Real s846 = -(s188*s215*s225*s839);
				const Real s847 = s605 + s606 + s845 + s846;
				const Real s848 = -(s201*s210*s225*s839);
				const Real s849 = -(s191*s223*s225*s839);
				const Real s850 = -(s188*s191*s218);
				const Real s851 = s848 + s849 + s850;
				const Real s852 = -(s201*s223*s225*s839);
				const Real s853 = -(s191*s215*s225*s839);
				const Real s854 = s667 + s668 + s852 + s853;
				const Real s855 = -(s203*s210*s225*s839);
				const Real s856 = -(s194*s223*s225*s839);
				const Real s857 = -(s188*s194*s218);
				const Real s858 = s855 + s856 + s857;
				const Real s859 = -(s203*s223*s225*s839);
				const Real s860 = -(s194*s215*s225*s839);
				const Real s861 = s726 + s727 + s859 + s860;
				const Real s862 = 2*s2*s3*s8;
				const Real s863 = -2*s19*s8;
				const Real s864 = -2*s0*s26;
				const Real s865 = -2*s0*s21;
				const Real s866 = -2*s0*s23;
				const Real s867 = 2*s0*s2*s3;
				const Real s868 = 2*s10*s4*s8;
				const Real s869 = 2*s0*s10*s4;
				const Real s870 = -2*s5*s8;
				const Real s871 = 2*s11*s6*s8;
				const Real s872 = 2*s0*s11*s6;
				const Real s873 = -2*s7*s8;
				const Real s874 = -2*s3*s8*s9;
				const Real s875 = -2*s0*s2*s9;
				const Real s876 = 2*s12*s26;
				const Real s877 = 2*s12*s21;
				const Real s878 = 2*s12*s23;
				const Real s879 = -4*s12*s2*s3;
				const Real s880 = 2*s12*s19;
				const Real s881 = -4*s10*s12*s4;
				const Real s882 = 2*s12*s5;
				const Real s883 = -4*s11*s12*s6;
				const Real s884 = 2*s12*s7;
				const Real s885 = -2*s10*s14*s8;
				const Real s886 = -2*s0*s14*s4;
				const Real s887 = -2*s11*s16*s8;
				const Real s888 = -2*s0*s16*s6;
				const Real s889 = s296 + s308 + s315 + s638 + s649 + s656 + s862 + s863 + s864 + s865 + s866 + s867 + s868 + s869 + s870 + s871 + s872 + s873 + s874 + s875 + s876 + s877 + s878 + s879 + s880 + s881 + s882 + s883 + s884 + s885 + s886 + s887 + s888;
				const Real s890 = 2*s201*s210;
				const Real s891 = -2*s191*s207;
				const Real s892 = s890 + s891;
				const Real s893 = -(s199*s210*s225*s892);
				const Real s894 = -(s188*s223*s225*s892);
				const Real s895 = s850 + s893 + s894;
				const Real s896 = -(s199*s223*s225*s892);
				const Real s897 = -(s188*s215*s225*s892);
				const Real s898 = s615 + s616 + s896 + s897;
				const Real s899 = -(s201*s210*s225*s892);
				const Real s900 = -(s191*s223*s225*s892);
				const Real s901 = -(s192*s218);
				const Real s902 = s843 + s899 + s900 + s901;
				const Real s903 = -(s201*s223*s225*s892);
				const Real s904 = -(s191*s215*s225*s892);
				const Real s905 = s606 + s675 + s903 + s904;
				const Real s906 = -(s203*s210*s225*s892);
				const Real s907 = -(s194*s223*s225*s892);
				const Real s908 = -(s191*s194*s218);
				const Real s909 = s906 + s907 + s908;
				const Real s910 = -(s203*s223*s225*s892);
				const Real s911 = -(s194*s215*s225*s892);
				const Real s912 = s734 + s735 + s910 + s911;
				const Real s913 = 2*s10*s2*s3;
				const Real s914 = -2*s10*s19;
				const Real s915 = 2*s0*s10*s8;
				const Real s916 = -2*s1*s10;
				const Real s917 = -2*s26*s4;
				const Real s918 = -2*s18*s4;
				const Real s919 = -2*s23*s4;
				const Real s920 = 2*s2*s3*s4;
				const Real s921 = 2*s0*s4*s8;
				const Real s922 = 2*s10*s11*s6;
				const Real s923 = 2*s11*s4*s6;
				const Real s924 = -2*s10*s7;
				const Real s925 = -2*s10*s3*s9;
				const Real s926 = -2*s2*s4*s9;
				const Real s927 = -2*s10*s12*s8;
				const Real s928 = -2*s0*s12*s4;
				const Real s929 = 2*s14*s26;
				const Real s930 = 2*s14*s18;
				const Real s931 = 2*s14*s23;
				const Real s932 = -4*s14*s2*s3;
				const Real s933 = 2*s14*s19;
				const Real s934 = -4*s0*s14*s8;
				const Real s935 = 2*s1*s14;
				const Real s936 = -4*s11*s14*s6;
				const Real s937 = 2*s14*s7;
				const Real s938 = -2*s10*s11*s16;
				const Real s939 = -2*s16*s4*s6;
				const Real s940 = s424 + s429 + s443 + s697 + s702 + s715 + s913 + s914 + s915 + s916 + s917 + s918 + s919 + s920 + s921 + s922 + s923 + s924 + s925 + s926 + s927 + s928 + s929 + s930 + s931 + s932 + s933 + s934 + s935 + s936 + s937 + s938 + s939;
				const Real s941 = 2*s203*s210;
				const Real s942 = -2*s194*s207;
				const Real s943 = s941 + s942;
				const Real s944 = -(s199*s210*s225*s943);
				const Real s945 = -(s188*s223*s225*s943);
				const Real s946 = s857 + s944 + s945;
				const Real s947 = -(s199*s223*s225*s943);
				const Real s948 = -(s188*s215*s225*s943);
				const Real s949 = s624 + s625 + s947 + s948;
				const Real s950 = -(s201*s210*s225*s943);
				const Real s951 = -(s191*s223*s225*s943);
				const Real s952 = s908 + s950 + s951;
				const Real s953 = -(s201*s223*s225*s943);
				const Real s954 = -(s191*s215*s225*s943);
				const Real s955 = s683 + s684 + s953 + s954;
				const Real s956 = -(s203*s210*s225*s943);
				const Real s957 = -(s194*s223*s225*s943);
				const Real s958 = -(s195*s218);
				const Real s959 = s843 + s956 + s957 + s958;
				const Real s960 = -(s203*s223*s225*s943);
				const Real s961 = -(s194*s215*s225*s943);
				const Real s962 = s606 + s742 + s960 + s961;
				const Real s963 = 2*s11*s2*s3;
				const Real s964 = -2*s11*s19;
				const Real s965 = 2*s0*s11*s8;
				const Real s966 = -2*s1*s11;
				const Real s967 = 2*s10*s11*s4;
				const Real s968 = -2*s11*s5;
				const Real s969 = -2*s26*s6;
				const Real s970 = -2*s18*s6;
				const Real s971 = -2*s21*s6;
				const Real s972 = 2*s2*s3*s6;
				const Real s973 = 2*s0*s6*s8;
				const Real s974 = 2*s10*s4*s6;
				const Real s975 = -2*s11*s3*s9;
				const Real s976 = -2*s2*s6*s9;
				const Real s977 = -2*s11*s12*s8;
				const Real s978 = -2*s0*s12*s6;
				const Real s979 = -2*s10*s11*s14;
				const Real s980 = -2*s14*s4*s6;
				const Real s981 = 2*s16*s26;
				const Real s982 = 2*s16*s18;
				const Real s983 = 2*s16*s21;
				const Real s984 = -4*s16*s2*s3;
				const Real s985 = 2*s16*s19;
				const Real s986 = -4*s0*s16*s8;
				const Real s987 = 2*s1*s16;
				const Real s988 = -4*s10*s16*s4;
				const Real s989 = 2*s16*s5;
				const Real s990 = s497 + s502 + s507 + s755 + s760 + s765 + s963 + s964 + s965 + s966 + s967 + s968 + s969 + s970 + s971 + s972 + s973 + s974 + s975 + s976 + s977 + s978 + s979 + s980 + s981 + s982 + s983 + s984 + s985 + s986 + s987 + s988 + s989;
				const Real s991 = 2*s205*s210;
				const Real s992 = -2*s197*s207;
				const Real s993 = s991 + s992;
				const Real s994 = -(s188*s197*s218);
				const Real s995 = -(s197*s199*s218);
				const Real s996 = 2*s188*s205*s218;
				const Real s997 = -(s199*s210*s225*s993);
				const Real s998 = -(s188*s223*s225*s993);
				const Real s999 = s994 + s997 + s998;
				const Real s1000 = -(s199*s223*s225*s993);
				const Real s1001 = -(s188*s215*s225*s993);
				const Real s1002 = s1000 + s1001 + s995 + s996;
				const Real s1003 = -(s191*s197*s218);
				const Real s1004 = -(s197*s201*s218);
				const Real s1005 = 2*s191*s205*s218;
				const Real s1006 = -(s201*s210*s225*s993);
				const Real s1007 = -(s191*s223*s225*s993);
				const Real s1008 = s1003 + s1006 + s1007;
				const Real s1009 = -(s201*s223*s225*s993);
				const Real s1010 = -(s191*s215*s225*s993);
				const Real s1011 = s1004 + s1005 + s1009 + s1010;
				const Real s1012 = -(s194*s197*s218);
				const Real s1013 = -(s197*s203*s218);
				const Real s1014 = 2*s194*s205*s218;
				const Real s1015 = -(s203*s210*s225*s993);
				const Real s1016 = -(s194*s223*s225*s993);
				const Real s1017 = s1012 + s1015 + s1016;
				const Real s1018 = -(s203*s223*s225*s993);
				const Real s1019 = -(s194*s215*s225*s993);
				const Real s1020 = s1013 + s1014 + s1018 + s1019;
				const Real s1021 = s197*s205*s218;
				buffer__[12*i+0] += (s152*s186*s286)/4. + (s152*s186*s323*s324)/12. + (s152*s186*s325*s326)/12. + (s152*s186*s327*s328)/12. + s287*((s152*s186*s288)/12. + s329) + s330*((s152*s186*s343)/4. + (s151*(s199*s242 + s188*s247 + s345 + s346 + s347 + s348))/2.) + s344*((s151*(s201*s242 + s191*s247))/2. + (s152*s186*s359)/4.) + s360*((s151*(s203*s242 + s194*s247))/2. + (s152*s186*s363)/4.) + s364*((s151*(s205*s242 + s197*s247))/2. + (s152*s186*s367)/4.) + s368*((s151*(s201*s258 + s191*s263))/2. + (s152*s186*s372)/4.) + s373*((s151*(s203*s258 + s194*s263))/2. + (s152*s186*s384)/4.) + s385*((s151*(s205*s258 + s197*s263))/2. + (s152*s186*s388)/4.) + s389*((s151*(s203*s274 + s194*s279))/2. + (s152*s186*s392)/4.) + s393*((s151*(s205*s274 + s197*s279))/2. + (s152*s186*s406)/4.) + s407*((s152*s186*s416)/4. + (s151*(s197*(-(s197*s215*s225*s233) - s205*s223*s225*s233 + s205*s218*s236 + s532) + s205*(-(s205*s210*s225*s233) - s197*s223*s225*s233 + s197*s218*s236 + s533)))/2.);
				buffer__[12*i+1] += (s152*s286*s322)/4. + (s152*s287*s288*s322)/12. + (s152*s322*s323*s324)/12. + (s152*s322*s325*s326)/12. + s327*((s152*s322*s328)/12. + s329) + s330*((s152*s322*s343)/4. + (s151*(s199*s352 + s188*s356))/2.) + s344*((s151*(s345 + s346 + s347 + s348 + s201*s352 + s191*s356))/2. + (s152*s322*s359)/4.) + s360*((s151*(s203*s352 + s194*s356))/2. + (s152*s322*s363)/4.) + s364*((s151*(s205*s352 + s197*s356))/2. + (s152*s322*s367)/4.) + s373*((s151*(s203*s377 + s194*s381))/2. + (s152*s322*s384)/4.) + s385*((s151*(s205*s377 + s197*s381))/2. + (s152*s322*s388)/4.) + s389*((s152*s322*s392)/4. + (s151*(s203*s398 + s194*s403))/2.) + s393*((s151*(s205*s398 + s197*s403))/2. + (s152*s322*s406)/4.) + s368*((s152*s322*s372)/4. + (s151*(s201*s377 + s191*s381 + s469 + s470 + s471 + s472))/2.) + s407*((s152*s322*s416)/4. + (s151*(s197*(-(s197*s215*s225*s338) - s205*s223*s225*s338 + s205*s218*s340 + s544) + s205*(-(s205*s210*s225*s338) - s197*s223*s225*s338 + s197*s218*s340 + s545)))/2.);
				buffer__[12*i+2] += (s152*s286*s450)/4. + (s152*s287*s288*s450)/12. + (s152*s325*s326*s450)/12. + (s152*s327*s328*s450)/12. + s323*(s329 + (s152*s324*s450)/12.) + s330*((s152*s343*s450)/4. + (s151*(s199*s464 + s188*s468))/2.) + s344*((s152*s359*s450)/4. + (s151*(s201*s464 + s191*s468))/2.) + s360*((s152*s363*s450)/4. + (s151*(s345 + s346 + s347 + s348 + s203*s464 + s194*s468))/2.) + s364*((s152*s367*s450)/4. + (s151*(s205*s464 + s197*s468))/2.) + s368*((s152*s372*s450)/4. + (s151*(s201*s476 + s191*s480))/2.) + s373*((s152*s384*s450)/4. + (s151*(s469 + s470 + s471 + s472 + s203*s476 + s194*s480))/2.) + s385*((s152*s388*s450)/4. + (s151*(s205*s476 + s197*s480))/2.) + s393*((s152*s406*s450)/4. + (s151*(s205*s485 + s197*s489))/2.) + s407*((s152*s416*s450)/4. + (s151*(s197*(-(s197*s215*s225*s458) - s205*s223*s225*s458 + s205*s218*s460 + s554) + s205*(-(s205*s210*s225*s458) - s197*s223*s225*s458 + s197*s218*s460 + s555)))/2.) + s389*((s152*s392*s450)/4. + (s151*(s203*s485 + s194*s489 + s556 + s557 + s558 + s559))/2.);
				buffer__[12*i+3] += (s152*s286*s523)/4. + (s152*s287*s288*s523)/12. + (s152*s323*s324*s523)/12. + (s152*s327*s328*s523)/12. + s325*(s329 + (s152*s326*s523)/12.) + s330*((s152*s343*s523)/4. + (s151*(s199*s539 + s188*s543))/2.) + s344*((s152*s359*s523)/4. + (s151*(s201*s539 + s191*s543))/2.) + s360*((s152*s363*s523)/4. + (s151*(s203*s539 + s194*s543))/2.) + s364*((s152*s367*s523)/4. + (s151*(s345 + s346 + s347 + s348 + s205*s539 + s197*s543))/2.) + s368*((s152*s372*s523)/4. + (s151*(s201*s549 + s191*s553))/2.) + s373*((s152*s384*s523)/4. + (s151*(s203*s549 + s194*s553))/2.) + s385*((s152*s388*s523)/4. + (s151*(s469 + s470 + s471 + s472 + s205*s549 + s197*s553))/2.) + s389*((s152*s392*s523)/4. + (s151*(s203*s563 + s194*s567))/2.) + s393*((s152*s406*s523)/4. + (s151*(s556 + s557 + s558 + s559 + s205*s563 + s197*s567))/2.) + s407*((s152*s416*s523)/4. + (s151*(-(s205*s210*s218) - s197*s215*s218 - s197*s218*s223 - s205*s218*s223 + s205*(s237 + s240 - s205*s210*s225*s531 - s197*s223*s225*s531 + s197*s218*s535 + s568) + s197*(s237 + s246 - s197*s215*s225*s531 - s205*s223*s225*s531 + s205*s218*s535 + s568)))/2.);
				buffer__[12*i+4] += (s152*s286*s599)/4. + (s152*s323*s324*s599)/12. + (s152*s325*s326*s599)/12. + (s152*s327*s328*s599)/12. + s287*(s329 + (s152*s288*s599)/12.) + s330*((s152*s343*s599)/4. + (s151*(s251 + s252 + s199*s607 + s188*s612))/2.) + s344*((s152*s359*s599)/4. + (s151*(s201*s607 + s191*s612))/2.) + s360*((s152*s363*s599)/4. + (s151*(s203*s607 + s194*s612))/2.) + s364*((s152*s367*s599)/4. + (s151*(s205*s607 + s197*s612))/2.) + s368*((s152*s372*s599)/4. + (s151*(s201*s617 + s191*s621))/2.) + s373*((s152*s384*s599)/4. + (s151*(s203*s617 + s194*s621))/2.) + s385*((s152*s388*s599)/4. + (s151*(s205*s617 + s197*s621))/2.) + s389*((s152*s392*s599)/4. + (s151*(s203*s626 + s194*s630))/2.) + s393*((s152*s406*s599)/4. + (s151*(s205*s626 + s197*s630))/2.) + s407*((s152*s416*s599)/4. + (s151*(s197*(-(s197*s215*s225*s602) - s205*s223*s225*s602 + s782) + s205*(-(s205*s210*s225*s602) - s197*s223*s225*s602 + s995 + s996)))/2.);
				buffer__[12*i+5] += (s152*s286*s661)/4. + (s152*s287*s288*s661)/12. + (s152*s323*s324*s661)/12. + (s152*s325*s326*s661)/12. + s327*(s329 + (s152*s328*s661)/12.) + s330*((s152*s343*s661)/4. + (s151*(s199*s669 + s188*s672))/2.) + s344*((s152*s359*s661)/4. + (s151*(s251 + s252 + s201*s669 + s191*s672))/2.) + s360*((s152*s363*s661)/4. + (s151*(s203*s669 + s194*s672))/2.) + s364*((s152*s367*s661)/4. + (s151*(s205*s669 + s197*s672))/2.) + s368*((s152*s372*s661)/4. + (s151*(s267 + s268 + s201*s676 + s191*s680))/2.) + s373*((s152*s384*s661)/4. + (s151*(s203*s676 + s194*s680))/2.) + s385*((s152*s388*s661)/4. + (s151*(s205*s676 + s197*s680))/2.) + s389*((s152*s392*s661)/4. + (s151*(s203*s685 + s194*s689))/2.) + s393*((s152*s406*s661)/4. + (s151*(s205*s685 + s197*s689))/2.) + s407*((s152*s416*s661)/4. + (s151*(s205*(s1004 + s1005 - s205*s210*s225*s664 - s197*s223*s225*s664) + s197*(-(s197*s215*s225*s664) - s205*s223*s225*s664 + s791)))/2.);
				buffer__[12*i+6] += (s152*s286*s720)/4. + (s152*s287*s288*s720)/12. + (s152*s325*s326*s720)/12. + (s152*s327*s328*s720)/12. + s323*(s329 + (s152*s324*s720)/12.) + s330*((s152*s343*s720)/4. + (s151*(s199*s728 + s188*s731))/2.) + s344*((s152*s359*s720)/4. + (s151*(s201*s728 + s191*s731))/2.) + s360*((s152*s363*s720)/4. + (s151*(s251 + s252 + s203*s728 + s194*s731))/2.) + s364*((s152*s367*s720)/4. + (s151*(s205*s728 + s197*s731))/2.) + s368*((s152*s372*s720)/4. + (s151*(s201*s736 + s191*s739))/2.) + s373*((s152*s384*s720)/4. + (s151*(s267 + s268 + s203*s736 + s194*s739))/2.) + s385*((s152*s388*s720)/4. + (s151*(s205*s736 + s197*s739))/2.) + s389*((s152*s392*s720)/4. + (s151*(s283 + s284 + s203*s743 + s194*s747))/2.) + s393*((s152*s406*s720)/4. + (s151*(s205*s743 + s197*s747))/2.) + s407*((s152*s416*s720)/4. + (s151*(s205*(s1013 + s1014 - s205*s210*s225*s723 - s197*s223*s225*s723) + s197*(-(s197*s215*s225*s723) - s205*s223*s225*s723 + s800)))/2.);
				buffer__[12*i+7] += (s152*s286*s778)/4. + (s152*s287*s288*s778)/12. + (s152*s323*s324*s778)/12. + (s152*s327*s328*s778)/12. + s325*(s329 + (s152*s326*s778)/12.) + s407*((s152*s416*s778)/4. + (s151*(s412 + s413 + s205*(s1021 + s606 - s205*s210*s225*s781 - s197*s223*s225*s781) + s197*(-(s214*s218) + s611 - s197*s215*s225*s781 - s205*s223*s225*s781)))/2.) + s330*((s152*s343*s778)/4. + (s151*(s199*s787 + s188*s790))/2.) + s344*((s152*s359*s778)/4. + (s151*(s201*s787 + s191*s790))/2.) + s360*((s152*s363*s778)/4. + (s151*(s203*s787 + s194*s790))/2.) + s364*((s152*s367*s778)/4. + (s151*(s251 + s252 + s205*s787 + s197*s790))/2.) + s368*((s152*s372*s778)/4. + (s151*(s201*s796 + s191*s799))/2.) + s373*((s152*s384*s778)/4. + (s151*(s203*s796 + s194*s799))/2.) + s385*((s152*s388*s778)/4. + (s151*(s267 + s268 + s205*s796 + s197*s799))/2.) + s389*((s152*s392*s778)/4. + (s151*(s203*s805 + s194*s808))/2.) + s393*((s152*s406*s778)/4. + (s151*(s283 + s284 + s205*s805 + s197*s808))/2.);
				buffer__[12*i+8] += (s152*s286*s836)/4. + (s152*s323*s324*s836)/12. + (s152*s325*s326*s836)/12. + (s152*s327*s328*s836)/12. + s287*(s329 + (s152*s288*s836)/12.) + s330*((s152*s343*s836)/4. + (s151*(s248 + s249 + s199*s844 + s188*s847))/2.) + s344*((s152*s359*s836)/4. + (s151*(s201*s844 + s191*s847))/2.) + s360*((s152*s363*s836)/4. + (s151*(s203*s844 + s194*s847))/2.) + s364*((s152*s367*s836)/4. + (s151*(s205*s844 + s197*s847))/2.) + s368*((s152*s372*s836)/4. + (s151*(s201*s851 + s191*s854))/2.) + s373*((s152*s384*s836)/4. + (s151*(s203*s851 + s194*s854))/2.) + s385*((s152*s388*s836)/4. + (s151*(s205*s851 + s197*s854))/2.) + s389*((s152*s392*s836)/4. + (s151*(s203*s858 + s194*s861))/2.) + s393*((s152*s406*s836)/4. + (s151*(s205*s858 + s197*s861))/2.) + s407*((s152*s416*s836)/4. + (s151*(s197*(s785 + s786 - s197*s215*s225*s839 - s205*s223*s225*s839) + s205*(-(s205*s210*s225*s839) - s197*s223*s225*s839 + s994)))/2.);
				buffer__[12*i+9] += (s152*s286*s889)/4. + (s152*s287*s288*s889)/12. + (s152*s323*s324*s889)/12. + (s152*s325*s326*s889)/12. + s327*(s329 + (s152*s328*s889)/12.) + s407*((s152*s416*s889)/4. + (s151*(s205*(s1003 - s205*s210*s225*s892 - s197*s223*s225*s892) + s197*(s794 + s795 - s197*s215*s225*s892 - s205*s223*s225*s892)))/2.) + s330*((s152*s343*s889)/4. + (s151*(s199*s895 + s188*s898))/2.) + s344*((s152*s359*s889)/4. + (s151*(s248 + s249 + s201*s895 + s191*s898))/2.) + s360*((s152*s363*s889)/4. + (s151*(s203*s895 + s194*s898))/2.) + s364*((s152*s367*s889)/4. + (s151*(s205*s895 + s197*s898))/2.) + s368*((s152*s372*s889)/4. + (s151*(s264 + s265 + s201*s902 + s191*s905))/2.) + s373*((s152*s384*s889)/4. + (s151*(s203*s902 + s194*s905))/2.) + s385*((s152*s388*s889)/4. + (s151*(s205*s902 + s197*s905))/2.) + s389*((s152*s392*s889)/4. + (s151*(s203*s909 + s194*s912))/2.) + s393*((s152*s406*s889)/4. + (s151*(s205*s909 + s197*s912))/2.);
				buffer__[12*i+10] += (s152*s286*s940)/4. + (s152*s287*s288*s940)/12. + (s152*s325*s326*s940)/12. + (s152*s327*s328*s940)/12. + s323*(s329 + (s152*s324*s940)/12.) + s407*((s152*s416*s940)/4. + (s151*(s205*(s1012 - s205*s210*s225*s943 - s197*s223*s225*s943) + s197*(s803 + s804 - s197*s215*s225*s943 - s205*s223*s225*s943)))/2.) + s330*((s152*s343*s940)/4. + (s151*(s199*s946 + s188*s949))/2.) + s344*((s152*s359*s940)/4. + (s151*(s201*s946 + s191*s949))/2.) + s360*((s152*s363*s940)/4. + (s151*(s248 + s249 + s203*s946 + s194*s949))/2.) + s364*((s152*s367*s940)/4. + (s151*(s205*s946 + s197*s949))/2.) + s368*((s152*s372*s940)/4. + (s151*(s201*s952 + s191*s955))/2.) + s373*((s152*s384*s940)/4. + (s151*(s264 + s265 + s203*s952 + s194*s955))/2.) + s385*((s152*s388*s940)/4. + (s151*(s205*s952 + s197*s955))/2.) + s389*((s152*s392*s940)/4. + (s151*(s280 + s281 + s203*s959 + s194*s962))/2.) + s393*((s152*s406*s940)/4. + (s151*(s205*s959 + s197*s962))/2.);
				buffer__[12*i+11] += (s152*s286*s990)/4. + (s152*s287*s288*s990)/12. + (s152*s323*s324*s990)/12. + (s152*s327*s328*s990)/12. + s325*(s329 + (s152*s326*s990)/12.) + s368*((s151*(s1011*s191 + s1008*s201))/2. + (s152*s372*s990)/4.) + s373*((s151*(s1011*s194 + s1008*s203))/2. + (s152*s384*s990)/4.) + s385*((s151*(s1011*s197 + s1008*s205 + s264 + s265))/2. + (s152*s388*s990)/4.) + s389*((s151*(s1020*s194 + s1017*s203))/2. + (s152*s392*s990)/4.) + s393*((s151*(s1020*s197 + s1017*s205 + s280 + s281))/2. + (s152*s406*s990)/4.) + s407*((s152*s416*s990)/4. + (s151*(s408 + s409 + s205*(-(s198*s218) + s843 - s205*s210*s225*s993 - s197*s223*s225*s993) + s197*(s1021 + s606 - s197*s215*s225*s993 - s205*s223*s225*s993)))/2.) + s330*((s152*s343*s990)/4. + (s151*(s1002*s188 + s199*s999))/2.) + s344*((s152*s359*s990)/4. + (s151*(s1002*s191 + s201*s999))/2.) + s360*((s152*s363*s990)/4. + (s151*(s1002*s194 + s203*s999))/2.) + s364*((s152*s367*s990)/4. + (s151*(s1002*s197 + s248 + s249 + s205*s999))/2.);
			}
		}
		else
		{
			#pragma omp parallel for num_threads( ThreadCount() )
			for( Int i = 0; i < simplices.Dimension(0); ++i )
			{
				const Real s0 = V_coords__[4*simplices__[3*i+1]+1];
				const Real s1 = s0*s0;
				const Real s2 = V_coords__[4*simplices__[3*i+1]+0];
				const Real s3 = V_coords__[4*simplices__[3*i+0]+0];
				const Real s4 = V_coords__[4*simplices__[3*i+1]+2];
				const Real s5 = s4*s4;
				const Real s6 = V_coords__[4*simplices__[3*i+1]+3];
				const Real s7 = s6*s6;
				const Real s8 = V_coords__[4*simplices__[3*i+0]+1];
				const Real s9 = V_coords__[4*simplices__[3*i+2]+0];
				const Real s10 = V_coords__[4*simplices__[3*i+0]+2];
				const Real s11 = V_coords__[4*simplices__[3*i+0]+3];
				const Real s12 = V_coords__[4*simplices__[3*i+2]+1];
				const Real s13 = s12*s12;
				const Real s14 = V_coords__[4*simplices__[3*i+2]+2];
				const Real s15 = s14*s14;
				const Real s16 = V_coords__[4*simplices__[3*i+2]+3];
				const Real s17 = s16*s16;
				const Real s18 = s8*s8;
				const Real s19 = s2*s2;
				const Real s20 = s18*s19;
				const Real s21 = s10*s10;
				const Real s22 = s19*s21;
				const Real s23 = s11*s11;
				const Real s24 = s19*s23;
				const Real s25 = -2*s0*s2*s3*s8;
				const Real s26 = s3*s3;
				const Real s27 = s1*s26;
				const Real s28 = s1*s21;
				const Real s29 = s1*s23;
				const Real s30 = -2*s10*s2*s3*s4;
				const Real s31 = -2*s0*s10*s4*s8;
				const Real s32 = s26*s5;
				const Real s33 = s18*s5;
				const Real s34 = s23*s5;
				const Real s35 = -2*s11*s2*s3*s6;
				const Real s36 = -2*s0*s11*s6*s8;
				const Real s37 = -2*s10*s11*s4*s6;
				const Real s38 = s26*s7;
				const Real s39 = s18*s7;
				const Real s40 = s21*s7;
				const Real s41 = -2*s18*s2*s9;
				const Real s42 = -2*s2*s21*s9;
				const Real s43 = -2*s2*s23*s9;
				const Real s44 = 2*s0*s3*s8*s9;
				const Real s45 = 2*s0*s2*s8*s9;
				const Real s46 = -2*s1*s3*s9;
				const Real s47 = 2*s10*s3*s4*s9;
				const Real s48 = 2*s10*s2*s4*s9;
				const Real s49 = -2*s3*s5*s9;
				const Real s50 = 2*s11*s3*s6*s9;
				const Real s51 = 2*s11*s2*s6*s9;
				const Real s52 = -2*s3*s7*s9;
				const Real s53 = s9*s9;
				const Real s54 = s18*s53;
				const Real s55 = s21*s53;
				const Real s56 = s23*s53;
				const Real s57 = -2*s0*s53*s8;
				const Real s58 = s1*s53;
				const Real s59 = -2*s10*s4*s53;
				const Real s60 = s5*s53;
				const Real s61 = -2*s11*s53*s6;
				const Real s62 = s53*s7;
				const Real s63 = 2*s12*s2*s3*s8;
				const Real s64 = -2*s12*s19*s8;
				const Real s65 = -2*s0*s12*s26;
				const Real s66 = -2*s0*s12*s21;
				const Real s67 = -2*s0*s12*s23;
				const Real s68 = 2*s0*s12*s2*s3;
				const Real s69 = 2*s10*s12*s4*s8;
				const Real s70 = 2*s0*s10*s12*s4;
				const Real s71 = -2*s12*s5*s8;
				const Real s72 = 2*s11*s12*s6*s8;
				const Real s73 = 2*s0*s11*s12*s6;
				const Real s74 = -2*s12*s7*s8;
				const Real s75 = -2*s12*s3*s8*s9;
				const Real s76 = 2*s12*s2*s8*s9;
				const Real s77 = 2*s0*s12*s3*s9;
				const Real s78 = -2*s0*s12*s2*s9;
				const Real s79 = s13*s26;
				const Real s80 = s13*s21;
				const Real s81 = s13*s23;
				const Real s82 = -2*s13*s2*s3;
				const Real s83 = s13*s19;
				const Real s84 = -2*s10*s13*s4;
				const Real s85 = s13*s5;
				const Real s86 = -2*s11*s13*s6;
				const Real s87 = s13*s7;
				const Real s88 = 2*s10*s14*s2*s3;
				const Real s89 = -2*s10*s14*s19;
				const Real s90 = 2*s0*s10*s14*s8;
				const Real s91 = -2*s1*s10*s14;
				const Real s92 = -2*s14*s26*s4;
				const Real s93 = -2*s14*s18*s4;
				const Real s94 = -2*s14*s23*s4;
				const Real s95 = 2*s14*s2*s3*s4;
				const Real s96 = 2*s0*s14*s4*s8;
				const Real s97 = 2*s10*s11*s14*s6;
				const Real s98 = 2*s11*s14*s4*s6;
				const Real s99 = -2*s10*s14*s7;
				const Real s100 = -2*s10*s14*s3*s9;
				const Real s101 = 2*s10*s14*s2*s9;
				const Real s102 = 2*s14*s3*s4*s9;
				const Real s103 = -2*s14*s2*s4*s9;
				const Real s104 = -2*s10*s12*s14*s8;
				const Real s105 = 2*s0*s10*s12*s14;
				const Real s106 = 2*s12*s14*s4*s8;
				const Real s107 = -2*s0*s12*s14*s4;
				const Real s108 = s15*s26;
				const Real s109 = s15*s18;
				const Real s110 = s15*s23;
				const Real s111 = -2*s15*s2*s3;
				const Real s112 = s15*s19;
				const Real s113 = -2*s0*s15*s8;
				const Real s114 = s1*s15;
				const Real s115 = -2*s11*s15*s6;
				const Real s116 = s15*s7;
				const Real s117 = 2*s11*s16*s2*s3;
				const Real s118 = -2*s11*s16*s19;
				const Real s119 = 2*s0*s11*s16*s8;
				const Real s120 = -2*s1*s11*s16;
				const Real s121 = 2*s10*s11*s16*s4;
				const Real s122 = -2*s11*s16*s5;
				const Real s123 = -2*s16*s26*s6;
				const Real s124 = -2*s16*s18*s6;
				const Real s125 = -2*s16*s21*s6;
				const Real s126 = 2*s16*s2*s3*s6;
				const Real s127 = 2*s0*s16*s6*s8;
				const Real s128 = 2*s10*s16*s4*s6;
				const Real s129 = -2*s11*s16*s3*s9;
				const Real s130 = 2*s11*s16*s2*s9;
				const Real s131 = 2*s16*s3*s6*s9;
				const Real s132 = -2*s16*s2*s6*s9;
				const Real s133 = -2*s11*s12*s16*s8;
				const Real s134 = 2*s0*s11*s12*s16;
				const Real s135 = 2*s12*s16*s6*s8;
				const Real s136 = -2*s0*s12*s16*s6;
				const Real s137 = -2*s10*s11*s14*s16;
				const Real s138 = 2*s11*s14*s16*s4;
				const Real s139 = 2*s10*s14*s16*s6;
				const Real s140 = -2*s14*s16*s4*s6;
				const Real s141 = s17*s26;
				const Real s142 = s17*s18;
				const Real s143 = s17*s21;
				const Real s144 = -2*s17*s2*s3;
				const Real s145 = s17*s19;
				const Real s146 = -2*s0*s17*s8;
				const Real s147 = s1*s17;
				const Real s148 = -2*s10*s17*s4;
				const Real s149 = s17*s5;
				const Real s150 = s100 + s101 + s102 + s103 + s104 + s105 + s106 + s107 + s108 + s109 + s110 + s111 + s112 + s113 + s114 + s115 + s116 + s117 + s118 + s119 + s120 + s121 + s122 + s123 + s124 + s125 + s126 + s127 + s128 + s129 + s130 + s131 + s132 + s133 + s134 + s135 + s136 + s137 + s138 + s139 + s140 + s141 + s142 + s143 + s144 + s145 + s146 + s147 + s148 + s149 + s20 + s22 + s24 + s25 + s27 + s28 + s29 + s30 + s31 + s32 + s33 + s34 + s35 + s36 + s37 + s38 + s39 + s40 + s41 + s42 + s43 + s44 + s45 + s46 + s47 + s48 + s49 + s50 + s51 + s52 + s54 + s55 + s56 + s57 + s58 + s59 + s60 + s61 + s62 + s63 + s64 + s65 + s66 + s67 + s68 + s69 + s70 + s71 + s72 + s73 + s74 + s75 + s76 + s77 + s78 + s79 + s80 + s81 + s82 + s83 + s84 + s85 + s86 + s87 + s88 + s89 + s90 + s91 + s92 + s93 + s94 + s95 + s96 + s97 + s98 + s99;
				const Real s151 = sqrt(s150);
				const Real s152 = 1/s151;
				const Real s153 = -2*s0*s2*s8;
				const Real s154 = 2*s1*s3;
				const Real s155 = -2*s10*s2*s4;
				const Real s156 = 2*s3*s5;
				const Real s157 = -2*s11*s2*s6;
				const Real s158 = 2*s3*s7;
				const Real s159 = 2*s0*s8*s9;
				const Real s160 = -2*s1*s9;
				const Real s161 = 2*s10*s4*s9;
				const Real s162 = -2*s5*s9;
				const Real s163 = 2*s11*s6*s9;
				const Real s164 = -2*s7*s9;
				const Real s165 = 2*s12*s2*s8;
				const Real s166 = -4*s0*s12*s3;
				const Real s167 = 2*s0*s12*s2;
				const Real s168 = -2*s12*s8*s9;
				const Real s169 = 2*s0*s12*s9;
				const Real s170 = 2*s13*s3;
				const Real s171 = -2*s13*s2;
				const Real s172 = 2*s10*s14*s2;
				const Real s173 = -4*s14*s3*s4;
				const Real s174 = 2*s14*s2*s4;
				const Real s175 = -2*s10*s14*s9;
				const Real s176 = 2*s14*s4*s9;
				const Real s177 = 2*s15*s3;
				const Real s178 = -2*s15*s2;
				const Real s179 = 2*s11*s16*s2;
				const Real s180 = -4*s16*s3*s6;
				const Real s181 = 2*s16*s2*s6;
				const Real s182 = -2*s11*s16*s9;
				const Real s183 = 2*s16*s6*s9;
				const Real s184 = 2*s17*s3;
				const Real s185 = -2*s17*s2;
				const Real s186 = s153 + s154 + s155 + s156 + s157 + s158 + s159 + s160 + s161 + s162 + s163 + s164 + s165 + s166 + s167 + s168 + s169 + s170 + s171 + s172 + s173 + s174 + s175 + s176 + s177 + s178 + s179 + s180 + s181 + s182 + s183 + s184 + s185;
				const Real s187 = -s3;
				const Real s188 = s187 + s2;
				const Real s189 = s188*s188;
				const Real s190 = -s8;
				const Real s191 = s0 + s190;
				const Real s192 = s191*s191;
				const Real s193 = -s10;
				const Real s194 = s193 + s4;
				const Real s195 = s194*s194;
				const Real s196 = -s11;
				const Real s197 = s196 + s6;
				const Real s198 = s197*s197;
				const Real s199 = s187 + s9;
				const Real s200 = s188*s199;
				const Real s201 = s12 + s190;
				const Real s202 = s191*s201;
				const Real s203 = s14 + s193;
				const Real s204 = s194*s203;
				const Real s205 = s16 + s196;
				const Real s206 = s197*s205;
				const Real s207 = s200 + s202 + s204 + s206;
				const Real s208 = s207*s207;
				const Real s209 = -s208;
				const Real s210 = s189 + s192 + s195 + s198;
				const Real s211 = s199*s199;
				const Real s212 = s201*s201;
				const Real s213 = s203*s203;
				const Real s214 = s205*s205;
				const Real s215 = s211 + s212 + s213 + s214;
				const Real s216 = s210*s215;
				const Real s217 = s209 + s216;
				const Real s218 = 1/s217;
				const Real s219 = -(s188*s199);
				const Real s220 = -(s191*s201);
				const Real s221 = -(s194*s203);
				const Real s222 = -(s197*s205);
				const Real s223 = s219 + s220 + s221 + s222;
				const Real s224 = s217*s217;
				const Real s225 = 1/s224;
				const Real s226 = -2*s199*s210;
				const Real s227 = 2*s3;
				const Real s228 = -s2;
				const Real s229 = -s9;
				const Real s230 = s227 + s228 + s229;
				const Real s231 = -2*s207*s230;
				const Real s232 = -2*s188*s215;
				const Real s233 = s226 + s231 + s232;
				const Real s234 = -2*s188*s199*s218;
				const Real s235 = -2*s3;
				const Real s236 = s2 + s235 + s9;
				const Real s237 = -(s218*s223);
				const Real s238 = -(s199*s210*s225*s233);
				const Real s239 = -(s188*s223*s225*s233);
				const Real s240 = -(s210*s218);
				const Real s241 = s188*s218*s236;
				const Real s242 = s234 + s237 + s238 + s239 + s240 + s241;
				const Real s243 = -(s199*s223*s225*s233);
				const Real s244 = -(s188*s215*s225*s233);
				const Real s245 = s199*s218*s236;
				const Real s246 = -(s215*s218);
				const Real s247 = s234 + s237 + s243 + s244 + s245 + s246;
				const Real s248 = s199*s210*s218;
				const Real s249 = s188*s218*s223;
				const Real s250 = s248 + s249;
				const Real s251 = s199*s218*s223;
				const Real s252 = s188*s215*s218;
				const Real s253 = s251 + s252;
				const Real s254 = -(s201*s210*s225*s233);
				const Real s255 = -(s191*s223*s225*s233);
				const Real s256 = s191*s218*s236;
				const Real s257 = -2*s188*s201*s218;
				const Real s258 = s254 + s255 + s256 + s257;
				const Real s259 = -(s201*s223*s225*s233);
				const Real s260 = -(s191*s215*s225*s233);
				const Real s261 = -2*s191*s199*s218;
				const Real s262 = s201*s218*s236;
				const Real s263 = s259 + s260 + s261 + s262;
				const Real s264 = s201*s210*s218;
				const Real s265 = s191*s218*s223;
				const Real s266 = s264 + s265;
				const Real s267 = s201*s218*s223;
				const Real s268 = s191*s215*s218;
				const Real s269 = s267 + s268;
				const Real s270 = -(s203*s210*s225*s233);
				const Real s271 = -(s194*s223*s225*s233);
				const Real s272 = s194*s218*s236;
				const Real s273 = -2*s188*s203*s218;
				const Real s274 = s270 + s271 + s272 + s273;
				const Real s275 = -(s203*s223*s225*s233);
				const Real s276 = -(s194*s215*s225*s233);
				const Real s277 = -2*s194*s199*s218;
				const Real s278 = s203*s218*s236;
				const Real s279 = s275 + s276 + s277 + s278;
				const Real s280 = s203*s210*s218;
				const Real s281 = s194*s218*s223;
				const Real s282 = s280 + s281;
				const Real s283 = s203*s218*s223;
				const Real s284 = s194*s215*s218;
				const Real s285 = s283 + s284;
				const Real s286 = P_D_far__[15*i+0];
				const Real s287 = P_D_far__[15*i+1];
				const Real s288 = s2 + s3 + s9;
				const Real s289 = 2*s19*s8;
				const Real s290 = -2*s0*s2*s3;
				const Real s291 = -2*s0*s10*s4;
				const Real s292 = 2*s5*s8;
				const Real s293 = -2*s0*s11*s6;
				const Real s294 = 2*s7*s8;
				const Real s295 = -4*s2*s8*s9;
				const Real s296 = 2*s0*s3*s9;
				const Real s297 = 2*s0*s2*s9;
				const Real s298 = 2*s53*s8;
				const Real s299 = -2*s0*s53;
				const Real s300 = 2*s12*s2*s3;
				const Real s301 = -2*s12*s19;
				const Real s302 = 2*s10*s12*s4;
				const Real s303 = -2*s12*s5;
				const Real s304 = 2*s11*s12*s6;
				const Real s305 = -2*s12*s7;
				const Real s306 = -2*s12*s3*s9;
				const Real s307 = 2*s12*s2*s9;
				const Real s308 = 2*s0*s10*s14;
				const Real s309 = -4*s14*s4*s8;
				const Real s310 = 2*s0*s14*s4;
				const Real s311 = -2*s10*s12*s14;
				const Real s312 = 2*s12*s14*s4;
				const Real s313 = 2*s15*s8;
				const Real s314 = -2*s0*s15;
				const Real s315 = 2*s0*s11*s16;
				const Real s316 = -4*s16*s6*s8;
				const Real s317 = 2*s0*s16*s6;
				const Real s318 = -2*s11*s12*s16;
				const Real s319 = 2*s12*s16*s6;
				const Real s320 = 2*s17*s8;
				const Real s321 = -2*s0*s17;
				const Real s322 = s289 + s290 + s291 + s292 + s293 + s294 + s295 + s296 + s297 + s298 + s299 + s300 + s301 + s302 + s303 + s304 + s305 + s306 + s307 + s308 + s309 + s310 + s311 + s312 + s313 + s314 + s315 + s316 + s317 + s318 + s319 + s320 + s321;
				const Real s323 = P_D_far__[15*i+3];
				const Real s324 = s10 + s14 + s4;
				const Real s325 = P_D_far__[15*i+4];
				const Real s326 = s11 + s16 + s6;
				const Real s327 = P_D_far__[15*i+2];
				const Real s328 = s0 + s12 + s8;
				const Real s329 = s151/6.;
				const Real s330 = P_D_far__[15*i+5];
				const Real s331 = -2*s201*s210;
				const Real s332 = 2*s8;
				const Real s333 = -s0;
				const Real s334 = -s12;
				const Real s335 = s332 + s333 + s334;
				const Real s336 = -2*s207*s335;
				const Real s337 = -2*s191*s215;
				const Real s338 = s331 + s336 + s337;
				const Real s339 = -2*s8;
				const Real s340 = s0 + s12 + s339;
				const Real s341 = s199*s250;
				const Real s342 = s188*s253;
				const Real s343 = s341 + s342;
				const Real s344 = P_D_far__[15*i+6];
				const Real s345 = -(s199*s210*s218);
				const Real s346 = -(s188*s218*s223);
				const Real s347 = -(s199*s218*s223);
				const Real s348 = -(s188*s215*s218);
				const Real s349 = -(s199*s210*s225*s338);
				const Real s350 = -(s188*s223*s225*s338);
				const Real s351 = s188*s218*s340;
				const Real s352 = s261 + s349 + s350 + s351;
				const Real s353 = -(s199*s223*s225*s338);
				const Real s354 = -(s188*s215*s225*s338);
				const Real s355 = s199*s218*s340;
				const Real s356 = s257 + s353 + s354 + s355;
				const Real s357 = s201*s250;
				const Real s358 = s191*s253;
				const Real s359 = s357 + s358;
				const Real s360 = P_D_far__[15*i+7];
				const Real s361 = s203*s250;
				const Real s362 = s194*s253;
				const Real s363 = s361 + s362;
				const Real s364 = P_D_far__[15*i+8];
				const Real s365 = s205*s250;
				const Real s366 = s197*s253;
				const Real s367 = s365 + s366;
				const Real s368 = P_D_far__[15*i+9];
				const Real s369 = -2*s191*s201*s218;
				const Real s370 = s201*s266;
				const Real s371 = s191*s269;
				const Real s372 = s370 + s371;
				const Real s373 = P_D_far__[15*i+10];
				const Real s374 = -(s201*s210*s225*s338);
				const Real s375 = -(s191*s223*s225*s338);
				const Real s376 = s191*s218*s340;
				const Real s377 = s237 + s240 + s369 + s374 + s375 + s376;
				const Real s378 = -(s201*s223*s225*s338);
				const Real s379 = -(s191*s215*s225*s338);
				const Real s380 = s201*s218*s340;
				const Real s381 = s237 + s246 + s369 + s378 + s379 + s380;
				const Real s382 = s203*s266;
				const Real s383 = s194*s269;
				const Real s384 = s382 + s383;
				const Real s385 = P_D_far__[15*i+11];
				const Real s386 = s205*s266;
				const Real s387 = s197*s269;
				const Real s388 = s386 + s387;
				const Real s389 = P_D_far__[15*i+12];
				const Real s390 = s203*s282;
				const Real s391 = s194*s285;
				const Real s392 = s390 + s391;
				const Real s393 = P_D_far__[15*i+13];
				const Real s394 = -(s203*s210*s225*s338);
				const Real s395 = -(s194*s223*s225*s338);
				const Real s396 = s194*s218*s340;
				const Real s397 = -2*s191*s203*s218;
				const Real s398 = s394 + s395 + s396 + s397;
				const Real s399 = -(s203*s223*s225*s338);
				const Real s400 = -(s194*s215*s225*s338);
				const Real s401 = -2*s194*s201*s218;
				const Real s402 = s203*s218*s340;
				const Real s403 = s399 + s400 + s401 + s402;
				const Real s404 = s205*s282;
				const Real s405 = s197*s285;
				const Real s406 = s404 + s405;
				const Real s407 = P_D_far__[15*i+14];
				const Real s408 = s205*s210*s218;
				const Real s409 = s197*s218*s223;
				const Real s410 = s408 + s409;
				const Real s411 = s205*s410;
				const Real s412 = s205*s218*s223;
				const Real s413 = s197*s215*s218;
				const Real s414 = s412 + s413;
				const Real s415 = s197*s414;
				const Real s416 = s411 + s415;
				const Real s417 = 2*s10*s19;
				const Real s418 = 2*s1*s10;
				const Real s419 = -2*s2*s3*s4;
				const Real s420 = -2*s0*s4*s8;
				const Real s421 = -2*s11*s4*s6;
				const Real s422 = 2*s10*s7;
				const Real s423 = -4*s10*s2*s9;
				const Real s424 = 2*s3*s4*s9;
				const Real s425 = 2*s2*s4*s9;
				const Real s426 = 2*s10*s53;
				const Real s427 = -2*s4*s53;
				const Real s428 = -4*s0*s10*s12;
				const Real s429 = 2*s12*s4*s8;
				const Real s430 = 2*s0*s12*s4;
				const Real s431 = 2*s10*s13;
				const Real s432 = -2*s13*s4;
				const Real s433 = 2*s14*s2*s3;
				const Real s434 = -2*s14*s19;
				const Real s435 = 2*s0*s14*s8;
				const Real s436 = -2*s1*s14;
				const Real s437 = 2*s11*s14*s6;
				const Real s438 = -2*s14*s7;
				const Real s439 = -2*s14*s3*s9;
				const Real s440 = 2*s14*s2*s9;
				const Real s441 = -2*s12*s14*s8;
				const Real s442 = 2*s0*s12*s14;
				const Real s443 = 2*s11*s16*s4;
				const Real s444 = -4*s10*s16*s6;
				const Real s445 = 2*s16*s4*s6;
				const Real s446 = -2*s11*s14*s16;
				const Real s447 = 2*s14*s16*s6;
				const Real s448 = 2*s10*s17;
				const Real s449 = -2*s17*s4;
				const Real s450 = s417 + s418 + s419 + s420 + s421 + s422 + s423 + s424 + s425 + s426 + s427 + s428 + s429 + s430 + s431 + s432 + s433 + s434 + s435 + s436 + s437 + s438 + s439 + s440 + s441 + s442 + s443 + s444 + s445 + s446 + s447 + s448 + s449;
				const Real s451 = -2*s203*s210;
				const Real s452 = 2*s10;
				const Real s453 = -s4;
				const Real s454 = -s14;
				const Real s455 = s452 + s453 + s454;
				const Real s456 = -2*s207*s455;
				const Real s457 = -2*s194*s215;
				const Real s458 = s451 + s456 + s457;
				const Real s459 = -2*s10;
				const Real s460 = s14 + s4 + s459;
				const Real s461 = -(s199*s210*s225*s458);
				const Real s462 = -(s188*s223*s225*s458);
				const Real s463 = s188*s218*s460;
				const Real s464 = s277 + s461 + s462 + s463;
				const Real s465 = -(s199*s223*s225*s458);
				const Real s466 = -(s188*s215*s225*s458);
				const Real s467 = s199*s218*s460;
				const Real s468 = s273 + s465 + s466 + s467;
				const Real s469 = -(s201*s210*s218);
				const Real s470 = -(s191*s218*s223);
				const Real s471 = -(s201*s218*s223);
				const Real s472 = -(s191*s215*s218);
				const Real s473 = -(s201*s210*s225*s458);
				const Real s474 = -(s191*s223*s225*s458);
				const Real s475 = s191*s218*s460;
				const Real s476 = s401 + s473 + s474 + s475;
				const Real s477 = -(s201*s223*s225*s458);
				const Real s478 = -(s191*s215*s225*s458);
				const Real s479 = s201*s218*s460;
				const Real s480 = s397 + s477 + s478 + s479;
				const Real s481 = -2*s194*s203*s218;
				const Real s482 = -(s203*s210*s225*s458);
				const Real s483 = -(s194*s223*s225*s458);
				const Real s484 = s194*s218*s460;
				const Real s485 = s237 + s240 + s481 + s482 + s483 + s484;
				const Real s486 = -(s203*s223*s225*s458);
				const Real s487 = -(s194*s215*s225*s458);
				const Real s488 = s203*s218*s460;
				const Real s489 = s237 + s246 + s481 + s486 + s487 + s488;
				const Real s490 = 2*s11*s19;
				const Real s491 = 2*s1*s11;
				const Real s492 = 2*s11*s5;
				const Real s493 = -2*s2*s3*s6;
				const Real s494 = -2*s0*s6*s8;
				const Real s495 = -2*s10*s4*s6;
				const Real s496 = -4*s11*s2*s9;
				const Real s497 = 2*s3*s6*s9;
				const Real s498 = 2*s2*s6*s9;
				const Real s499 = 2*s11*s53;
				const Real s500 = -2*s53*s6;
				const Real s501 = -4*s0*s11*s12;
				const Real s502 = 2*s12*s6*s8;
				const Real s503 = 2*s0*s12*s6;
				const Real s504 = 2*s11*s13;
				const Real s505 = -2*s13*s6;
				const Real s506 = -4*s11*s14*s4;
				const Real s507 = 2*s10*s14*s6;
				const Real s508 = 2*s14*s4*s6;
				const Real s509 = 2*s11*s15;
				const Real s510 = -2*s15*s6;
				const Real s511 = 2*s16*s2*s3;
				const Real s512 = -2*s16*s19;
				const Real s513 = 2*s0*s16*s8;
				const Real s514 = -2*s1*s16;
				const Real s515 = 2*s10*s16*s4;
				const Real s516 = -2*s16*s5;
				const Real s517 = -2*s16*s3*s9;
				const Real s518 = 2*s16*s2*s9;
				const Real s519 = -2*s12*s16*s8;
				const Real s520 = 2*s0*s12*s16;
				const Real s521 = -2*s10*s14*s16;
				const Real s522 = 2*s14*s16*s4;
				const Real s523 = s490 + s491 + s492 + s493 + s494 + s495 + s496 + s497 + s498 + s499 + s500 + s501 + s502 + s503 + s504 + s505 + s506 + s507 + s508 + s509 + s510 + s511 + s512 + s513 + s514 + s515 + s516 + s517 + s518 + s519 + s520 + s521 + s522;
				const Real s524 = -2*s205*s210;
				const Real s525 = 2*s11;
				const Real s526 = -s6;
				const Real s527 = -s16;
				const Real s528 = s525 + s526 + s527;
				const Real s529 = -2*s207*s528;
				const Real s530 = -2*s197*s215;
				const Real s531 = s524 + s529 + s530;
				const Real s532 = -2*s197*s199*s218;
				const Real s533 = -2*s188*s205*s218;
				const Real s534 = -2*s11;
				const Real s535 = s16 + s534 + s6;
				const Real s536 = -(s199*s210*s225*s531);
				const Real s537 = -(s188*s223*s225*s531);
				const Real s538 = s188*s218*s535;
				const Real s539 = s532 + s536 + s537 + s538;
				const Real s540 = -(s199*s223*s225*s531);
				const Real s541 = -(s188*s215*s225*s531);
				const Real s542 = s199*s218*s535;
				const Real s543 = s533 + s540 + s541 + s542;
				const Real s544 = -2*s197*s201*s218;
				const Real s545 = -2*s191*s205*s218;
				const Real s546 = -(s201*s210*s225*s531);
				const Real s547 = -(s191*s223*s225*s531);
				const Real s548 = s191*s218*s535;
				const Real s549 = s544 + s546 + s547 + s548;
				const Real s550 = -(s201*s223*s225*s531);
				const Real s551 = -(s191*s215*s225*s531);
				const Real s552 = s201*s218*s535;
				const Real s553 = s545 + s550 + s551 + s552;
				const Real s554 = -2*s197*s203*s218;
				const Real s555 = -2*s194*s205*s218;
				const Real s556 = -(s203*s210*s218);
				const Real s557 = -(s194*s218*s223);
				const Real s558 = -(s203*s218*s223);
				const Real s559 = -(s194*s215*s218);
				const Real s560 = -(s203*s210*s225*s531);
				const Real s561 = -(s194*s223*s225*s531);
				const Real s562 = s194*s218*s535;
				const Real s563 = s554 + s560 + s561 + s562;
				const Real s564 = -(s203*s223*s225*s531);
				const Real s565 = -(s194*s215*s225*s531);
				const Real s566 = s203*s218*s535;
				const Real s567 = s555 + s564 + s565 + s566;
				const Real s568 = -2*s197*s205*s218;
				const Real s569 = 2*s18*s2;
				const Real s570 = 2*s2*s21;
				const Real s571 = 2*s2*s23;
				const Real s572 = -2*s0*s3*s8;
				const Real s573 = -2*s10*s3*s4;
				const Real s574 = -2*s11*s3*s6;
				const Real s575 = -2*s18*s9;
				const Real s576 = -2*s21*s9;
				const Real s577 = -2*s23*s9;
				const Real s578 = 2*s12*s3*s8;
				const Real s579 = -4*s12*s2*s8;
				const Real s580 = 2*s0*s12*s3;
				const Real s581 = 2*s12*s8*s9;
				const Real s582 = -2*s0*s12*s9;
				const Real s583 = -2*s13*s3;
				const Real s584 = 2*s13*s2;
				const Real s585 = 2*s10*s14*s3;
				const Real s586 = -4*s10*s14*s2;
				const Real s587 = 2*s14*s3*s4;
				const Real s588 = 2*s10*s14*s9;
				const Real s589 = -2*s14*s4*s9;
				const Real s590 = -2*s15*s3;
				const Real s591 = 2*s15*s2;
				const Real s592 = 2*s11*s16*s3;
				const Real s593 = -4*s11*s16*s2;
				const Real s594 = 2*s16*s3*s6;
				const Real s595 = 2*s11*s16*s9;
				const Real s596 = -2*s16*s6*s9;
				const Real s597 = -2*s17*s3;
				const Real s598 = 2*s17*s2;
				const Real s599 = s159 + s161 + s163 + s569 + s570 + s571 + s572 + s573 + s574 + s575 + s576 + s577 + s578 + s579 + s580 + s581 + s582 + s583 + s584 + s585 + s586 + s587 + s588 + s589 + s590 + s591 + s592 + s593 + s594 + s595 + s596 + s597 + s598;
				const Real s600 = -2*s199*s207;
				const Real s601 = 2*s188*s215;
				const Real s602 = s600 + s601;
				const Real s603 = -(s199*s210*s225*s602);
				const Real s604 = -(s188*s223*s225*s602);
				const Real s605 = s188*s199*s218;
				const Real s606 = s218*s223;
				const Real s607 = s603 + s604 + s605 + s606;
				const Real s608 = -(s199*s223*s225*s602);
				const Real s609 = -(s188*s215*s225*s602);
				const Real s610 = -(s211*s218);
				const Real s611 = s215*s218;
				const Real s612 = s608 + s609 + s610 + s611;
				const Real s613 = -(s201*s210*s225*s602);
				const Real s614 = -(s191*s223*s225*s602);
				const Real s615 = -(s191*s199*s218);
				const Real s616 = 2*s188*s201*s218;
				const Real s617 = s613 + s614 + s615 + s616;
				const Real s618 = -(s201*s223*s225*s602);
				const Real s619 = -(s191*s215*s225*s602);
				const Real s620 = -(s199*s201*s218);
				const Real s621 = s618 + s619 + s620;
				const Real s622 = -(s203*s210*s225*s602);
				const Real s623 = -(s194*s223*s225*s602);
				const Real s624 = -(s194*s199*s218);
				const Real s625 = 2*s188*s203*s218;
				const Real s626 = s622 + s623 + s624 + s625;
				const Real s627 = -(s203*s223*s225*s602);
				const Real s628 = -(s194*s215*s225*s602);
				const Real s629 = -(s199*s203*s218);
				const Real s630 = s627 + s628 + s629;
				const Real s631 = -2*s2*s3*s8;
				const Real s632 = 2*s0*s26;
				const Real s633 = 2*s0*s21;
				const Real s634 = 2*s0*s23;
				const Real s635 = -2*s10*s4*s8;
				const Real s636 = -2*s11*s6*s8;
				const Real s637 = 2*s3*s8*s9;
				const Real s638 = 2*s2*s8*s9;
				const Real s639 = -4*s0*s3*s9;
				const Real s640 = -2*s53*s8;
				const Real s641 = 2*s0*s53;
				const Real s642 = -2*s12*s26;
				const Real s643 = -2*s12*s21;
				const Real s644 = -2*s12*s23;
				const Real s645 = 2*s12*s3*s9;
				const Real s646 = -2*s12*s2*s9;
				const Real s647 = 2*s10*s14*s8;
				const Real s648 = -4*s0*s10*s14;
				const Real s649 = 2*s14*s4*s8;
				const Real s650 = 2*s10*s12*s14;
				const Real s651 = -2*s12*s14*s4;
				const Real s652 = -2*s15*s8;
				const Real s653 = 2*s0*s15;
				const Real s654 = 2*s11*s16*s8;
				const Real s655 = -4*s0*s11*s16;
				const Real s656 = 2*s16*s6*s8;
				const Real s657 = 2*s11*s12*s16;
				const Real s658 = -2*s12*s16*s6;
				const Real s659 = -2*s17*s8;
				const Real s660 = 2*s0*s17;
				const Real s661 = s300 + s302 + s304 + s631 + s632 + s633 + s634 + s635 + s636 + s637 + s638 + s639 + s640 + s641 + s642 + s643 + s644 + s645 + s646 + s647 + s648 + s649 + s650 + s651 + s652 + s653 + s654 + s655 + s656 + s657 + s658 + s659 + s660;
				const Real s662 = -2*s201*s207;
				const Real s663 = 2*s191*s215;
				const Real s664 = s662 + s663;
				const Real s665 = -(s199*s210*s225*s664);
				const Real s666 = -(s188*s223*s225*s664);
				const Real s667 = 2*s191*s199*s218;
				const Real s668 = -(s188*s201*s218);
				const Real s669 = s665 + s666 + s667 + s668;
				const Real s670 = -(s199*s223*s225*s664);
				const Real s671 = -(s188*s215*s225*s664);
				const Real s672 = s620 + s670 + s671;
				const Real s673 = -(s201*s210*s225*s664);
				const Real s674 = -(s191*s223*s225*s664);
				const Real s675 = s191*s201*s218;
				const Real s676 = s606 + s673 + s674 + s675;
				const Real s677 = -(s201*s223*s225*s664);
				const Real s678 = -(s191*s215*s225*s664);
				const Real s679 = -(s212*s218);
				const Real s680 = s611 + s677 + s678 + s679;
				const Real s681 = -(s203*s210*s225*s664);
				const Real s682 = -(s194*s223*s225*s664);
				const Real s683 = -(s194*s201*s218);
				const Real s684 = 2*s191*s203*s218;
				const Real s685 = s681 + s682 + s683 + s684;
				const Real s686 = -(s203*s223*s225*s664);
				const Real s687 = -(s194*s215*s225*s664);
				const Real s688 = -(s201*s203*s218);
				const Real s689 = s686 + s687 + s688;
				const Real s690 = -2*s10*s2*s3;
				const Real s691 = -2*s0*s10*s8;
				const Real s692 = 2*s26*s4;
				const Real s693 = 2*s18*s4;
				const Real s694 = 2*s23*s4;
				const Real s695 = -2*s10*s11*s6;
				const Real s696 = 2*s10*s3*s9;
				const Real s697 = 2*s10*s2*s9;
				const Real s698 = -4*s3*s4*s9;
				const Real s699 = -2*s10*s53;
				const Real s700 = 2*s4*s53;
				const Real s701 = 2*s10*s12*s8;
				const Real s702 = 2*s0*s10*s12;
				const Real s703 = -4*s12*s4*s8;
				const Real s704 = -2*s10*s13;
				const Real s705 = 2*s13*s4;
				const Real s706 = -2*s14*s26;
				const Real s707 = -2*s14*s18;
				const Real s708 = -2*s14*s23;
				const Real s709 = 2*s14*s3*s9;
				const Real s710 = -2*s14*s2*s9;
				const Real s711 = 2*s12*s14*s8;
				const Real s712 = -2*s0*s12*s14;
				const Real s713 = 2*s10*s11*s16;
				const Real s714 = -4*s11*s16*s4;
				const Real s715 = 2*s10*s16*s6;
				const Real s716 = 2*s11*s14*s16;
				const Real s717 = -2*s14*s16*s6;
				const Real s718 = -2*s10*s17;
				const Real s719 = 2*s17*s4;
				const Real s720 = s433 + s435 + s437 + s690 + s691 + s692 + s693 + s694 + s695 + s696 + s697 + s698 + s699 + s700 + s701 + s702 + s703 + s704 + s705 + s706 + s707 + s708 + s709 + s710 + s711 + s712 + s713 + s714 + s715 + s716 + s717 + s718 + s719;
				const Real s721 = -2*s203*s207;
				const Real s722 = 2*s194*s215;
				const Real s723 = s721 + s722;
				const Real s724 = -(s199*s210*s225*s723);
				const Real s725 = -(s188*s223*s225*s723);
				const Real s726 = 2*s194*s199*s218;
				const Real s727 = -(s188*s203*s218);
				const Real s728 = s724 + s725 + s726 + s727;
				const Real s729 = -(s199*s223*s225*s723);
				const Real s730 = -(s188*s215*s225*s723);
				const Real s731 = s629 + s729 + s730;
				const Real s732 = -(s201*s210*s225*s723);
				const Real s733 = -(s191*s223*s225*s723);
				const Real s734 = 2*s194*s201*s218;
				const Real s735 = -(s191*s203*s218);
				const Real s736 = s732 + s733 + s734 + s735;
				const Real s737 = -(s201*s223*s225*s723);
				const Real s738 = -(s191*s215*s225*s723);
				const Real s739 = s688 + s737 + s738;
				const Real s740 = -(s203*s210*s225*s723);
				const Real s741 = -(s194*s223*s225*s723);
				const Real s742 = s194*s203*s218;
				const Real s743 = s606 + s740 + s741 + s742;
				const Real s744 = -(s203*s223*s225*s723);
				const Real s745 = -(s194*s215*s225*s723);
				const Real s746 = -(s213*s218);
				const Real s747 = s611 + s744 + s745 + s746;
				const Real s748 = -2*s11*s2*s3;
				const Real s749 = -2*s0*s11*s8;
				const Real s750 = -2*s10*s11*s4;
				const Real s751 = 2*s26*s6;
				const Real s752 = 2*s18*s6;
				const Real s753 = 2*s21*s6;
				const Real s754 = 2*s11*s3*s9;
				const Real s755 = 2*s11*s2*s9;
				const Real s756 = -4*s3*s6*s9;
				const Real s757 = -2*s11*s53;
				const Real s758 = 2*s53*s6;
				const Real s759 = 2*s11*s12*s8;
				const Real s760 = 2*s0*s11*s12;
				const Real s761 = -4*s12*s6*s8;
				const Real s762 = -2*s11*s13;
				const Real s763 = 2*s13*s6;
				const Real s764 = 2*s10*s11*s14;
				const Real s765 = 2*s11*s14*s4;
				const Real s766 = -4*s10*s14*s6;
				const Real s767 = -2*s11*s15;
				const Real s768 = 2*s15*s6;
				const Real s769 = -2*s16*s26;
				const Real s770 = -2*s16*s18;
				const Real s771 = -2*s16*s21;
				const Real s772 = 2*s16*s3*s9;
				const Real s773 = -2*s16*s2*s9;
				const Real s774 = 2*s12*s16*s8;
				const Real s775 = -2*s0*s12*s16;
				const Real s776 = 2*s10*s14*s16;
				const Real s777 = -2*s14*s16*s4;
				const Real s778 = s511 + s513 + s515 + s748 + s749 + s750 + s751 + s752 + s753 + s754 + s755 + s756 + s757 + s758 + s759 + s760 + s761 + s762 + s763 + s764 + s765 + s766 + s767 + s768 + s769 + s770 + s771 + s772 + s773 + s774 + s775 + s776 + s777;
				const Real s779 = -2*s205*s207;
				const Real s780 = 2*s197*s215;
				const Real s781 = s779 + s780;
				const Real s782 = -(s199*s205*s218);
				const Real s783 = -(s199*s210*s225*s781);
				const Real s784 = -(s188*s223*s225*s781);
				const Real s785 = 2*s197*s199*s218;
				const Real s786 = -(s188*s205*s218);
				const Real s787 = s783 + s784 + s785 + s786;
				const Real s788 = -(s199*s223*s225*s781);
				const Real s789 = -(s188*s215*s225*s781);
				const Real s790 = s782 + s788 + s789;
				const Real s791 = -(s201*s205*s218);
				const Real s792 = -(s201*s210*s225*s781);
				const Real s793 = -(s191*s223*s225*s781);
				const Real s794 = 2*s197*s201*s218;
				const Real s795 = -(s191*s205*s218);
				const Real s796 = s792 + s793 + s794 + s795;
				const Real s797 = -(s201*s223*s225*s781);
				const Real s798 = -(s191*s215*s225*s781);
				const Real s799 = s791 + s797 + s798;
				const Real s800 = -(s203*s205*s218);
				const Real s801 = -(s203*s210*s225*s781);
				const Real s802 = -(s194*s223*s225*s781);
				const Real s803 = 2*s197*s203*s218;
				const Real s804 = -(s194*s205*s218);
				const Real s805 = s801 + s802 + s803 + s804;
				const Real s806 = -(s203*s223*s225*s781);
				const Real s807 = -(s194*s215*s225*s781);
				const Real s808 = s800 + s806 + s807;
				const Real s809 = -2*s18*s2;
				const Real s810 = -2*s2*s21;
				const Real s811 = -2*s2*s23;
				const Real s812 = 2*s0*s3*s8;
				const Real s813 = 2*s0*s2*s8;
				const Real s814 = -2*s1*s3;
				const Real s815 = 2*s10*s3*s4;
				const Real s816 = 2*s10*s2*s4;
				const Real s817 = -2*s3*s5;
				const Real s818 = 2*s11*s3*s6;
				const Real s819 = 2*s11*s2*s6;
				const Real s820 = -2*s3*s7;
				const Real s821 = 2*s18*s9;
				const Real s822 = 2*s21*s9;
				const Real s823 = 2*s23*s9;
				const Real s824 = -4*s0*s8*s9;
				const Real s825 = 2*s1*s9;
				const Real s826 = -4*s10*s4*s9;
				const Real s827 = 2*s5*s9;
				const Real s828 = -4*s11*s6*s9;
				const Real s829 = 2*s7*s9;
				const Real s830 = -2*s12*s3*s8;
				const Real s831 = -2*s0*s12*s2;
				const Real s832 = -2*s10*s14*s3;
				const Real s833 = -2*s14*s2*s4;
				const Real s834 = -2*s11*s16*s3;
				const Real s835 = -2*s16*s2*s6;
				const Real s836 = s165 + s172 + s179 + s580 + s587 + s594 + s809 + s810 + s811 + s812 + s813 + s814 + s815 + s816 + s817 + s818 + s819 + s820 + s821 + s822 + s823 + s824 + s825 + s826 + s827 + s828 + s829 + s830 + s831 + s832 + s833 + s834 + s835;
				const Real s837 = 2*s199*s210;
				const Real s838 = -2*s188*s207;
				const Real s839 = s837 + s838;
				const Real s840 = -(s199*s210*s225*s839);
				const Real s841 = -(s188*s223*s225*s839);
				const Real s842 = -(s189*s218);
				const Real s843 = s210*s218;
				const Real s844 = s840 + s841 + s842 + s843;
				const Real s845 = -(s199*s223*s225*s839);
				const Real s846 = -(s188*s215*s225*s839);
				const Real s847 = s605 + s606 + s845 + s846;
				const Real s848 = -(s201*s210*s225*s839);
				const Real s849 = -(s191*s223*s225*s839);
				const Real s850 = -(s188*s191*s218);
				const Real s851 = s848 + s849 + s850;
				const Real s852 = -(s201*s223*s225*s839);
				const Real s853 = -(s191*s215*s225*s839);
				const Real s854 = s667 + s668 + s852 + s853;
				const Real s855 = -(s203*s210*s225*s839);
				const Real s856 = -(s194*s223*s225*s839);
				const Real s857 = -(s188*s194*s218);
				const Real s858 = s855 + s856 + s857;
				const Real s859 = -(s203*s223*s225*s839);
				const Real s860 = -(s194*s215*s225*s839);
				const Real s861 = s726 + s727 + s859 + s860;
				const Real s862 = 2*s2*s3*s8;
				const Real s863 = -2*s19*s8;
				const Real s864 = -2*s0*s26;
				const Real s865 = -2*s0*s21;
				const Real s866 = -2*s0*s23;
				const Real s867 = 2*s0*s2*s3;
				const Real s868 = 2*s10*s4*s8;
				const Real s869 = 2*s0*s10*s4;
				const Real s870 = -2*s5*s8;
				const Real s871 = 2*s11*s6*s8;
				const Real s872 = 2*s0*s11*s6;
				const Real s873 = -2*s7*s8;
				const Real s874 = -2*s3*s8*s9;
				const Real s875 = -2*s0*s2*s9;
				const Real s876 = 2*s12*s26;
				const Real s877 = 2*s12*s21;
				const Real s878 = 2*s12*s23;
				const Real s879 = -4*s12*s2*s3;
				const Real s880 = 2*s12*s19;
				const Real s881 = -4*s10*s12*s4;
				const Real s882 = 2*s12*s5;
				const Real s883 = -4*s11*s12*s6;
				const Real s884 = 2*s12*s7;
				const Real s885 = -2*s10*s14*s8;
				const Real s886 = -2*s0*s14*s4;
				const Real s887 = -2*s11*s16*s8;
				const Real s888 = -2*s0*s16*s6;
				const Real s889 = s296 + s308 + s315 + s638 + s649 + s656 + s862 + s863 + s864 + s865 + s866 + s867 + s868 + s869 + s870 + s871 + s872 + s873 + s874 + s875 + s876 + s877 + s878 + s879 + s880 + s881 + s882 + s883 + s884 + s885 + s886 + s887 + s888;
				const Real s890 = 2*s201*s210;
				const Real s891 = -2*s191*s207;
				const Real s892 = s890 + s891;
				const Real s893 = -(s199*s210*s225*s892);
				const Real s894 = -(s188*s223*s225*s892);
				const Real s895 = s850 + s893 + s894;
				const Real s896 = -(s199*s223*s225*s892);
				const Real s897 = -(s188*s215*s225*s892);
				const Real s898 = s615 + s616 + s896 + s897;
				const Real s899 = -(s201*s210*s225*s892);
				const Real s900 = -(s191*s223*s225*s892);
				const Real s901 = -(s192*s218);
				const Real s902 = s843 + s899 + s900 + s901;
				const Real s903 = -(s201*s223*s225*s892);
				const Real s904 = -(s191*s215*s225*s892);
				const Real s905 = s606 + s675 + s903 + s904;
				const Real s906 = -(s203*s210*s225*s892);
				const Real s907 = -(s194*s223*s225*s892);
				const Real s908 = -(s191*s194*s218);
				const Real s909 = s906 + s907 + s908;
				const Real s910 = -(s203*s223*s225*s892);
				const Real s911 = -(s194*s215*s225*s892);
				const Real s912 = s734 + s735 + s910 + s911;
				const Real s913 = 2*s10*s2*s3;
				const Real s914 = -2*s10*s19;
				const Real s915 = 2*s0*s10*s8;
				const Real s916 = -2*s1*s10;
				const Real s917 = -2*s26*s4;
				const Real s918 = -2*s18*s4;
				const Real s919 = -2*s23*s4;
				const Real s920 = 2*s2*s3*s4;
				const Real s921 = 2*s0*s4*s8;
				const Real s922 = 2*s10*s11*s6;
				const Real s923 = 2*s11*s4*s6;
				const Real s924 = -2*s10*s7;
				const Real s925 = -2*s10*s3*s9;
				const Real s926 = -2*s2*s4*s9;
				const Real s927 = -2*s10*s12*s8;
				const Real s928 = -2*s0*s12*s4;
				const Real s929 = 2*s14*s26;
				const Real s930 = 2*s14*s18;
				const Real s931 = 2*s14*s23;
				const Real s932 = -4*s14*s2*s3;
				const Real s933 = 2*s14*s19;
				const Real s934 = -4*s0*s14*s8;
				const Real s935 = 2*s1*s14;
				const Real s936 = -4*s11*s14*s6;
				const Real s937 = 2*s14*s7;
				const Real s938 = -2*s10*s11*s16;
				const Real s939 = -2*s16*s4*s6;
				const Real s940 = s424 + s429 + s443 + s697 + s702 + s715 + s913 + s914 + s915 + s916 + s917 + s918 + s919 + s920 + s921 + s922 + s923 + s924 + s925 + s926 + s927 + s928 + s929 + s930 + s931 + s932 + s933 + s934 + s935 + s936 + s937 + s938 + s939;
				const Real s941 = 2*s203*s210;
				const Real s942 = -2*s194*s207;
				const Real s943 = s941 + s942;
				const Real s944 = -(s199*s210*s225*s943);
				const Real s945 = -(s188*s223*s225*s943);
				const Real s946 = s857 + s944 + s945;
				const Real s947 = -(s199*s223*s225*s943);
				const Real s948 = -(s188*s215*s225*s943);
				const Real s949 = s624 + s625 + s947 + s948;
				const Real s950 = -(s201*s210*s225*s943);
				const Real s951 = -(s191*s223*s225*s943);
				const Real s952 = s908 + s950 + s951;
				const Real s953 = -(s201*s223*s225*s943);
				const Real s954 = -(s191*s215*s225*s943);
				const Real s955 = s683 + s684 + s953 + s954;
				const Real s956 = -(s203*s210*s225*s943);
				const Real s957 = -(s194*s223*s225*s943);
				const Real s958 = -(s195*s218);
				const Real s959 = s843 + s956 + s957 + s958;
				const Real s960 = -(s203*s223*s225*s943);
				const Real s961 = -(s194*s215*s225*s943);
				const Real s962 = s606 + s742 + s960 + s961;
				const Real s963 = 2*s11*s2*s3;
				const Real s964 = -2*s11*s19;
				const Real s965 = 2*s0*s11*s8;
				const Real s966 = -2*s1*s11;
				const Real s967 = 2*s10*s11*s4;
				const Real s968 = -2*s11*s5;
				const Real s969 = -2*s26*s6;
				const Real s970 = -2*s18*s6;
				const Real s971 = -2*s21*s6;
				const Real s972 = 2*s2*s3*s6;
				const Real s973 = 2*s0*s6*s8;
				const Real s974 = 2*s10*s4*s6;
				const Real s975 = -2*s11*s3*s9;
				const Real s976 = -2*s2*s6*s9;
				const Real s977 = -2*s11*s12*s8;
				const Real s978 = -2*s0*s12*s6;
				const Real s979 = -2*s10*s11*s14;
				const Real s980 = -2*s14*s4*s6;
				const Real s981 = 2*s16*s26;
				const Real s982 = 2*s16*s18;
				const Real s983 = 2*s16*s21;
				const Real s984 = -4*s16*s2*s3;
				const Real s985 = 2*s16*s19;
				const Real s986 = -4*s0*s16*s8;
				const Real s987 = 2*s1*s16;
				const Real s988 = -4*s10*s16*s4;
				const Real s989 = 2*s16*s5;
				const Real s990 = s497 + s502 + s507 + s755 + s760 + s765 + s963 + s964 + s965 + s966 + s967 + s968 + s969 + s970 + s971 + s972 + s973 + s974 + s975 + s976 + s977 + s978 + s979 + s980 + s981 + s982 + s983 + s984 + s985 + s986 + s987 + s988 + s989;
				const Real s991 = 2*s205*s210;
				const Real s992 = -2*s197*s207;
				const Real s993 = s991 + s992;
				const Real s994 = -(s188*s197*s218);
				const Real s995 = -(s197*s199*s218);
				const Real s996 = 2*s188*s205*s218;
				const Real s997 = -(s199*s210*s225*s993);
				const Real s998 = -(s188*s223*s225*s993);
				const Real s999 = s994 + s997 + s998;
				const Real s1000 = -(s199*s223*s225*s993);
				const Real s1001 = -(s188*s215*s225*s993);
				const Real s1002 = s1000 + s1001 + s995 + s996;
				const Real s1003 = -(s191*s197*s218);
				const Real s1004 = -(s197*s201*s218);
				const Real s1005 = 2*s191*s205*s218;
				const Real s1006 = -(s201*s210*s225*s993);
				const Real s1007 = -(s191*s223*s225*s993);
				const Real s1008 = s1003 + s1006 + s1007;
				const Real s1009 = -(s201*s223*s225*s993);
				const Real s1010 = -(s191*s215*s225*s993);
				const Real s1011 = s1004 + s1005 + s1009 + s1010;
				const Real s1012 = -(s194*s197*s218);
				const Real s1013 = -(s197*s203*s218);
				const Real s1014 = 2*s194*s205*s218;
				const Real s1015 = -(s203*s210*s225*s993);
				const Real s1016 = -(s194*s223*s225*s993);
				const Real s1017 = s1012 + s1015 + s1016;
				const Real s1018 = -(s203*s223*s225*s993);
				const Real s1019 = -(s194*s215*s225*s993);
				const Real s1020 = s1013 + s1014 + s1018 + s1019;
				const Real s1021 = s197*s205*s218;
				buffer__[12*i+0] = (s152*s186*s286)/4. + (s152*s186*s323*s324)/12. + (s152*s186*s325*s326)/12. + (s152*s186*s327*s328)/12. + s287*((s152*s186*s288)/12. + s329) + s330*((s152*s186*s343)/4. + (s151*(s199*s242 + s188*s247 + s345 + s346 + s347 + s348))/2.) + s344*((s151*(s201*s242 + s191*s247))/2. + (s152*s186*s359)/4.) + s360*((s151*(s203*s242 + s194*s247))/2. + (s152*s186*s363)/4.) + s364*((s151*(s205*s242 + s197*s247))/2. + (s152*s186*s367)/4.) + s368*((s151*(s201*s258 + s191*s263))/2. + (s152*s186*s372)/4.) + s373*((s151*(s203*s258 + s194*s263))/2. + (s152*s186*s384)/4.) + s385*((s151*(s205*s258 + s197*s263))/2. + (s152*s186*s388)/4.) + s389*((s151*(s203*s274 + s194*s279))/2. + (s152*s186*s392)/4.) + s393*((s151*(s205*s274 + s197*s279))/2. + (s152*s186*s406)/4.) + s407*((s152*s186*s416)/4. + (s151*(s197*(-(s197*s215*s225*s233) - s205*s223*s225*s233 + s205*s218*s236 + s532) + s205*(-(s205*s210*s225*s233) - s197*s223*s225*s233 + s197*s218*s236 + s533)))/2.);
				buffer__[12*i+1] = (s152*s286*s322)/4. + (s152*s287*s288*s322)/12. + (s152*s322*s323*s324)/12. + (s152*s322*s325*s326)/12. + s327*((s152*s322*s328)/12. + s329) + s330*((s152*s322*s343)/4. + (s151*(s199*s352 + s188*s356))/2.) + s344*((s151*(s345 + s346 + s347 + s348 + s201*s352 + s191*s356))/2. + (s152*s322*s359)/4.) + s360*((s151*(s203*s352 + s194*s356))/2. + (s152*s322*s363)/4.) + s364*((s151*(s205*s352 + s197*s356))/2. + (s152*s322*s367)/4.) + s373*((s151*(s203*s377 + s194*s381))/2. + (s152*s322*s384)/4.) + s385*((s151*(s205*s377 + s197*s381))/2. + (s152*s322*s388)/4.) + s389*((s152*s322*s392)/4. + (s151*(s203*s398 + s194*s403))/2.) + s393*((s151*(s205*s398 + s197*s403))/2. + (s152*s322*s406)/4.) + s368*((s152*s322*s372)/4. + (s151*(s201*s377 + s191*s381 + s469 + s470 + s471 + s472))/2.) + s407*((s152*s322*s416)/4. + (s151*(s197*(-(s197*s215*s225*s338) - s205*s223*s225*s338 + s205*s218*s340 + s544) + s205*(-(s205*s210*s225*s338) - s197*s223*s225*s338 + s197*s218*s340 + s545)))/2.);
				buffer__[12*i+2] = (s152*s286*s450)/4. + (s152*s287*s288*s450)/12. + (s152*s325*s326*s450)/12. + (s152*s327*s328*s450)/12. + s323*(s329 + (s152*s324*s450)/12.) + s330*((s152*s343*s450)/4. + (s151*(s199*s464 + s188*s468))/2.) + s344*((s152*s359*s450)/4. + (s151*(s201*s464 + s191*s468))/2.) + s360*((s152*s363*s450)/4. + (s151*(s345 + s346 + s347 + s348 + s203*s464 + s194*s468))/2.) + s364*((s152*s367*s450)/4. + (s151*(s205*s464 + s197*s468))/2.) + s368*((s152*s372*s450)/4. + (s151*(s201*s476 + s191*s480))/2.) + s373*((s152*s384*s450)/4. + (s151*(s469 + s470 + s471 + s472 + s203*s476 + s194*s480))/2.) + s385*((s152*s388*s450)/4. + (s151*(s205*s476 + s197*s480))/2.) + s393*((s152*s406*s450)/4. + (s151*(s205*s485 + s197*s489))/2.) + s407*((s152*s416*s450)/4. + (s151*(s197*(-(s197*s215*s225*s458) - s205*s223*s225*s458 + s205*s218*s460 + s554) + s205*(-(s205*s210*s225*s458) - s197*s223*s225*s458 + s197*s218*s460 + s555)))/2.) + s389*((s152*s392*s450)/4. + (s151*(s203*s485 + s194*s489 + s556 + s557 + s558 + s559))/2.);
				buffer__[12*i+3] = (s152*s286*s523)/4. + (s152*s287*s288*s523)/12. + (s152*s323*s324*s523)/12. + (s152*s327*s328*s523)/12. + s325*(s329 + (s152*s326*s523)/12.) + s330*((s152*s343*s523)/4. + (s151*(s199*s539 + s188*s543))/2.) + s344*((s152*s359*s523)/4. + (s151*(s201*s539 + s191*s543))/2.) + s360*((s152*s363*s523)/4. + (s151*(s203*s539 + s194*s543))/2.) + s364*((s152*s367*s523)/4. + (s151*(s345 + s346 + s347 + s348 + s205*s539 + s197*s543))/2.) + s368*((s152*s372*s523)/4. + (s151*(s201*s549 + s191*s553))/2.) + s373*((s152*s384*s523)/4. + (s151*(s203*s549 + s194*s553))/2.) + s385*((s152*s388*s523)/4. + (s151*(s469 + s470 + s471 + s472 + s205*s549 + s197*s553))/2.) + s389*((s152*s392*s523)/4. + (s151*(s203*s563 + s194*s567))/2.) + s393*((s152*s406*s523)/4. + (s151*(s556 + s557 + s558 + s559 + s205*s563 + s197*s567))/2.) + s407*((s152*s416*s523)/4. + (s151*(-(s205*s210*s218) - s197*s215*s218 - s197*s218*s223 - s205*s218*s223 + s205*(s237 + s240 - s205*s210*s225*s531 - s197*s223*s225*s531 + s197*s218*s535 + s568) + s197*(s237 + s246 - s197*s215*s225*s531 - s205*s223*s225*s531 + s205*s218*s535 + s568)))/2.);
				buffer__[12*i+4] = (s152*s286*s599)/4. + (s152*s323*s324*s599)/12. + (s152*s325*s326*s599)/12. + (s152*s327*s328*s599)/12. + s287*(s329 + (s152*s288*s599)/12.) + s330*((s152*s343*s599)/4. + (s151*(s251 + s252 + s199*s607 + s188*s612))/2.) + s344*((s152*s359*s599)/4. + (s151*(s201*s607 + s191*s612))/2.) + s360*((s152*s363*s599)/4. + (s151*(s203*s607 + s194*s612))/2.) + s364*((s152*s367*s599)/4. + (s151*(s205*s607 + s197*s612))/2.) + s368*((s152*s372*s599)/4. + (s151*(s201*s617 + s191*s621))/2.) + s373*((s152*s384*s599)/4. + (s151*(s203*s617 + s194*s621))/2.) + s385*((s152*s388*s599)/4. + (s151*(s205*s617 + s197*s621))/2.) + s389*((s152*s392*s599)/4. + (s151*(s203*s626 + s194*s630))/2.) + s393*((s152*s406*s599)/4. + (s151*(s205*s626 + s197*s630))/2.) + s407*((s152*s416*s599)/4. + (s151*(s197*(-(s197*s215*s225*s602) - s205*s223*s225*s602 + s782) + s205*(-(s205*s210*s225*s602) - s197*s223*s225*s602 + s995 + s996)))/2.);
				buffer__[12*i+5] = (s152*s286*s661)/4. + (s152*s287*s288*s661)/12. + (s152*s323*s324*s661)/12. + (s152*s325*s326*s661)/12. + s327*(s329 + (s152*s328*s661)/12.) + s330*((s152*s343*s661)/4. + (s151*(s199*s669 + s188*s672))/2.) + s344*((s152*s359*s661)/4. + (s151*(s251 + s252 + s201*s669 + s191*s672))/2.) + s360*((s152*s363*s661)/4. + (s151*(s203*s669 + s194*s672))/2.) + s364*((s152*s367*s661)/4. + (s151*(s205*s669 + s197*s672))/2.) + s368*((s152*s372*s661)/4. + (s151*(s267 + s268 + s201*s676 + s191*s680))/2.) + s373*((s152*s384*s661)/4. + (s151*(s203*s676 + s194*s680))/2.) + s385*((s152*s388*s661)/4. + (s151*(s205*s676 + s197*s680))/2.) + s389*((s152*s392*s661)/4. + (s151*(s203*s685 + s194*s689))/2.) + s393*((s152*s406*s661)/4. + (s151*(s205*s685 + s197*s689))/2.) + s407*((s152*s416*s661)/4. + (s151*(s205*(s1004 + s1005 - s205*s210*s225*s664 - s197*s223*s225*s664) + s197*(-(s197*s215*s225*s664) - s205*s223*s225*s664 + s791)))/2.);
				buffer__[12*i+6] = (s152*s286*s720)/4. + (s152*s287*s288*s720)/12. + (s152*s325*s326*s720)/12. + (s152*s327*s328*s720)/12. + s323*(s329 + (s152*s324*s720)/12.) + s330*((s152*s343*s720)/4. + (s151*(s199*s728 + s188*s731))/2.) + s344*((s152*s359*s720)/4. + (s151*(s201*s728 + s191*s731))/2.) + s360*((s152*s363*s720)/4. + (s151*(s251 + s252 + s203*s728 + s194*s731))/2.) + s364*((s152*s367*s720)/4. + (s151*(s205*s728 + s197*s731))/2.) + s368*((s152*s372*s720)/4. + (s151*(s201*s736 + s191*s739))/2.) + s373*((s152*s384*s720)/4. + (s151*(s267 + s268 + s203*s736 + s194*s739))/2.) + s385*((s152*s388*s720)/4. + (s151*(s205*s736 + s197*s739))/2.) + s389*((s152*s392*s720)/4. + (s151*(s283 + s284 + s203*s743 + s194*s747))/2.) + s393*((s152*s406*s720)/4. + (s151*(s205*s743 + s197*s747))/2.) + s407*((s152*s416*s720)/4. + (s151*(s205*(s1013 + s1014 - s205*s210*s225*s723 - s197*s223*s225*s723) + s197*(-(s197*s215*s225*s723) - s205*s223*s225*s723 + s800)))/2.);
				buffer__[12*i+7] = (s152*s286*s778)/4. + (s152*s287*s288*s778)/12. + (s152*s323*s324*s778)/12. + (s152*s327*s328*s778)/12. + s325*(s329 + (s152*s326*s778)/12.) + s407*((s152*s416*s778)/4. + (s151*(s412 + s413 + s205*(s1021 + s606 - s205*s210*s225*s781 - s197*s223*s225*s781) + s197*(-(s214*s218) + s611 - s197*s215*s225*s781 - s205*s223*s225*s781)))/2.) + s330*((s152*s343*s778)/4. + (s151*(s199*s787 + s188*s790))/2.) + s344*((s152*s359*s778)/4. + (s151*(s201*s787 + s191*s790))/2.) + s360*((s152*s363*s778)/4. + (s151*(s203*s787 + s194*s790))/2.) + s364*((s152*s367*s778)/4. + (s151*(s251 + s252 + s205*s787 + s197*s790))/2.) + s368*((s152*s372*s778)/4. + (s151*(s201*s796 + s191*s799))/2.) + s373*((s152*s384*s778)/4. + (s151*(s203*s796 + s194*s799))/2.) + s385*((s152*s388*s778)/4. + (s151*(s267 + s268 + s205*s796 + s197*s799))/2.) + s389*((s152*s392*s778)/4. + (s151*(s203*s805 + s194*s808))/2.) + s393*((s152*s406*s778)/4. + (s151*(s283 + s284 + s205*s805 + s197*s808))/2.);
				buffer__[12*i+8] = (s152*s286*s836)/4. + (s152*s323*s324*s836)/12. + (s152*s325*s326*s836)/12. + (s152*s327*s328*s836)/12. + s287*(s329 + (s152*s288*s836)/12.) + s330*((s152*s343*s836)/4. + (s151*(s248 + s249 + s199*s844 + s188*s847))/2.) + s344*((s152*s359*s836)/4. + (s151*(s201*s844 + s191*s847))/2.) + s360*((s152*s363*s836)/4. + (s151*(s203*s844 + s194*s847))/2.) + s364*((s152*s367*s836)/4. + (s151*(s205*s844 + s197*s847))/2.) + s368*((s152*s372*s836)/4. + (s151*(s201*s851 + s191*s854))/2.) + s373*((s152*s384*s836)/4. + (s151*(s203*s851 + s194*s854))/2.) + s385*((s152*s388*s836)/4. + (s151*(s205*s851 + s197*s854))/2.) + s389*((s152*s392*s836)/4. + (s151*(s203*s858 + s194*s861))/2.) + s393*((s152*s406*s836)/4. + (s151*(s205*s858 + s197*s861))/2.) + s407*((s152*s416*s836)/4. + (s151*(s197*(s785 + s786 - s197*s215*s225*s839 - s205*s223*s225*s839) + s205*(-(s205*s210*s225*s839) - s197*s223*s225*s839 + s994)))/2.);
				buffer__[12*i+9] = (s152*s286*s889)/4. + (s152*s287*s288*s889)/12. + (s152*s323*s324*s889)/12. + (s152*s325*s326*s889)/12. + s327*(s329 + (s152*s328*s889)/12.) + s407*((s152*s416*s889)/4. + (s151*(s205*(s1003 - s205*s210*s225*s892 - s197*s223*s225*s892) + s197*(s794 + s795 - s197*s215*s225*s892 - s205*s223*s225*s892)))/2.) + s330*((s152*s343*s889)/4. + (s151*(s199*s895 + s188*s898))/2.) + s344*((s152*s359*s889)/4. + (s151*(s248 + s249 + s201*s895 + s191*s898))/2.) + s360*((s152*s363*s889)/4. + (s151*(s203*s895 + s194*s898))/2.) + s364*((s152*s367*s889)/4. + (s151*(s205*s895 + s197*s898))/2.) + s368*((s152*s372*s889)/4. + (s151*(s264 + s265 + s201*s902 + s191*s905))/2.) + s373*((s152*s384*s889)/4. + (s151*(s203*s902 + s194*s905))/2.) + s385*((s152*s388*s889)/4. + (s151*(s205*s902 + s197*s905))/2.) + s389*((s152*s392*s889)/4. + (s151*(s203*s909 + s194*s912))/2.) + s393*((s152*s406*s889)/4. + (s151*(s205*s909 + s197*s912))/2.);
				buffer__[12*i+10] = (s152*s286*s940)/4. + (s152*s287*s288*s940)/12. + (s152*s325*s326*s940)/12. + (s152*s327*s328*s940)/12. + s323*(s329 + (s152*s324*s940)/12.) + s407*((s152*s416*s940)/4. + (s151*(s205*(s1012 - s205*s210*s225*s943 - s197*s223*s225*s943) + s197*(s803 + s804 - s197*s215*s225*s943 - s205*s223*s225*s943)))/2.) + s330*((s152*s343*s940)/4. + (s151*(s199*s946 + s188*s949))/2.) + s344*((s152*s359*s940)/4. + (s151*(s201*s946 + s191*s949))/2.) + s360*((s152*s363*s940)/4. + (s151*(s248 + s249 + s203*s946 + s194*s949))/2.) + s364*((s152*s367*s940)/4. + (s151*(s205*s946 + s197*s949))/2.) + s368*((s152*s372*s940)/4. + (s151*(s201*s952 + s191*s955))/2.) + s373*((s152*s384*s940)/4. + (s151*(s264 + s265 + s203*s952 + s194*s955))/2.) + s385*((s152*s388*s940)/4. + (s151*(s205*s952 + s197*s955))/2.) + s389*((s152*s392*s940)/4. + (s151*(s280 + s281 + s203*s959 + s194*s962))/2.) + s393*((s152*s406*s940)/4. + (s151*(s205*s959 + s197*s962))/2.);
				buffer__[12*i+11] = (s152*s286*s990)/4. + (s152*s287*s288*s990)/12. + (s152*s323*s324*s990)/12. + (s152*s327*s328*s990)/12. + s325*(s329 + (s152*s326*s990)/12.) + s368*((s151*(s1011*s191 + s1008*s201))/2. + (s152*s372*s990)/4.) + s373*((s151*(s1011*s194 + s1008*s203))/2. + (s152*s384*s990)/4.) + s385*((s151*(s1011*s197 + s1008*s205 + s264 + s265))/2. + (s152*s388*s990)/4.) + s389*((s151*(s1020*s194 + s1017*s203))/2. + (s152*s392*s990)/4.) + s393*((s151*(s1020*s197 + s1017*s205 + s280 + s281))/2. + (s152*s406*s990)/4.) + s407*((s152*s416*s990)/4. + (s151*(s408 + s409 + s205*(-(s198*s218) + s843 - s205*s210*s225*s993 - s197*s223*s225*s993) + s197*(s1021 + s606 - s197*s215*s225*s993 - s205*s223*s225*s993)))/2.) + s330*((s152*s343*s990)/4. + (s151*(s1002*s188 + s199*s999))/2.) + s344*((s152*s359*s990)/4. + (s151*(s1002*s191 + s201*s999))/2.) + s360*((s152*s363*s990)/4. + (s151*(s1002*s194 + s203*s999))/2.) + s364*((s152*s367*s990)/4. + (s151*(s1002*s197 + s248 + s249 + s205*s999))/2.);
			}
		}

        ptoc(ClassName()+"::DFarToHulls");
        
    }

	}; // SimplicialMeshDetails<2,4,Real,Int>

} // namespace Repulsor