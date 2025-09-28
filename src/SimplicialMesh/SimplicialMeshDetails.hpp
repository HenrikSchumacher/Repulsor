// Machine generated code. Don't edit this file!

#pragma once

namespace Repulsor
{
	template<int DOM_DIM, int AMB_DIM, typename Real, typename Int, typename LInt>
    struct SimplicialMeshDetails
    {
	public:

		using SparseMatrix_T = Sparse::MatrixCSR<Real,Int,LInt>;

		std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<")+ToString(DOM_DIM)+","+ToString(AMB_DIM)+","+TypeName<Real>+","+TypeName<Int>+","+TypeName<LInt>+">";
        }

		
	void DNearFarToHulls( 
		cref<Tensor2<Real,Int>> V_coords, 
		cref<Tensor2<Int ,Int>> simplices, 
		cref<Tensor2<Real,Int>> P_D_near, 
		cref<Tensor2<Real,Int>> P_D_far, 
		cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
		mref<Tensor3<Real,Int>> buffer
	) const
    {
        TOOLS_PTIMER(timer,ClassName()+"::DNearFarToHulls");

		(void)V_coords;
		(void)simplices;
		(void)P_D_near;
		(void)P_D_far;
		(void)V_charges;

        eprint(ClassName()+"::DNearFarToHulls not implemented. Returning 0.");
		
		buffer.Fill(static_cast<Real>(0));
    }

	}; // SimplicialMeshDetails<DOM_DIM,AMB,Real,Int,LInt>

//------------------------------------------------------------

    
    template<typename Real, typename Int, typename LInt>
    struct SimplicialMeshDetails<0,1,Real,Int,LInt>
    {
	public:

		using SparseMatrix_T = Sparse::MatrixCSR<Real,Int,LInt>;

	private:

		const Int thread_count = 1;

		static constexpr Int  SIZE = 1;
		static constexpr Real nth  = Inv<Real>( 1 );

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(Ramp(thread_count_))
		{}
	
		Int ThreadCount() const
		{
			return thread_count;
		}
	
		std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<0,1,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void DNearFarToHulls( 
		cref<Tensor2<Real,Int>> V_coords, 
		cref<Tensor2<Int ,Int>> simplices, 
		cref<Tensor2<Real,Int>> P_D_near,
		cref<Tensor2<Real,Int>> P_D_far, 
		cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
		mref<Tensor3<Real,Int>> buffer
	) const
    {
        TOOLS_PTIMER(timer,ClassName()+"::DNearFarToHulls");

        if( P_D_near.Dim(1) != 3 )
        {
            eprint("in DNearFarToHulls: P_D_near.Dim(1) != 3. Aborting");
        }
		
		(void)V_coords;

		//cptr<Real> X = V_coords.data();
		cptr<Int>  S = simplices.data();
		cptr<Real> N = P_D_near.data();
		cptr<Real> F = P_D_far.data();
		mptr<Real> B = buffer.data();
        
		ParallelDo(
			[=]( const Int i )
			{
				Real charge = 0;
				for( Int k = 0; k < SIZE; ++k )
				{
					charge += V_charges[S[SIZE * i + k]];
				}
				charge *= nth;
				B[1*i+0] = charge*(F[3*i+1] + N[3*i+1]);
			},
			simplices.Dim(0),
			ThreadCount()
		);
    }

	}; // SimplicialMeshDetails<0,1,Real,Int>

//------------------------------------------------------------

    template<typename Real, typename Int, typename LInt>
    struct SimplicialMeshDetails<0,2,Real,Int,LInt>
    {
	public:

		using SparseMatrix_T = Sparse::MatrixCSR<Real,Int,LInt>;

	private:

		const Int thread_count = 1;

		static constexpr Int  SIZE = 1;
		static constexpr Real nth  = Inv<Real>( 1 );

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(Ramp(thread_count_))
		{}
	
		Int ThreadCount() const
		{
			return thread_count;
		}
	
		std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<0,2,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void DNearFarToHulls( 
		cref<Tensor2<Real,Int>> V_coords, 
		cref<Tensor2<Int ,Int>> simplices, 
		cref<Tensor2<Real,Int>> P_D_near,
		cref<Tensor2<Real,Int>> P_D_far, 
		cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
		mref<Tensor3<Real,Int>> buffer
	) const
    {
        TOOLS_PTIMER(timer,ClassName()+"::DNearFarToHulls");

        if( P_D_near.Dim(1) != 6 )
        {
            eprint("in DNearFarToHulls: P_D_near.Dim(1) != 6. Aborting");
        }
		
		(void)V_coords;

		//cptr<Real> X = V_coords.data();
		cptr<Int>  S = simplices.data();
		cptr<Real> N = P_D_near.data();
		cptr<Real> F = P_D_far.data();
		mptr<Real> B = buffer.data();
        
		ParallelDo(
			[=]( const Int i )
			{
				Real charge = 0;
				for( Int k = 0; k < SIZE; ++k )
				{
					charge += V_charges[S[SIZE * i + k]];
				}
				charge *= nth;
				B[2*i+0] = charge*(F[6*i+1] + N[6*i+1]);
				B[2*i+1] = charge*(F[6*i+2] + N[6*i+2]);
			},
			simplices.Dim(0),
			ThreadCount()
		);
    }

	}; // SimplicialMeshDetails<0,2,Real,Int>

//------------------------------------------------------------

    template<typename Real, typename Int, typename LInt>
    struct SimplicialMeshDetails<0,3,Real,Int,LInt>
    {
	public:

		using SparseMatrix_T = Sparse::MatrixCSR<Real,Int,LInt>;

	private:

		const Int thread_count = 1;

		static constexpr Int  SIZE = 1;
		static constexpr Real nth  = Inv<Real>( 1 );

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(Ramp(thread_count_))
		{}
	
		Int ThreadCount() const
		{
			return thread_count;
		}
	
		std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<0,3,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void DNearFarToHulls( 
		cref<Tensor2<Real,Int>> V_coords, 
		cref<Tensor2<Int ,Int>> simplices, 
		cref<Tensor2<Real,Int>> P_D_near,
		cref<Tensor2<Real,Int>> P_D_far, 
		cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
		mref<Tensor3<Real,Int>> buffer
	) const
    {
        TOOLS_PTIMER(timer,ClassName()+"::DNearFarToHulls");

        if( P_D_near.Dim(1) != 10 )
        {
            eprint("in DNearFarToHulls: P_D_near.Dim(1) != 10. Aborting");
        }
		
		(void)V_coords;

		//cptr<Real> X = V_coords.data();
		cptr<Int>  S = simplices.data();
		cptr<Real> N = P_D_near.data();
		cptr<Real> F = P_D_far.data();
		mptr<Real> B = buffer.data();
        
		ParallelDo(
			[=]( const Int i )
			{
				Real charge = 0;
				for( Int k = 0; k < SIZE; ++k )
				{
					charge += V_charges[S[SIZE * i + k]];
				}
				charge *= nth;
				B[3*i+0] = charge*(F[10*i+1] + N[10*i+1]);
				B[3*i+1] = charge*(F[10*i+2] + N[10*i+2]);
				B[3*i+2] = charge*(F[10*i+3] + N[10*i+3]);
			},
			simplices.Dim(0),
			ThreadCount()
		);
    }

	}; // SimplicialMeshDetails<0,3,Real,Int>

//------------------------------------------------------------

    template<typename Real, typename Int, typename LInt>
    struct SimplicialMeshDetails<0,4,Real,Int,LInt>
    {
	public:

		using SparseMatrix_T = Sparse::MatrixCSR<Real,Int,LInt>;

	private:

		const Int thread_count = 1;

		static constexpr Int  SIZE = 1;
		static constexpr Real nth  = Inv<Real>( 1 );

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(Ramp(thread_count_))
		{}
	
		Int ThreadCount() const
		{
			return thread_count;
		}
	
		std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<0,4,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void DNearFarToHulls( 
		cref<Tensor2<Real,Int>> V_coords, 
		cref<Tensor2<Int ,Int>> simplices, 
		cref<Tensor2<Real,Int>> P_D_near,
		cref<Tensor2<Real,Int>> P_D_far, 
		cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
		mref<Tensor3<Real,Int>> buffer
	) const
    {
        TOOLS_PTIMER(timer,ClassName()+"::DNearFarToHulls");

        if( P_D_near.Dim(1) != 15 )
        {
            eprint("in DNearFarToHulls: P_D_near.Dim(1) != 15. Aborting");
        }
		
		(void)V_coords;

		//cptr<Real> X = V_coords.data();
		cptr<Int>  S = simplices.data();
		cptr<Real> N = P_D_near.data();
		cptr<Real> F = P_D_far.data();
		mptr<Real> B = buffer.data();
        
		ParallelDo(
			[=]( const Int i )
			{
				Real charge = 0;
				for( Int k = 0; k < SIZE; ++k )
				{
					charge += V_charges[S[SIZE * i + k]];
				}
				charge *= nth;
				B[4*i+0] = charge*(F[15*i+1] + N[15*i+1]);
				B[4*i+1] = charge*(F[15*i+2] + N[15*i+2]);
				B[4*i+2] = charge*(F[15*i+3] + N[15*i+3]);
				B[4*i+3] = charge*(F[15*i+4] + N[15*i+4]);
			},
			simplices.Dim(0),
			ThreadCount()
		);
    }

	}; // SimplicialMeshDetails<0,4,Real,Int>

//------------------------------------------------------------

    template<typename Real, typename Int, typename LInt>
    struct SimplicialMeshDetails<1,2,Real,Int,LInt>
    {
	public:

		using SparseMatrix_T = Sparse::MatrixCSR<Real,Int,LInt>;

	private:

		const Int thread_count = 1;

		static constexpr Int  SIZE = 2;
		static constexpr Real nth  = Inv<Real>( 2 );

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(Ramp(thread_count_))
		{}
	
		Int ThreadCount() const
		{
			return thread_count;
		}
	
		std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<1,2,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void DNearFarToHulls( 
		cref<Tensor2<Real,Int>> V_coords, 
		cref<Tensor2<Int ,Int>> simplices, 
		cref<Tensor2<Real,Int>> P_D_near, 
		cref<Tensor2<Real,Int>> P_D_far, 
		cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
		mref<Tensor3<Real,Int>> buffer
	) const
    {
        TOOLS_PTIMER(timer,ClassName()+"::DNearFarToHulls");

        if( P_D_near.Dim(1) != 8 )
        {
            eprint("in DNearFarToHulls: P_D_near.Dim(1) != 8. Aborting");
        }

		cptr<Real> X = V_coords.data();
		cptr<Int>  S = simplices.data();
		cptr<Real> N = P_D_near.data();
		cptr<Real> F = P_D_far.data();
		mptr<Real> B = buffer.data();
        
		ParallelDo(
			[=]( const Int i )
			{
				Real charge = 0;
				for( Int k = 0; k < SIZE; ++k )
				{
					charge += V_charges[S[SIZE * i + k]];
				}
				charge *= nth;
				const Real s0 = X[2*S[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = X[2*S[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = X[2*S[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = X[2*S[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = s4 + s9;
				const Real s11 = std::sqrt(s10);
				const Real s12 = 1/s11;
				const Real s13 = s10*s11;
				const Real s14 = 1/s13;
				const Real s15 = -(s14*s4*s8);
				const Real s16 = s12*s8;
				const Real s17 = s15 + s16;
				const Real s18 = s3*s4;
				const Real s19 = s10*s10;
				const Real s20 = 1/s19;
				const Real s21 = 1/s10;
				const Real s22 = -2*s18*s20;
				const Real s23 = 2*s21*s3;
				const Real s24 = s22 + s23;
				const Real s25 = s11*s24;
				const Real s26 = -(s21*s4);
				const Real s27 = 1 + s26;
				const Real s28 = -(s12*s27*s3);
				const Real s29 = s25 + s28;
				const Real s30 = -2*s14*s3*s9;
				const Real s31 = -(s21*s9);
				const Real s32 = 1 + s31;
				const Real s33 = -(s12*s3*s32);
				const Real s34 = s30 + s33;
				const Real s35 = F[6*i+0];
				const Real s36 = N[8*i+0];
				const Real s37 = N[8*i+1];
				const Real s38 = N[8*i+3];
				const Real s39 = F[6*i+1];
				const Real s40 = s0 + s2;
				const Real s41 = N[8*i+4];
				const Real s42 = F[6*i+4];
				const Real s43 = N[8*i+6];
				const Real s44 = -(s14*s3*s9);
				const Real s45 = s12*s3;
				const Real s46 = s44 + s45;
				const Real s47 = F[6*i+2];
				const Real s48 = s5 + s7;
				const Real s49 = s11/2.;
				const Real s50 = N[8*i+2];
				const Real s51 = F[6*i+3];
				const Real s52 = N[8*i+5];
				const Real s53 = -2*s14*s4*s8;
				const Real s54 = -(s12*s27*s8);
				const Real s55 = s53 + s54;
				const Real s56 = F[6*i+5];
				const Real s57 = s8*s9;
				const Real s58 = N[8*i+7];
				const Real s59 = -2*s20*s57;
				const Real s60 = 2*s21*s8;
				const Real s61 = s59 + s60;
				const Real s62 = s11*s61;
				const Real s63 = -(s12*s32*s8);
				const Real s64 = s62 + s63;
				const Real s65 = s14*s4*s8;
				const Real s66 = -(s12*s8);
				const Real s67 = s65 + s66;
				const Real s68 = 2*s18*s20;
				const Real s69 = -2*s21*s3;
				const Real s70 = s68 + s69;
				const Real s71 = s11*s70;
				const Real s72 = s12*s27*s3;
				const Real s73 = s71 + s72;
				const Real s74 = 2*s14*s3*s9;
				const Real s75 = s12*s3*s32;
				const Real s76 = s74 + s75;
				const Real s77 = s14*s3*s9;
				const Real s78 = -(s12*s3);
				const Real s79 = s77 + s78;
				const Real s80 = 2*s14*s4*s8;
				const Real s81 = s12*s27*s8;
				const Real s82 = s80 + s81;
				const Real s83 = 2*s20*s57;
				const Real s84 = -2*s21*s8;
				const Real s85 = s83 + s84;
				const Real s86 = s11*s85;
				const Real s87 = s12*s32*s8;
				const Real s88 = s86 + s87;
				B[4*i+0] = charge*(-(s12*s3*s35) - s12*s3*s36 + (s11 - s0*s12*s3)*s37 - s12*s2*s3*s38 + s17*s42 + s17*s43 - (s12*s3*s47*s48)/2. + s39*(-0.5*(s12*s3*s40) + s49) - s12*s3*s5*s50 + s29*s51 + s29*s52 + s34*s56 + s34*s58 - s12*s3*s41*s7);
				B[4*i+1] = charge*(s42*s46 + s43*s46 + s51*s55 + s52*s55 + s56*s64 + s58*s64 - s12*s35*s8 - s12*s36*s8 - s0*s12*s37*s8 - s12*s2*s38*s8 - (s12*s39*s40*s8)/2. - s12*s41*s7*s8 + s47*(s49 - (s12*s48*s8)/2.) + s50*(s11 - s12*s5*s8));
				B[4*i+2] = charge*(s12*s3*s35 + s12*s3*s36 + s0*s12*s3*s37 + (s11 + s12*s2*s3)*s38 + (s12*s3*s47*s48)/2. + s39*((s12*s3*s40)/2. + s49) + s12*s3*s5*s50 + s42*s67 + s43*s67 + s12*s3*s41*s7 + s51*s73 + s52*s73 + s56*s76 + s58*s76);
				B[4*i+3] = charge*(s42*s79 + s43*s79 + s12*s35*s8 + s12*s36*s8 + s0*s12*s37*s8 + s12*s2*s38*s8 + (s12*s39*s40*s8)/2. + s12*s5*s50*s8 + s47*(s49 + (s12*s48*s8)/2.) + s41*(s11 + s12*s7*s8) + s51*s82 + s52*s82 + s56*s88 + s58*s88);
			},
			simplices.Dim(0),
			ThreadCount()
		);
    }

	}; // SimplicialMeshDetails<1,2,Real,Int>

//------------------------------------------------------------

    template<typename Real, typename Int, typename LInt>
    struct SimplicialMeshDetails<1,3,Real,Int,LInt>
    {
	public:

		using SparseMatrix_T = Sparse::MatrixCSR<Real,Int,LInt>;

	private:

		const Int thread_count = 1;

		static constexpr Int  SIZE = 2;
		static constexpr Real nth  = Inv<Real>( 2 );

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(Ramp(thread_count_))
		{}
	
		Int ThreadCount() const
		{
			return thread_count;
		}
	
		std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<1,3,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void DNearFarToHulls( 
		cref<Tensor2<Real,Int>> V_coords, 
		cref<Tensor2<Int ,Int>> simplices, 
		cref<Tensor2<Real,Int>> P_D_near, 
		cref<Tensor2<Real,Int>> P_D_far, 
		cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
		mref<Tensor3<Real,Int>> buffer
	) const
    {
        TOOLS_PTIMER(timer,ClassName()+"::DNearFarToHulls");

        if( P_D_near.Dim(1) != 13 )
        {
            eprint("in DNearFarToHulls: P_D_near.Dim(1) != 13. Aborting");
        }

		cptr<Real> X = V_coords.data();
		cptr<Int>  S = simplices.data();
		cptr<Real> N = P_D_near.data();
		cptr<Real> F = P_D_far.data();
		mptr<Real> B = buffer.data();
        
		ParallelDo(
			[=]( const Int i )
			{
				Real charge = 0;
				for( Int k = 0; k < SIZE; ++k )
				{
					charge += V_charges[S[SIZE * i + k]];
				}
				charge *= nth;
				const Real s0 = X[3*S[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = X[3*S[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = X[3*S[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = X[3*S[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = X[3*S[2*i+0]+2];
				const Real s11 = -s10;
				const Real s12 = X[3*S[2*i+1]+2];
				const Real s13 = s11 + s12;
				const Real s14 = s13*s13;
				const Real s15 = s14 + s4 + s9;
				const Real s16 = std::sqrt(s15);
				const Real s17 = s15*s16;
				const Real s18 = 1/s17;
				const Real s19 = 1/s16;
				const Real s20 = -(s18*s4*s8);
				const Real s21 = s19*s8;
				const Real s22 = s20 + s21;
				const Real s23 = -(s13*s18*s4);
				const Real s24 = s13*s19;
				const Real s25 = s23 + s24;
				const Real s26 = s3*s4;
				const Real s27 = s15*s15;
				const Real s28 = 1/s27;
				const Real s29 = 1/s15;
				const Real s30 = -2*s26*s28;
				const Real s31 = 2*s29*s3;
				const Real s32 = s30 + s31;
				const Real s33 = s16*s32;
				const Real s34 = -(s29*s4);
				const Real s35 = 1 + s34;
				const Real s36 = -(s19*s3*s35);
				const Real s37 = s33 + s36;
				const Real s38 = -2*s18*s3*s9;
				const Real s39 = -(s29*s9);
				const Real s40 = 1 + s39;
				const Real s41 = -(s19*s3*s40);
				const Real s42 = s38 + s41;
				const Real s43 = -2*s14*s18*s3;
				const Real s44 = -(s14*s29);
				const Real s45 = 1 + s44;
				const Real s46 = -(s19*s3*s45);
				const Real s47 = s43 + s46;
				const Real s48 = F[10*i+6];
				const Real s49 = N[13*i+9];
				const Real s50 = F[10*i+0];
				const Real s51 = N[13*i+0];
				const Real s52 = N[13*i+1];
				const Real s53 = N[13*i+3];
				const Real s54 = N[13*i+4];
				const Real s55 = F[10*i+1];
				const Real s56 = s0 + s2;
				const Real s57 = N[13*i+5];
				const Real s58 = N[13*i+6];
				const Real s59 = F[10*i+3];
				const Real s60 = s10 + s12;
				const Real s61 = F[10*i+5];
				const Real s62 = N[13*i+8];
				const Real s63 = -(s18*s3*s9);
				const Real s64 = s19*s3;
				const Real s65 = s63 + s64;
				const Real s66 = F[10*i+8];
				const Real s67 = N[13*i+11];
				const Real s68 = -(s13*s18*s9);
				const Real s69 = s24 + s68;
				const Real s70 = F[10*i+2];
				const Real s71 = s5 + s7;
				const Real s72 = s16/2.;
				const Real s73 = N[13*i+2];
				const Real s74 = F[10*i+4];
				const Real s75 = N[13*i+7];
				const Real s76 = -2*s18*s4*s8;
				const Real s77 = -(s19*s35*s8);
				const Real s78 = s76 + s77;
				const Real s79 = F[10*i+7];
				const Real s80 = s8*s9;
				const Real s81 = N[13*i+10];
				const Real s82 = -2*s28*s80;
				const Real s83 = 2*s29*s8;
				const Real s84 = s82 + s83;
				const Real s85 = s16*s84;
				const Real s86 = -(s19*s40*s8);
				const Real s87 = s85 + s86;
				const Real s88 = F[10*i+9];
				const Real s89 = N[13*i+12];
				const Real s90 = -2*s14*s18*s8;
				const Real s91 = -(s19*s45*s8);
				const Real s92 = s90 + s91;
				const Real s93 = -(s14*s18*s3);
				const Real s94 = s64 + s93;
				const Real s95 = -(s14*s18*s8);
				const Real s96 = s21 + s95;
				const Real s97 = -2*s13*s18*s4;
				const Real s98 = -(s13*s19*s35);
				const Real s99 = s97 + s98;
				const Real s100 = -2*s13*s18*s9;
				const Real s101 = -(s13*s19*s40);
				const Real s102 = s100 + s101;
				const Real s103 = s13*s14;
				const Real s104 = -2*s103*s28;
				const Real s105 = 2*s13*s29;
				const Real s106 = s104 + s105;
				const Real s107 = s106*s16;
				const Real s108 = -(s13*s19*s45);
				const Real s109 = s107 + s108;
				const Real s110 = s18*s4*s8;
				const Real s111 = -(s19*s8);
				const Real s112 = s110 + s111;
				const Real s113 = s13*s18*s4;
				const Real s114 = -(s13*s19);
				const Real s115 = s113 + s114;
				const Real s116 = 2*s26*s28;
				const Real s117 = -2*s29*s3;
				const Real s118 = s116 + s117;
				const Real s119 = s118*s16;
				const Real s120 = s19*s3*s35;
				const Real s121 = s119 + s120;
				const Real s122 = 2*s18*s3*s9;
				const Real s123 = s19*s3*s40;
				const Real s124 = s122 + s123;
				const Real s125 = 2*s14*s18*s3;
				const Real s126 = s19*s3*s45;
				const Real s127 = s125 + s126;
				const Real s128 = s18*s3*s9;
				const Real s129 = -(s19*s3);
				const Real s130 = s128 + s129;
				const Real s131 = s13*s18*s9;
				const Real s132 = s114 + s131;
				const Real s133 = 2*s18*s4*s8;
				const Real s134 = s19*s35*s8;
				const Real s135 = s133 + s134;
				const Real s136 = 2*s28*s80;
				const Real s137 = -2*s29*s8;
				const Real s138 = s136 + s137;
				const Real s139 = s138*s16;
				const Real s140 = s19*s40*s8;
				const Real s141 = s139 + s140;
				const Real s142 = 2*s14*s18*s8;
				const Real s143 = s19*s45*s8;
				const Real s144 = s142 + s143;
				const Real s145 = s14*s18*s3;
				const Real s146 = s129 + s145;
				const Real s147 = s14*s18*s8;
				const Real s148 = s111 + s147;
				const Real s149 = 2*s13*s18*s4;
				const Real s150 = s13*s19*s35;
				const Real s151 = s149 + s150;
				const Real s152 = 2*s13*s18*s9;
				const Real s153 = s13*s19*s40;
				const Real s154 = s152 + s153;
				const Real s155 = 2*s103*s28;
				const Real s156 = -2*s13*s29;
				const Real s157 = s155 + s156;
				const Real s158 = s157*s16;
				const Real s159 = s13*s19*s45;
				const Real s160 = s158 + s159;
				B[6*i+0] = charge*(s25*s48 + s25*s49 - s19*s3*s50 - s19*s3*s51 + (s16 - s0*s19*s3)*s52 - s10*s19*s3*s53 - s19*s2*s3*s54 - s12*s19*s3*s58 - (s19*s3*s59*s60)/2. + s22*s61 + s22*s62 - s19*s3*s57*s7 - (s19*s3*s70*s71)/2. + s55*(-0.5*(s19*s3*s56) + s72) - s19*s3*s5*s73 + s37*s74 + s37*s75 + s42*s79 - s13*s18*s3*s66*s8 - s13*s18*s3*s67*s8 + s42*s81 + s47*s88 + s47*s89);
				B[6*i+1] = charge*(s61*s65 + s62*s65 + s66*s69 + s67*s69 + s74*s78 + s75*s78 - s13*s18*s3*s48*s8 - s13*s18*s3*s49*s8 - s19*s50*s8 - s19*s51*s8 - s0*s19*s52*s8 - s10*s19*s53*s8 - s19*s2*s54*s8 - (s19*s55*s56*s8)/2. - s12*s19*s58*s8 - (s19*s59*s60*s8)/2. - s19*s57*s7*s8 + s73*(s16 - s19*s5*s8) + s70*(s72 - (s19*s71*s8)/2.) + s79*s87 + s81*s87 + s88*s92 + s89*s92);
				B[6*i+2] = charge*(-(s13*s19*s50) - s13*s19*s51 - s0*s13*s19*s52 + (s16 - s10*s13*s19)*s53 - s13*s19*s2*s54 - (s13*s19*s55*s56)/2. - s12*s13*s19*s58 - s13*s19*s57*s7 - (s13*s19*s70*s71)/2. + s59*(-0.5*(s13*s19*s60) + s72) - s13*s19*s5*s73 + s102*s79 - s13*s18*s3*s61*s8 - s13*s18*s3*s62*s8 + s102*s81 + s109*s88 + s109*s89 + s48*s94 + s49*s94 + s66*s96 + s67*s96 + s74*s99 + s75*s99);
				B[6*i+3] = charge*(s115*s48 + s115*s49 + s19*s3*s50 + s19*s3*s51 + s0*s19*s3*s52 + s10*s19*s3*s53 + (s16 + s19*s2*s3)*s54 + s12*s19*s3*s58 + (s19*s3*s59*s60)/2. + s112*s61 + s112*s62 + s19*s3*s57*s7 + (s19*s3*s70*s71)/2. + s55*((s19*s3*s56)/2. + s72) + s19*s3*s5*s73 + s121*s74 + s121*s75 + s124*s79 + s13*s18*s3*s66*s8 + s13*s18*s3*s67*s8 + s124*s81 + s127*s88 + s127*s89);
				B[6*i+4] = charge*(s130*s61 + s130*s62 + s132*s66 + s132*s67 + s135*s74 + s135*s75 + s141*s79 + s13*s18*s3*s48*s8 + s13*s18*s3*s49*s8 + s19*s50*s8 + s19*s51*s8 + s0*s19*s52*s8 + s10*s19*s53*s8 + s19*s2*s54*s8 + (s19*s55*s56*s8)/2. + s12*s19*s58*s8 + (s19*s59*s60*s8)/2. + s19*s5*s73*s8 + s57*(s16 + s19*s7*s8) + s70*(s72 + (s19*s71*s8)/2.) + s141*s81 + s144*s88 + s144*s89);
				B[6*i+5] = charge*(s146*s48 + s146*s49 + s13*s19*s50 + s13*s19*s51 + s0*s13*s19*s52 + s10*s13*s19*s53 + s13*s19*s2*s54 + (s13*s19*s55*s56)/2. + (s16 + s12*s13*s19)*s58 + s148*s66 + s148*s67 + s13*s19*s57*s7 + (s13*s19*s70*s71)/2. + s59*((s13*s19*s60)/2. + s72) + s13*s19*s5*s73 + s151*s74 + s151*s75 + s154*s79 + s13*s18*s3*s61*s8 + s13*s18*s3*s62*s8 + s154*s81 + s160*s88 + s160*s89);
			},
			simplices.Dim(0),
			ThreadCount()
		);
    }

	}; // SimplicialMeshDetails<1,3,Real,Int>

//------------------------------------------------------------

    template<typename Real, typename Int, typename LInt>
    struct SimplicialMeshDetails<1,4,Real,Int,LInt>
    {
	public:

		using SparseMatrix_T = Sparse::MatrixCSR<Real,Int,LInt>;

	private:

		const Int thread_count = 1;

		static constexpr Int  SIZE = 2;
		static constexpr Real nth  = Inv<Real>( 2 );

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(Ramp(thread_count_))
		{}
	
		Int ThreadCount() const
		{
			return thread_count;
		}
	
		std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<1,4,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
	void DNearFarToHulls( 
		cref<Tensor2<Real,Int>> V_coords, 
		cref<Tensor2<Int ,Int>> simplices, 
		cref<Tensor2<Real,Int>> P_D_near, 
		cref<Tensor2<Real,Int>> P_D_far, 
		cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
		mref<Tensor3<Real,Int>> buffer
	) const
    {
        TOOLS_PTIMER(timer,ClassName()+"::DNearFarToHulls");

        if( P_D_near.Dim(1) != 19 )
        {
            eprint("in DNearFarToHulls: P_D_near.Dim(1) != 19. Aborting");
        }

		cptr<Real> X = V_coords.data();
		cptr<Int>  S = simplices.data();
		cptr<Real> N = P_D_near.data();
		cptr<Real> F = P_D_far.data();
		mptr<Real> B = buffer.data();
        
		ParallelDo(
			[=]( const Int i )
			{
				Real charge = 0;
				for( Int k = 0; k < SIZE; ++k )
				{
					charge += V_charges[S[SIZE * i + k]];
				}
				charge *= nth;
				const Real s0 = X[4*S[2*i+0]+0];
				const Real s1 = -s0;
				const Real s2 = X[4*S[2*i+1]+0];
				const Real s3 = s1 + s2;
				const Real s4 = s3*s3;
				const Real s5 = X[4*S[2*i+0]+1];
				const Real s6 = -s5;
				const Real s7 = X[4*S[2*i+1]+1];
				const Real s8 = s6 + s7;
				const Real s9 = s8*s8;
				const Real s10 = X[4*S[2*i+0]+2];
				const Real s11 = -s10;
				const Real s12 = X[4*S[2*i+1]+2];
				const Real s13 = s11 + s12;
				const Real s14 = s13*s13;
				const Real s15 = X[4*S[2*i+0]+3];
				const Real s16 = -s15;
				const Real s17 = X[4*S[2*i+1]+3];
				const Real s18 = s16 + s17;
				const Real s19 = s18*s18;
				const Real s20 = s14 + s19 + s4 + s9;
				const Real s21 = std::sqrt(s20);
				const Real s22 = s20*s21;
				const Real s23 = 1/s22;
				const Real s24 = 1/s21;
				const Real s25 = -(s23*s4*s8);
				const Real s26 = s24*s8;
				const Real s27 = s25 + s26;
				const Real s28 = -(s13*s23*s4);
				const Real s29 = s13*s24;
				const Real s30 = s28 + s29;
				const Real s31 = -(s18*s23*s4);
				const Real s32 = s18*s24;
				const Real s33 = s31 + s32;
				const Real s34 = s3*s4;
				const Real s35 = s20*s20;
				const Real s36 = 1/s35;
				const Real s37 = 1/s20;
				const Real s38 = -2*s34*s36;
				const Real s39 = 2*s3*s37;
				const Real s40 = s38 + s39;
				const Real s41 = s21*s40;
				const Real s42 = -(s37*s4);
				const Real s43 = 1 + s42;
				const Real s44 = -(s24*s3*s43);
				const Real s45 = s41 + s44;
				const Real s46 = -2*s23*s3*s9;
				const Real s47 = -(s37*s9);
				const Real s48 = 1 + s47;
				const Real s49 = -(s24*s3*s48);
				const Real s50 = s46 + s49;
				const Real s51 = -2*s14*s23*s3;
				const Real s52 = -(s14*s37);
				const Real s53 = 1 + s52;
				const Real s54 = -(s24*s3*s53);
				const Real s55 = s51 + s54;
				const Real s56 = -2*s19*s23*s3;
				const Real s57 = -(s19*s37);
				const Real s58 = 1 + s57;
				const Real s59 = -(s24*s3*s58);
				const Real s60 = s56 + s59;
				const Real s61 = F[15*i+7];
				const Real s62 = N[19*i+11];
				const Real s63 = F[15*i+8];
				const Real s64 = N[19*i+12];
				const Real s65 = F[15*i+13];
				const Real s66 = N[19*i+17];
				const Real s67 = F[15*i+0];
				const Real s68 = N[19*i+0];
				const Real s69 = N[19*i+1];
				const Real s70 = N[19*i+3];
				const Real s71 = N[19*i+4];
				const Real s72 = N[19*i+5];
				const Real s73 = F[15*i+1];
				const Real s74 = s0 + s2;
				const Real s75 = N[19*i+6];
				const Real s76 = N[19*i+7];
				const Real s77 = F[15*i+3];
				const Real s78 = s10 + s12;
				const Real s79 = N[19*i+8];
				const Real s80 = F[15*i+4];
				const Real s81 = s15 + s17;
				const Real s82 = F[15*i+6];
				const Real s83 = N[19*i+10];
				const Real s84 = -(s23*s3*s9);
				const Real s85 = s24*s3;
				const Real s86 = s84 + s85;
				const Real s87 = F[15*i+10];
				const Real s88 = N[19*i+14];
				const Real s89 = -(s13*s23*s9);
				const Real s90 = s29 + s89;
				const Real s91 = F[15*i+11];
				const Real s92 = N[19*i+15];
				const Real s93 = -(s18*s23*s9);
				const Real s94 = s32 + s93;
				const Real s95 = F[15*i+2];
				const Real s96 = s5 + s7;
				const Real s97 = s21/2.;
				const Real s98 = N[19*i+2];
				const Real s99 = F[15*i+5];
				const Real s100 = N[19*i+9];
				const Real s101 = -2*s23*s4*s8;
				const Real s102 = -(s24*s43*s8);
				const Real s103 = s101 + s102;
				const Real s104 = F[15*i+9];
				const Real s105 = s8*s9;
				const Real s106 = N[19*i+13];
				const Real s107 = -2*s105*s36;
				const Real s108 = 2*s37*s8;
				const Real s109 = s107 + s108;
				const Real s110 = s109*s21;
				const Real s111 = -(s24*s48*s8);
				const Real s112 = s110 + s111;
				const Real s113 = F[15*i+12];
				const Real s114 = N[19*i+16];
				const Real s115 = -2*s14*s23*s8;
				const Real s116 = -(s24*s53*s8);
				const Real s117 = s115 + s116;
				const Real s118 = F[15*i+14];
				const Real s119 = N[19*i+18];
				const Real s120 = -2*s19*s23*s8;
				const Real s121 = -(s24*s58*s8);
				const Real s122 = s120 + s121;
				const Real s123 = -(s14*s23*s3);
				const Real s124 = s123 + s85;
				const Real s125 = -(s14*s23*s8);
				const Real s126 = s125 + s26;
				const Real s127 = -(s14*s18*s23);
				const Real s128 = s127 + s32;
				const Real s129 = -2*s13*s23*s4;
				const Real s130 = -(s13*s24*s43);
				const Real s131 = s129 + s130;
				const Real s132 = -2*s13*s23*s9;
				const Real s133 = -(s13*s24*s48);
				const Real s134 = s132 + s133;
				const Real s135 = s13*s14;
				const Real s136 = -2*s135*s36;
				const Real s137 = 2*s13*s37;
				const Real s138 = s136 + s137;
				const Real s139 = s138*s21;
				const Real s140 = -(s13*s24*s53);
				const Real s141 = s139 + s140;
				const Real s142 = -2*s13*s19*s23;
				const Real s143 = -(s13*s24*s58);
				const Real s144 = s142 + s143;
				const Real s145 = -(s19*s23*s3);
				const Real s146 = s145 + s85;
				const Real s147 = -(s19*s23*s8);
				const Real s148 = s147 + s26;
				const Real s149 = -(s13*s19*s23);
				const Real s150 = s149 + s29;
				const Real s151 = -2*s18*s23*s4;
				const Real s152 = -(s18*s24*s43);
				const Real s153 = s151 + s152;
				const Real s154 = -2*s18*s23*s9;
				const Real s155 = -(s18*s24*s48);
				const Real s156 = s154 + s155;
				const Real s157 = -2*s14*s18*s23;
				const Real s158 = -(s18*s24*s53);
				const Real s159 = s157 + s158;
				const Real s160 = s18*s19;
				const Real s161 = -2*s160*s36;
				const Real s162 = 2*s18*s37;
				const Real s163 = s161 + s162;
				const Real s164 = s163*s21;
				const Real s165 = -(s18*s24*s58);
				const Real s166 = s164 + s165;
				const Real s167 = s23*s4*s8;
				const Real s168 = -(s24*s8);
				const Real s169 = s167 + s168;
				const Real s170 = s13*s23*s4;
				const Real s171 = -(s13*s24);
				const Real s172 = s170 + s171;
				const Real s173 = s18*s23*s4;
				const Real s174 = -(s18*s24);
				const Real s175 = s173 + s174;
				const Real s176 = 2*s34*s36;
				const Real s177 = -2*s3*s37;
				const Real s178 = s176 + s177;
				const Real s179 = s178*s21;
				const Real s180 = s24*s3*s43;
				const Real s181 = s179 + s180;
				const Real s182 = 2*s23*s3*s9;
				const Real s183 = s24*s3*s48;
				const Real s184 = s182 + s183;
				const Real s185 = 2*s14*s23*s3;
				const Real s186 = s24*s3*s53;
				const Real s187 = s185 + s186;
				const Real s188 = 2*s19*s23*s3;
				const Real s189 = s24*s3*s58;
				const Real s190 = s188 + s189;
				const Real s191 = s23*s3*s9;
				const Real s192 = -(s24*s3);
				const Real s193 = s191 + s192;
				const Real s194 = s13*s23*s9;
				const Real s195 = s171 + s194;
				const Real s196 = s18*s23*s9;
				const Real s197 = s174 + s196;
				const Real s198 = 2*s23*s4*s8;
				const Real s199 = s24*s43*s8;
				const Real s200 = s198 + s199;
				const Real s201 = 2*s105*s36;
				const Real s202 = -2*s37*s8;
				const Real s203 = s201 + s202;
				const Real s204 = s203*s21;
				const Real s205 = s24*s48*s8;
				const Real s206 = s204 + s205;
				const Real s207 = 2*s14*s23*s8;
				const Real s208 = s24*s53*s8;
				const Real s209 = s207 + s208;
				const Real s210 = 2*s19*s23*s8;
				const Real s211 = s24*s58*s8;
				const Real s212 = s210 + s211;
				const Real s213 = s14*s23*s3;
				const Real s214 = s192 + s213;
				const Real s215 = s14*s23*s8;
				const Real s216 = s168 + s215;
				const Real s217 = s14*s18*s23;
				const Real s218 = s174 + s217;
				const Real s219 = 2*s13*s23*s4;
				const Real s220 = s13*s24*s43;
				const Real s221 = s219 + s220;
				const Real s222 = 2*s13*s23*s9;
				const Real s223 = s13*s24*s48;
				const Real s224 = s222 + s223;
				const Real s225 = 2*s135*s36;
				const Real s226 = -2*s13*s37;
				const Real s227 = s225 + s226;
				const Real s228 = s21*s227;
				const Real s229 = s13*s24*s53;
				const Real s230 = s228 + s229;
				const Real s231 = 2*s13*s19*s23;
				const Real s232 = s13*s24*s58;
				const Real s233 = s231 + s232;
				const Real s234 = s19*s23*s3;
				const Real s235 = s192 + s234;
				const Real s236 = s19*s23*s8;
				const Real s237 = s168 + s236;
				const Real s238 = s13*s19*s23;
				const Real s239 = s171 + s238;
				const Real s240 = 2*s18*s23*s4;
				const Real s241 = s18*s24*s43;
				const Real s242 = s240 + s241;
				const Real s243 = 2*s18*s23*s9;
				const Real s244 = s18*s24*s48;
				const Real s245 = s243 + s244;
				const Real s246 = 2*s14*s18*s23;
				const Real s247 = s18*s24*s53;
				const Real s248 = s246 + s247;
				const Real s249 = 2*s160*s36;
				const Real s250 = -2*s18*s37;
				const Real s251 = s249 + s250;
				const Real s252 = s21*s251;
				const Real s253 = s18*s24*s58;
				const Real s254 = s252 + s253;
				B[8*i+0] = charge*(s100*s45 + s104*s50 + s106*s50 + s113*s55 + s114*s55 + s118*s60 + s119*s60 + s30*s61 + s30*s62 + s33*s63 + s33*s64 - s13*s18*s23*s3*s65 - s13*s18*s23*s3*s66 - s24*s3*s67 - s24*s3*s68 + (s21 - s0*s24*s3)*s69 - s10*s24*s3*s70 - s15*s24*s3*s71 - s2*s24*s3*s72 - s24*s3*s7*s75 - s12*s24*s3*s76 - (s24*s3*s77*s78)/2. - s17*s24*s3*s79 - (s24*s3*s80*s81)/2. + s27*s82 + s27*s83 - s13*s23*s3*s8*s87 - s13*s23*s3*s8*s88 - s18*s23*s3*s8*s91 - s18*s23*s3*s8*s92 - (s24*s3*s95*s96)/2. + s73*(-0.5*(s24*s3*s74) + s97) - s24*s3*s5*s98 + s45*s99);
				B[8*i+1] = charge*(s100*s103 + s104*s112 + s106*s112 + s113*s117 + s114*s117 + s118*s122 + s119*s122 - s13*s23*s3*s61*s8 - s13*s23*s3*s62*s8 - s18*s23*s3*s63*s8 - s18*s23*s3*s64*s8 - s13*s18*s23*s65*s8 - s13*s18*s23*s66*s8 - s24*s67*s8 - s24*s68*s8 - s0*s24*s69*s8 - s10*s24*s70*s8 - s15*s24*s71*s8 - s2*s24*s72*s8 - (s24*s73*s74*s8)/2. - s24*s7*s75*s8 - s12*s24*s76*s8 - (s24*s77*s78*s8)/2. - s17*s24*s79*s8 - (s24*s8*s80*s81)/2. + s82*s86 + s83*s86 + s87*s90 + s88*s90 + s91*s94 + s92*s94 + s95*(-0.5*(s24*s8*s96) + s97) + (s21 - s24*s5*s8)*s98 + s103*s99);
				B[8*i+2] = charge*(s100*s131 + s104*s134 + s106*s134 + s113*s141 + s114*s141 + s118*s144 + s119*s144 + s124*s61 + s124*s62 - s13*s18*s23*s3*s63 - s13*s18*s23*s3*s64 + s128*s65 + s128*s66 - s13*s24*s67 - s13*s24*s68 - s0*s13*s24*s69 + (s21 - s10*s13*s24)*s70 - s13*s15*s24*s71 - s13*s2*s24*s72 - (s13*s24*s73*s74)/2. - s13*s24*s7*s75 - s12*s13*s24*s76 - s13*s17*s24*s79 - (s13*s24*s80*s81)/2. - s13*s23*s3*s8*s82 - s13*s23*s3*s8*s83 + s126*s87 + s126*s88 - s13*s18*s23*s8*s91 - s13*s18*s23*s8*s92 - (s13*s24*s95*s96)/2. + s77*(-0.5*(s13*s24*s78) + s97) - s13*s24*s5*s98 + s131*s99);
				B[8*i+3] = charge*(s100*s153 + s104*s156 + s106*s156 + s113*s159 + s114*s159 + s118*s166 + s119*s166 - s13*s18*s23*s3*s61 - s13*s18*s23*s3*s62 + s146*s63 + s146*s64 + s150*s65 + s150*s66 - s18*s24*s67 - s18*s24*s68 - s0*s18*s24*s69 - s10*s18*s24*s70 + (s21 - s15*s18*s24)*s71 - s18*s2*s24*s72 - (s18*s24*s73*s74)/2. - s18*s24*s7*s75 - s12*s18*s24*s76 - (s18*s24*s77*s78)/2. - s17*s18*s24*s79 - s18*s23*s3*s8*s82 - s18*s23*s3*s8*s83 - s13*s18*s23*s8*s87 - s13*s18*s23*s8*s88 + s148*s91 + s148*s92 - (s18*s24*s95*s96)/2. + s80*(-0.5*(s18*s24*s81) + s97) - s18*s24*s5*s98 + s153*s99);
				B[8*i+4] = charge*(s100*s181 + s104*s184 + s106*s184 + s113*s187 + s114*s187 + s118*s190 + s119*s190 + s172*s61 + s172*s62 + s175*s63 + s175*s64 + s13*s18*s23*s3*s65 + s13*s18*s23*s3*s66 + s24*s3*s67 + s24*s3*s68 + s0*s24*s3*s69 + s10*s24*s3*s70 + s15*s24*s3*s71 + (s21 + s2*s24*s3)*s72 + s24*s3*s7*s75 + s12*s24*s3*s76 + (s24*s3*s77*s78)/2. + s17*s24*s3*s79 + (s24*s3*s80*s81)/2. + s169*s82 + s169*s83 + s13*s23*s3*s8*s87 + s13*s23*s3*s8*s88 + s18*s23*s3*s8*s91 + s18*s23*s3*s8*s92 + (s24*s3*s95*s96)/2. + s73*((s24*s3*s74)/2. + s97) + s24*s3*s5*s98 + s181*s99);
				B[8*i+5] = charge*(s100*s200 + s104*s206 + s106*s206 + s113*s209 + s114*s209 + s118*s212 + s119*s212 + s13*s23*s3*s61*s8 + s13*s23*s3*s62*s8 + s18*s23*s3*s63*s8 + s18*s23*s3*s64*s8 + s13*s18*s23*s65*s8 + s13*s18*s23*s66*s8 + s24*s67*s8 + s24*s68*s8 + s0*s24*s69*s8 + s10*s24*s70*s8 + s15*s24*s71*s8 + s2*s24*s72*s8 + (s24*s73*s74*s8)/2. + s12*s24*s76*s8 + (s24*s77*s78*s8)/2. + s17*s24*s79*s8 + s75*(s21 + s24*s7*s8) + (s24*s8*s80*s81)/2. + s193*s82 + s193*s83 + s195*s87 + s195*s88 + s197*s91 + s197*s92 + s95*((s24*s8*s96)/2. + s97) + s24*s5*s8*s98 + s200*s99);
				B[8*i+6] = charge*(s100*s221 + s104*s224 + s106*s224 + s113*s230 + s114*s230 + s118*s233 + s119*s233 + s214*s61 + s214*s62 + s13*s18*s23*s3*s63 + s13*s18*s23*s3*s64 + s218*s65 + s218*s66 + s13*s24*s67 + s13*s24*s68 + s0*s13*s24*s69 + s10*s13*s24*s70 + s13*s15*s24*s71 + s13*s2*s24*s72 + (s13*s24*s73*s74)/2. + s13*s24*s7*s75 + (s21 + s12*s13*s24)*s76 + s13*s17*s24*s79 + (s13*s24*s80*s81)/2. + s13*s23*s3*s8*s82 + s13*s23*s3*s8*s83 + s216*s87 + s216*s88 + s13*s18*s23*s8*s91 + s13*s18*s23*s8*s92 + (s13*s24*s95*s96)/2. + s77*((s13*s24*s78)/2. + s97) + s13*s24*s5*s98 + s221*s99);
				B[8*i+7] = charge*(s100*s242 + s104*s245 + s106*s245 + s113*s248 + s114*s248 + s118*s254 + s119*s254 + s13*s18*s23*s3*s61 + s13*s18*s23*s3*s62 + s235*s63 + s235*s64 + s239*s65 + s239*s66 + s18*s24*s67 + s18*s24*s68 + s0*s18*s24*s69 + s10*s18*s24*s70 + s15*s18*s24*s71 + s18*s2*s24*s72 + (s18*s24*s73*s74)/2. + s18*s24*s7*s75 + s12*s18*s24*s76 + (s18*s24*s77*s78)/2. + (s21 + s17*s18*s24)*s79 + s18*s23*s3*s8*s82 + s18*s23*s3*s8*s83 + s13*s18*s23*s8*s87 + s13*s18*s23*s8*s88 + s237*s91 + s237*s92 + (s18*s24*s95*s96)/2. + s80*((s18*s24*s81)/2. + s97) + s18*s24*s5*s98 + s242*s99);
			},
			simplices.Dim(0),
			ThreadCount()
		);
    }

	}; // SimplicialMeshDetails<1,4,Real,Int>

//------------------------------------------------------------

    template<typename Real, typename Int, typename LInt>
    struct SimplicialMeshDetails<2,3,Real,Int,LInt>
    {
	public:

		using SparseMatrix_T = Sparse::MatrixCSR<Real,Int,LInt>;

	private:

		const Int thread_count = 1;

		static constexpr Int  SIZE = 3;
		static constexpr Real nth  = Inv<Real>( 3 );

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(Ramp(thread_count_))
		{}
	
		Int ThreadCount() const
		{
			return thread_count;
		}
	
		std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<2,3,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
    void DNearFarToHulls( 
		cref<Tensor2<Real,Int>> V_coords, 
		cref<Tensor2<Int ,Int>> simplices, 
		cref<Tensor2<Real,Int>> P_D_near, 
		cref<Tensor2<Real,Int>> P_D_far, 
		cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
		mref<Tensor3<Real,Int>> buffer
	) const
    {
        TOOLS_PTIMER(timer,ClassName()+"::DNearFarToHulls");

        if( P_D_near.Dim(1) != 16 )
        {
            eprint("in DNearFarToHulls: P_D_near.Dim(1) != 16. Aborting");
        }
        
		cptr<Real> X = V_coords.data();
        cptr<Int>  S = simplices.data();
		cptr<Real> N = P_D_near.data();
		cptr<Real> F = P_D_far.data();
        mptr<Real> B = buffer.data();

		ParallelDo(
			[=]( const Int i )
			{
				Real charge = 0;
				for( Int k = 0; k < SIZE; ++k )
				{
					charge += V_charges[S[SIZE * i + k]];
				}
				charge *= nth;

				const Real s0 = X[3*S[3*i+0]+2];
				const Real s1 = X[3*S[3*i+1]+1];
				const Real s2 = -(s0*s1);
				const Real s3 = X[3*S[3*i+0]+1];
				const Real s4 = X[3*S[3*i+1]+2];
				const Real s5 = s3*s4;
				const Real s6 = X[3*S[3*i+2]+1];
				const Real s7 = s0*s6;
				const Real s8 = -(s4*s6);
				const Real s9 = X[3*S[3*i+2]+2];
				const Real s10 = -(s3*s9);
				const Real s11 = s1*s9;
				const Real s12 = s10 + s11 + s2 + s5 + s7 + s8;
				const Real s13 = s12*s12;
				const Real s14 = X[3*S[3*i+2]+0];
				const Real s15 = X[3*S[3*i+0]+0];
				const Real s16 = X[3*S[3*i+1]+0];
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
				const Real s34 = std::sqrt(s33);
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
				const Real s45 = -0.25*(s24*s36*s43);
				const Real s46 = s23*s38*s44;
				const Real s47 = s45 + s46;
				const Real s48 = -0.25*(s23*s31*s36*s43);
				const Real s49 = (s23*s41*s44)/2.;
				const Real s50 = (s31*s38*s44)/2.;
				const Real s51 = s48 + s49 + s50;
				const Real s52 = -0.25*(s32*s36*s43);
				const Real s53 = s31*s41*s44;
				const Real s54 = s52 + s53;
				const Real s55 = -0.25*(s12*s23*s36*s43);
				const Real s56 = (s12*s38*s44)/2.;
				const Real s57 = s55 + s56;
				const Real s58 = -0.25*(s12*s31*s36*s43);
				const Real s59 = (s12*s41*s44)/2.;
				const Real s60 = s58 + s59;
				const Real s61 = F[10*i+7];
				const Real s62 = N[16*i+13];
				const Real s63 = -s16;
				const Real s64 = s14 + s63;
				const Real s65 = 2*s23*s64;
				const Real s66 = -s9;
				const Real s67 = s4 + s66;
				const Real s68 = 2*s12*s67;
				const Real s69 = s65 + s68;
				const Real s70 = F[10*i+0];
				const Real s71 = N[16*i+0];
				const Real s72 = N[16*i+1];
				const Real s73 = N[16*i+3];
				const Real s74 = N[16*i+4];
				const Real s75 = N[16*i+5];
				const Real s76 = N[16*i+6];
				const Real s77 = N[16*i+7];
				const Real s78 = F[10*i+1];
				const Real s79 = s14 + s15 + s16;
				const Real s80 = N[16*i+8];
				const Real s81 = N[16*i+9];
				const Real s82 = F[10*i+3];
				const Real s83 = s0 + s4 + s9;
				const Real s84 = F[10*i+9];
				const Real s85 = N[16*i+15];
				const Real s86 = -0.25*(s24*s36*s69);
				const Real s87 = s23*s44*s64;
				const Real s88 = s86 + s87;
				const Real s89 = F[10*i+8];
				const Real s90 = N[16*i+14];
				const Real s91 = -0.25*(s23*s31*s36*s69);
				const Real s92 = (s31*s44*s64)/2.;
				const Real s93 = s91 + s92;
				const Real s94 = F[10*i+5];
				const Real s95 = N[16*i+11];
				const Real s96 = -0.25*(s12*s31*s36*s69);
				const Real s97 = (s31*s44*s67)/2.;
				const Real s98 = s96 + s97;
				const Real s99 = F[10*i+6];
				const Real s100 = N[16*i+12];
				const Real s101 = -0.25*(s12*s23*s36*s69);
				const Real s102 = (s23*s44*s67)/2.;
				const Real s103 = (s12*s44*s64)/2.;
				const Real s104 = s101 + s102 + s103;
				const Real s105 = F[10*i+4];
				const Real s106 = N[16*i+10];
				const Real s107 = -0.25*(s13*s36*s69);
				const Real s108 = s12*s44*s67;
				const Real s109 = s107 + s108;
				const Real s110 = F[10*i+2];
				const Real s111 = s1 + s3 + s6;
				const Real s112 = s34/6.;
				const Real s113 = N[16*i+2];
				const Real s114 = s34/2.;
				const Real s115 = -s14;
				const Real s116 = s115 + s16;
				const Real s117 = 2*s116*s31;
				const Real s118 = -s1;
				const Real s119 = s118 + s6;
				const Real s120 = 2*s119*s12;
				const Real s121 = s117 + s120;
				const Real s122 = -0.25*(s121*s23*s31*s36);
				const Real s123 = (s116*s23*s44)/2.;
				const Real s124 = s122 + s123;
				const Real s125 = -0.25*(s12*s121*s23*s36);
				const Real s126 = (s119*s23*s44)/2.;
				const Real s127 = s125 + s126;
				const Real s128 = -0.25*(s121*s32*s36);
				const Real s129 = s116*s31*s44;
				const Real s130 = s128 + s129;
				const Real s131 = -0.25*(s12*s121*s31*s36);
				const Real s132 = (s119*s31*s44)/2.;
				const Real s133 = (s116*s12*s44)/2.;
				const Real s134 = s131 + s132 + s133;
				const Real s135 = -0.25*(s121*s13*s36);
				const Real s136 = s119*s12*s44;
				const Real s137 = s135 + s136;
				const Real s138 = -s3;
				const Real s139 = s138 + s6;
				const Real s140 = 2*s139*s23;
				const Real s141 = s0 + s66;
				const Real s142 = 2*s141*s31;
				const Real s143 = s140 + s142;
				const Real s144 = -0.25*(s143*s24*s36);
				const Real s145 = s139*s23*s44;
				const Real s146 = s144 + s145;
				const Real s147 = -0.25*(s143*s23*s31*s36);
				const Real s148 = (s141*s23*s44)/2.;
				const Real s149 = (s139*s31*s44)/2.;
				const Real s150 = s147 + s148 + s149;
				const Real s151 = -0.25*(s143*s32*s36);
				const Real s152 = s141*s31*s44;
				const Real s153 = s151 + s152;
				const Real s154 = -0.25*(s12*s143*s23*s36);
				const Real s155 = (s12*s139*s44)/2.;
				const Real s156 = s154 + s155;
				const Real s157 = -0.25*(s12*s143*s31*s36);
				const Real s158 = (s12*s141*s44)/2.;
				const Real s159 = s157 + s158;
				const Real s160 = s115 + s15;
				const Real s161 = 2*s160*s23;
				const Real s162 = -s0;
				const Real s163 = s162 + s9;
				const Real s164 = 2*s12*s163;
				const Real s165 = s161 + s164;
				const Real s166 = -0.25*(s165*s24*s36);
				const Real s167 = s160*s23*s44;
				const Real s168 = s166 + s167;
				const Real s169 = -0.25*(s165*s23*s31*s36);
				const Real s170 = (s160*s31*s44)/2.;
				const Real s171 = s169 + s170;
				const Real s172 = -0.25*(s12*s165*s31*s36);
				const Real s173 = (s163*s31*s44)/2.;
				const Real s174 = s172 + s173;
				const Real s175 = -0.25*(s12*s165*s23*s36);
				const Real s176 = (s163*s23*s44)/2.;
				const Real s177 = (s12*s160*s44)/2.;
				const Real s178 = s175 + s176 + s177;
				const Real s179 = -0.25*(s13*s165*s36);
				const Real s180 = s12*s163*s44;
				const Real s181 = s179 + s180;
				const Real s182 = -s15;
				const Real s183 = s14 + s182;
				const Real s184 = 2*s183*s31;
				const Real s185 = s3 + s37;
				const Real s186 = 2*s12*s185;
				const Real s187 = s184 + s186;
				const Real s188 = -0.25*(s187*s23*s31*s36);
				const Real s189 = (s183*s23*s44)/2.;
				const Real s190 = s188 + s189;
				const Real s191 = -0.25*(s12*s187*s23*s36);
				const Real s192 = (s185*s23*s44)/2.;
				const Real s193 = s191 + s192;
				const Real s194 = -0.25*(s187*s32*s36);
				const Real s195 = s183*s31*s44;
				const Real s196 = s194 + s195;
				const Real s197 = -0.25*(s12*s187*s31*s36);
				const Real s198 = (s185*s31*s44)/2.;
				const Real s199 = (s12*s183*s44)/2.;
				const Real s200 = s197 + s198 + s199;
				const Real s201 = -0.25*(s13*s187*s36);
				const Real s202 = s12*s185*s44;
				const Real s203 = s201 + s202;
				const Real s204 = s118 + s3;
				const Real s205 = 2*s204*s23;
				const Real s206 = s162 + s4;
				const Real s207 = 2*s206*s31;
				const Real s208 = s205 + s207;
				const Real s209 = -0.25*(s208*s24*s36);
				const Real s210 = s204*s23*s44;
				const Real s211 = s209 + s210;
				const Real s212 = -0.25*(s208*s23*s31*s36);
				const Real s213 = (s206*s23*s44)/2.;
				const Real s214 = (s204*s31*s44)/2.;
				const Real s215 = s212 + s213 + s214;
				const Real s216 = -0.25*(s208*s32*s36);
				const Real s217 = s206*s31*s44;
				const Real s218 = s216 + s217;
				const Real s219 = -0.25*(s12*s208*s23*s36);
				const Real s220 = (s12*s204*s44)/2.;
				const Real s221 = s219 + s220;
				const Real s222 = -0.25*(s12*s208*s31*s36);
				const Real s223 = (s12*s206*s44)/2.;
				const Real s224 = s222 + s223;
				const Real s225 = s16 + s182;
				const Real s226 = 2*s225*s23;
				const Real s227 = s0 + s40;
				const Real s228 = 2*s12*s227;
				const Real s229 = s226 + s228;
				const Real s230 = -0.25*(s229*s24*s36);
				const Real s231 = s225*s23*s44;
				const Real s232 = s230 + s231;
				const Real s233 = -0.25*(s229*s23*s31*s36);
				const Real s234 = (s225*s31*s44)/2.;
				const Real s235 = s233 + s234;
				const Real s236 = -0.25*(s12*s229*s31*s36);
				const Real s237 = (s227*s31*s44)/2.;
				const Real s238 = s236 + s237;
				const Real s239 = -0.25*(s12*s229*s23*s36);
				const Real s240 = (s227*s23*s44)/2.;
				const Real s241 = (s12*s225*s44)/2.;
				const Real s242 = s239 + s240 + s241;
				const Real s243 = -0.25*(s13*s229*s36);
				const Real s244 = s12*s227*s44;
				const Real s245 = s243 + s244;
				const Real s246 = s15 + s63;
				const Real s247 = 2*s246*s31;
				const Real s248 = s1 + s138;
				const Real s249 = 2*s12*s248;
				const Real s250 = s247 + s249;
				const Real s251 = -0.25*(s23*s250*s31*s36);
				const Real s252 = (s23*s246*s44)/2.;
				const Real s253 = s251 + s252;
				const Real s254 = -0.25*(s12*s23*s250*s36);
				const Real s255 = (s23*s248*s44)/2.;
				const Real s256 = s254 + s255;
				const Real s257 = -0.25*(s250*s32*s36);
				const Real s258 = s246*s31*s44;
				const Real s259 = s257 + s258;
				const Real s260 = -0.25*(s12*s250*s31*s36);
				const Real s261 = (s248*s31*s44)/2.;
				const Real s262 = (s12*s246*s44)/2.;
				const Real s263 = s260 + s261 + s262;
				const Real s264 = -0.25*(s13*s250*s36);
				const Real s265 = s12*s248*s44;
				const Real s266 = s264 + s265;
				B[9*i+0] = charge*(-0.25*(s105*s13*s36*s43) - (s106*s13*s36*s43)/4. + (s110*s111*s43*s44)/12. + (s113*s3*s43*s44)/4. + s100*s57 + s54*s61 + s54*s62 + (s43*s44*s70)/4. + (s43*s44*s71)/4. + (s114 + (s15*s43*s44)/4.)*s72 + (s0*s43*s44*s73)/4. + (s16*s43*s44*s74)/4. + (s1*s43*s44*s75)/4. + (s4*s43*s44*s76)/4. + (s14*s43*s44*s77)/4. + s78*(s112 + (s43*s44*s79)/12.) + (s43*s44*s6*s80)/4. + (s43*s44*s82*s83)/12. + s47*s84 + s47*s85 + s51*s89 + (s43*s44*s81*s9)/4. + s51*s90 + s60*s94 + s60*s95 + s57*s99);
				B[9*i+1] = charge*(s100*s104 + s105*s109 + s106*s109 - (s32*s36*s61*s69)/4. - (s32*s36*s62*s69)/4. + s110*(s112 + (s111*s44*s69)/12.) + s113*(s114 + (s3*s44*s69)/4.) + (s44*s69*s70)/4. + (s44*s69*s71)/4. + (s15*s44*s69*s72)/4. + (s0*s44*s69*s73)/4. + (s16*s44*s69*s74)/4. + (s1*s44*s69*s75)/4. + (s4*s44*s69*s76)/4. + (s14*s44*s69*s77)/4. + (s44*s69*s78*s79)/12. + (s44*s6*s69*s80)/4. + (s44*s69*s82*s83)/12. + s84*s88 + s85*s88 + (s44*s69*s81*s9)/4. + s89*s93 + s90*s93 + s94*s98 + s95*s98 + s104*s99);
				B[9*i+2] = charge*(s100*s127 + s105*s137 + s106*s137 + (s110*s111*s121*s44)/12. + (s113*s121*s3*s44)/4. + s130*s61 + s130*s62 + (s121*s44*s70)/4. + (s121*s44*s71)/4. + (s121*s15*s44*s72)/4. + (s114 + (s0*s121*s44)/4.)*s73 + (s121*s16*s44*s74)/4. + (s1*s121*s44*s75)/4. + (s121*s4*s44*s76)/4. + (s121*s14*s44*s77)/4. + (s121*s44*s78*s79)/12. + (s121*s44*s6*s80)/4. + s82*(s112 + (s121*s44*s83)/12.) - (s121*s24*s36*s84)/4. - (s121*s24*s36*s85)/4. + s124*s89 + (s121*s44*s81*s9)/4. + s124*s90 + s134*s94 + s134*s95 + s127*s99);
				B[9*i+3] = charge*(s100*s156 - (s105*s13*s143*s36)/4. - (s106*s13*s143*s36)/4. + (s110*s111*s143*s44)/12. + (s113*s143*s3*s44)/4. + s153*s61 + s153*s62 + (s143*s44*s70)/4. + (s143*s44*s71)/4. + (s143*s15*s44*s72)/4. + (s0*s143*s44*s73)/4. + (s114 + (s143*s16*s44)/4.)*s74 + (s1*s143*s44*s75)/4. + (s143*s4*s44*s76)/4. + (s14*s143*s44*s77)/4. + s78*(s112 + (s143*s44*s79)/12.) + (s143*s44*s6*s80)/4. + (s143*s44*s82*s83)/12. + s146*s84 + s146*s85 + s150*s89 + (s143*s44*s81*s9)/4. + s150*s90 + s159*s94 + s159*s95 + s156*s99);
				B[9*i+4] = charge*(s100*s178 + s105*s181 + s106*s181 + (s113*s165*s3*s44)/4. + s110*(s112 + (s111*s165*s44)/12.) - (s165*s32*s36*s61)/4. - (s165*s32*s36*s62)/4. + (s165*s44*s70)/4. + (s165*s44*s71)/4. + (s15*s165*s44*s72)/4. + (s0*s165*s44*s73)/4. + (s16*s165*s44*s74)/4. + (s114 + (s1*s165*s44)/4.)*s75 + (s165*s4*s44*s76)/4. + (s14*s165*s44*s77)/4. + (s165*s44*s78*s79)/12. + (s165*s44*s6*s80)/4. + (s165*s44*s82*s83)/12. + s168*s84 + s168*s85 + s171*s89 + (s165*s44*s81*s9)/4. + s171*s90 + s174*s94 + s174*s95 + s178*s99);
				B[9*i+5] = charge*(s100*s193 + s105*s203 + s106*s203 + (s110*s111*s187*s44)/12. + (s113*s187*s3*s44)/4. + s196*s61 + s196*s62 + (s187*s44*s70)/4. + (s187*s44*s71)/4. + (s15*s187*s44*s72)/4. + (s0*s187*s44*s73)/4. + (s16*s187*s44*s74)/4. + (s1*s187*s44*s75)/4. + (s114 + (s187*s4*s44)/4.)*s76 + (s14*s187*s44*s77)/4. + (s187*s44*s78*s79)/12. + (s187*s44*s6*s80)/4. + s82*(s112 + (s187*s44*s83)/12.) - (s187*s24*s36*s84)/4. - (s187*s24*s36*s85)/4. + s190*s89 + (s187*s44*s81*s9)/4. + s190*s90 + s200*s94 + s200*s95 + s193*s99);
				B[9*i+6] = charge*(s100*s221 - (s105*s13*s208*s36)/4. - (s106*s13*s208*s36)/4. + (s110*s111*s208*s44)/12. + (s113*s208*s3*s44)/4. + s218*s61 + s218*s62 + (s208*s44*s70)/4. + (s208*s44*s71)/4. + (s15*s208*s44*s72)/4. + (s0*s208*s44*s73)/4. + (s16*s208*s44*s74)/4. + (s1*s208*s44*s75)/4. + (s208*s4*s44*s76)/4. + (s114 + (s14*s208*s44)/4.)*s77 + s78*(s112 + (s208*s44*s79)/12.) + (s208*s44*s6*s80)/4. + (s208*s44*s82*s83)/12. + s211*s84 + s211*s85 + s215*s89 + (s208*s44*s81*s9)/4. + s215*s90 + s224*s94 + s224*s95 + s221*s99);
				B[9*i+7] = charge*(s100*s242 + s105*s245 + s106*s245 + (s113*s229*s3*s44)/4. + s110*(s112 + (s111*s229*s44)/12.) - (s229*s32*s36*s61)/4. - (s229*s32*s36*s62)/4. + (s229*s44*s70)/4. + (s229*s44*s71)/4. + (s15*s229*s44*s72)/4. + (s0*s229*s44*s73)/4. + (s16*s229*s44*s74)/4. + (s1*s229*s44*s75)/4. + (s229*s4*s44*s76)/4. + (s14*s229*s44*s77)/4. + (s229*s44*s78*s79)/12. + (s114 + (s229*s44*s6)/4.)*s80 + (s229*s44*s82*s83)/12. + s232*s84 + s232*s85 + s235*s89 + (s229*s44*s81*s9)/4. + s235*s90 + s238*s94 + s238*s95 + s242*s99);
				B[9*i+8] = charge*(s100*s256 + s105*s266 + s106*s266 + (s110*s111*s250*s44)/12. + (s113*s250*s3*s44)/4. + s259*s61 + s259*s62 + (s250*s44*s70)/4. + (s250*s44*s71)/4. + (s15*s250*s44*s72)/4. + (s0*s250*s44*s73)/4. + (s16*s250*s44*s74)/4. + (s1*s250*s44*s75)/4. + (s250*s4*s44*s76)/4. + (s14*s250*s44*s77)/4. + (s250*s44*s78*s79)/12. + (s250*s44*s6*s80)/4. + s82*(s112 + (s250*s44*s83)/12.) - (s24*s250*s36*s84)/4. - (s24*s250*s36*s85)/4. + s253*s89 + s81*(s114 + (s250*s44*s9)/4.) + s253*s90 + s263*s94 + s263*s95 + s256*s99);
			},
			simplices.Dim(0),
			ThreadCount()
		);
    }

	}; // SimplicialMeshDetails<2,3,Real,Int>

//------------------------------------------------------------

    template<typename Real, typename Int, typename LInt>
    struct SimplicialMeshDetails<2,4,Real,Int,LInt>
    {
	public:

		using SparseMatrix_T = Sparse::MatrixCSR<Real,Int,LInt>;

	private:

		const Int thread_count = 1;

		static constexpr Int  SIZE = 3;
		static constexpr Real nth  = Inv<Real>( 3 );

	public:

		explicit SimplicialMeshDetails( const Int thread_count_ = 1 ) 
		:
			thread_count(Ramp(thread_count_))
		{}
	
		Int ThreadCount() const
		{
			return thread_count;
		}
	
		std::string ClassName() const
        {
            return std::string("SimplicialMeshDetails<2,4,")+TypeName<Real>+","+TypeName<Int>+">";
        }
		
    void DNearFarToHulls( 
		cref<Tensor2<Real,Int>> V_coords, 
		cref<Tensor2<Int ,Int>> simplices, 
		cref<Tensor2<Real,Int>> P_D_near, 
		cref<Tensor2<Real,Int>> P_D_far, 
		cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
		mref<Tensor3<Real,Int>> buffer
	) const
    {
        TOOLS_PTIMER(timer,ClassName()+"::DNearFarToHulls");

        if( P_D_near.Dim(1) != 23 )
        {
            eprint("in DNearFarToHulls: P_D_near.Dim(1) != 23. Aborting");
        }
        
		cptr<Real> X = V_coords.data();
        cptr<Int>  S = simplices.data();
		cptr<Real> N = P_D_near.data();
		cptr<Real> F = P_D_far.data();
        mptr<Real> B = buffer.data();

		ParallelDo(
			[=]( const Int i )
			{
				Real charge = 0;
				for( Int k = 0; k < SIZE; ++k )
				{
					charge += V_charges[S[SIZE * i + k]];
				}
				charge *= nth;

				const Real s0 = X[4*S[3*i+1]+1];
				const Real s1 = s0*s0;
				const Real s2 = X[4*S[3*i+1]+0];
				const Real s3 = X[4*S[3*i+0]+0];
				const Real s4 = X[4*S[3*i+1]+2];
				const Real s5 = s4*s4;
				const Real s6 = X[4*S[3*i+1]+3];
				const Real s7 = s6*s6;
				const Real s8 = X[4*S[3*i+0]+1];
				const Real s9 = X[4*S[3*i+2]+0];
				const Real s10 = X[4*S[3*i+0]+2];
				const Real s11 = X[4*S[3*i+0]+3];
				const Real s12 = X[4*S[3*i+2]+1];
				const Real s13 = s12*s12;
				const Real s14 = X[4*S[3*i+2]+2];
				const Real s15 = s14*s14;
				const Real s16 = X[4*S[3*i+2]+3];
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
				const Real s151 = std::sqrt(s150);
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
				const Real s238 = -(s199*s210*s218);
				const Real s239 = -(s188*s218*s223);
				const Real s240 = -(s199*s218*s223);
				const Real s241 = -(s188*s215*s218);
				const Real s242 = -(s199*s210*s225*s233);
				const Real s243 = -(s188*s223*s225*s233);
				const Real s244 = -(s210*s218);
				const Real s245 = s188*s218*s236;
				const Real s246 = s234 + s237 + s242 + s243 + s244 + s245;
				const Real s247 = s199*s246;
				const Real s248 = -(s199*s223*s225*s233);
				const Real s249 = -(s188*s215*s225*s233);
				const Real s250 = s199*s218*s236;
				const Real s251 = -(s215*s218);
				const Real s252 = s234 + s237 + s248 + s249 + s250 + s251;
				const Real s253 = s188*s252;
				const Real s254 = s238 + s239 + s240 + s241 + s247 + s253;
				const Real s255 = (charge*s151*s254)/2.;
				const Real s256 = s199*s210*s218;
				const Real s257 = s188*s218*s223;
				const Real s258 = s256 + s257;
				const Real s259 = s199*s258;
				const Real s260 = s199*s218*s223;
				const Real s261 = s188*s215*s218;
				const Real s262 = s260 + s261;
				const Real s263 = s188*s262;
				const Real s264 = s259 + s263;
				const Real s265 = (charge*s152*s186*s264)/4.;
				const Real s266 = s255 + s265;
				const Real s267 = s201*s246;
				const Real s268 = s191*s252;
				const Real s269 = s267 + s268;
				const Real s270 = (charge*s151*s269)/2.;
				const Real s271 = s201*s258;
				const Real s272 = s191*s262;
				const Real s273 = s271 + s272;
				const Real s274 = (charge*s152*s186*s273)/4.;
				const Real s275 = s270 + s274;
				const Real s276 = s203*s246;
				const Real s277 = s194*s252;
				const Real s278 = s276 + s277;
				const Real s279 = (charge*s151*s278)/2.;
				const Real s280 = s203*s258;
				const Real s281 = s194*s262;
				const Real s282 = s280 + s281;
				const Real s283 = (charge*s152*s186*s282)/4.;
				const Real s284 = s279 + s283;
				const Real s285 = s205*s246;
				const Real s286 = s197*s252;
				const Real s287 = s285 + s286;
				const Real s288 = (charge*s151*s287)/2.;
				const Real s289 = s205*s258;
				const Real s290 = s197*s262;
				const Real s291 = s289 + s290;
				const Real s292 = (charge*s152*s186*s291)/4.;
				const Real s293 = s288 + s292;
				const Real s294 = -(s201*s210*s225*s233);
				const Real s295 = -(s191*s223*s225*s233);
				const Real s296 = s191*s218*s236;
				const Real s297 = -2*s188*s201*s218;
				const Real s298 = s294 + s295 + s296 + s297;
				const Real s299 = s201*s298;
				const Real s300 = -(s201*s223*s225*s233);
				const Real s301 = -(s191*s215*s225*s233);
				const Real s302 = -2*s191*s199*s218;
				const Real s303 = s201*s218*s236;
				const Real s304 = s300 + s301 + s302 + s303;
				const Real s305 = s191*s304;
				const Real s306 = s299 + s305;
				const Real s307 = (charge*s151*s306)/2.;
				const Real s308 = s201*s210*s218;
				const Real s309 = s191*s218*s223;
				const Real s310 = s308 + s309;
				const Real s311 = s201*s310;
				const Real s312 = s201*s218*s223;
				const Real s313 = s191*s215*s218;
				const Real s314 = s312 + s313;
				const Real s315 = s191*s314;
				const Real s316 = s311 + s315;
				const Real s317 = (charge*s152*s186*s316)/4.;
				const Real s318 = s307 + s317;
				const Real s319 = s203*s298;
				const Real s320 = s194*s304;
				const Real s321 = s319 + s320;
				const Real s322 = (charge*s151*s321)/2.;
				const Real s323 = s203*s310;
				const Real s324 = s194*s314;
				const Real s325 = s323 + s324;
				const Real s326 = (charge*s152*s186*s325)/4.;
				const Real s327 = s322 + s326;
				const Real s328 = s205*s298;
				const Real s329 = s197*s304;
				const Real s330 = s328 + s329;
				const Real s331 = (charge*s151*s330)/2.;
				const Real s332 = s205*s310;
				const Real s333 = s197*s314;
				const Real s334 = s332 + s333;
				const Real s335 = (charge*s152*s186*s334)/4.;
				const Real s336 = s331 + s335;
				const Real s337 = -(s203*s210*s225*s233);
				const Real s338 = -(s194*s223*s225*s233);
				const Real s339 = s194*s218*s236;
				const Real s340 = -2*s188*s203*s218;
				const Real s341 = s337 + s338 + s339 + s340;
				const Real s342 = s203*s341;
				const Real s343 = -(s203*s223*s225*s233);
				const Real s344 = -(s194*s215*s225*s233);
				const Real s345 = -2*s194*s199*s218;
				const Real s346 = s203*s218*s236;
				const Real s347 = s343 + s344 + s345 + s346;
				const Real s348 = s194*s347;
				const Real s349 = s342 + s348;
				const Real s350 = (charge*s151*s349)/2.;
				const Real s351 = s203*s210*s218;
				const Real s352 = s194*s218*s223;
				const Real s353 = s351 + s352;
				const Real s354 = s203*s353;
				const Real s355 = s203*s218*s223;
				const Real s356 = s194*s215*s218;
				const Real s357 = s355 + s356;
				const Real s358 = s194*s357;
				const Real s359 = s354 + s358;
				const Real s360 = (charge*s152*s186*s359)/4.;
				const Real s361 = s350 + s360;
				const Real s362 = s205*s341;
				const Real s363 = s197*s347;
				const Real s364 = s362 + s363;
				const Real s365 = (charge*s151*s364)/2.;
				const Real s366 = s205*s353;
				const Real s367 = s197*s357;
				const Real s368 = s366 + s367;
				const Real s369 = (charge*s152*s186*s368)/4.;
				const Real s370 = s365 + s369;
				const Real s371 = -(s205*s210*s225*s233);
				const Real s372 = -(s197*s223*s225*s233);
				const Real s373 = s197*s218*s236;
				const Real s374 = -2*s188*s205*s218;
				const Real s375 = s371 + s372 + s373 + s374;
				const Real s376 = s205*s375;
				const Real s377 = -(s205*s223*s225*s233);
				const Real s378 = -(s197*s215*s225*s233);
				const Real s379 = -2*s197*s199*s218;
				const Real s380 = s205*s218*s236;
				const Real s381 = s377 + s378 + s379 + s380;
				const Real s382 = s197*s381;
				const Real s383 = s376 + s382;
				const Real s384 = (charge*s151*s383)/2.;
				const Real s385 = s205*s210*s218;
				const Real s386 = s197*s218*s223;
				const Real s387 = s385 + s386;
				const Real s388 = s205*s387;
				const Real s389 = s205*s218*s223;
				const Real s390 = s197*s215*s218;
				const Real s391 = s389 + s390;
				const Real s392 = s197*s391;
				const Real s393 = s388 + s392;
				const Real s394 = (charge*s152*s186*s393)/4.;
				const Real s395 = s384 + s394;
				const Real s396 = F[15*i+0];
				const Real s397 = N[23*i+0];
				const Real s398 = 2*s19*s8;
				const Real s399 = -2*s0*s2*s3;
				const Real s400 = -2*s0*s10*s4;
				const Real s401 = 2*s5*s8;
				const Real s402 = -2*s0*s11*s6;
				const Real s403 = 2*s7*s8;
				const Real s404 = -4*s2*s8*s9;
				const Real s405 = 2*s0*s3*s9;
				const Real s406 = 2*s0*s2*s9;
				const Real s407 = 2*s53*s8;
				const Real s408 = -2*s0*s53;
				const Real s409 = 2*s12*s2*s3;
				const Real s410 = -2*s12*s19;
				const Real s411 = 2*s10*s12*s4;
				const Real s412 = -2*s12*s5;
				const Real s413 = 2*s11*s12*s6;
				const Real s414 = -2*s12*s7;
				const Real s415 = -2*s12*s3*s9;
				const Real s416 = 2*s12*s2*s9;
				const Real s417 = 2*s0*s10*s14;
				const Real s418 = -4*s14*s4*s8;
				const Real s419 = 2*s0*s14*s4;
				const Real s420 = -2*s10*s12*s14;
				const Real s421 = 2*s12*s14*s4;
				const Real s422 = 2*s15*s8;
				const Real s423 = -2*s0*s15;
				const Real s424 = 2*s0*s11*s16;
				const Real s425 = -4*s16*s6*s8;
				const Real s426 = 2*s0*s16*s6;
				const Real s427 = -2*s11*s12*s16;
				const Real s428 = 2*s12*s16*s6;
				const Real s429 = 2*s17*s8;
				const Real s430 = -2*s0*s17;
				const Real s431 = s398 + s399 + s400 + s401 + s402 + s403 + s404 + s405 + s406 + s407 + s408 + s409 + s410 + s411 + s412 + s413 + s414 + s415 + s416 + s417 + s418 + s419 + s420 + s421 + s422 + s423 + s424 + s425 + s426 + s427 + s428 + s429 + s430;
				const Real s432 = N[23*i+1];
				const Real s433 = N[23*i+3];
				const Real s434 = N[23*i+4];
				const Real s435 = N[23*i+5];
				const Real s436 = N[23*i+6];
				const Real s437 = N[23*i+7];
				const Real s438 = N[23*i+8];
				const Real s439 = N[23*i+9];
				const Real s440 = F[15*i+1];
				const Real s441 = s2 + s3 + s9;
				const Real s442 = N[23*i+10];
				const Real s443 = N[23*i+11];
				const Real s444 = F[15*i+3];
				const Real s445 = s10 + s14 + s4;
				const Real s446 = N[23*i+12];
				const Real s447 = F[15*i+4];
				const Real s448 = s11 + s16 + s6;
				const Real s449 = F[15*i+2];
				const Real s450 = s0 + s12 + s8;
				const Real s451 = (charge*s151)/6.;
				const Real s452 = N[23*i+2];
				const Real s453 = (charge*s151)/2.;
				const Real s454 = F[15*i+5];
				const Real s455 = -2*s201*s210;
				const Real s456 = 2*s8;
				const Real s457 = -s0;
				const Real s458 = -s12;
				const Real s459 = s456 + s457 + s458;
				const Real s460 = -2*s207*s459;
				const Real s461 = -2*s191*s215;
				const Real s462 = s455 + s460 + s461;
				const Real s463 = -2*s8;
				const Real s464 = s0 + s12 + s463;
				const Real s465 = N[23*i+13];
				const Real s466 = -(s199*s210*s225*s462);
				const Real s467 = -(s188*s223*s225*s462);
				const Real s468 = s188*s218*s464;
				const Real s469 = s302 + s466 + s467 + s468;
				const Real s470 = s199*s469;
				const Real s471 = -(s199*s223*s225*s462);
				const Real s472 = -(s188*s215*s225*s462);
				const Real s473 = s199*s218*s464;
				const Real s474 = s297 + s471 + s472 + s473;
				const Real s475 = s188*s474;
				const Real s476 = s470 + s475;
				const Real s477 = (charge*s151*s476)/2.;
				const Real s478 = (charge*s152*s264*s431)/4.;
				const Real s479 = s477 + s478;
				const Real s480 = F[15*i+6];
				const Real s481 = N[23*i+14];
				const Real s482 = s201*s469;
				const Real s483 = s191*s474;
				const Real s484 = s238 + s239 + s240 + s241 + s482 + s483;
				const Real s485 = (charge*s151*s484)/2.;
				const Real s486 = (charge*s152*s273*s431)/4.;
				const Real s487 = s485 + s486;
				const Real s488 = F[15*i+7];
				const Real s489 = N[23*i+15];
				const Real s490 = s203*s469;
				const Real s491 = s194*s474;
				const Real s492 = s490 + s491;
				const Real s493 = (charge*s151*s492)/2.;
				const Real s494 = (charge*s152*s282*s431)/4.;
				const Real s495 = s493 + s494;
				const Real s496 = F[15*i+8];
				const Real s497 = N[23*i+16];
				const Real s498 = s205*s469;
				const Real s499 = s197*s474;
				const Real s500 = s498 + s499;
				const Real s501 = (charge*s151*s500)/2.;
				const Real s502 = (charge*s152*s291*s431)/4.;
				const Real s503 = s501 + s502;
				const Real s504 = F[15*i+9];
				const Real s505 = -2*s191*s201*s218;
				const Real s506 = N[23*i+17];
				const Real s507 = -(s201*s210*s218);
				const Real s508 = -(s191*s218*s223);
				const Real s509 = -(s201*s218*s223);
				const Real s510 = -(s191*s215*s218);
				const Real s511 = -(s201*s210*s225*s462);
				const Real s512 = -(s191*s223*s225*s462);
				const Real s513 = s191*s218*s464;
				const Real s514 = s237 + s244 + s505 + s511 + s512 + s513;
				const Real s515 = s201*s514;
				const Real s516 = -(s201*s223*s225*s462);
				const Real s517 = -(s191*s215*s225*s462);
				const Real s518 = s201*s218*s464;
				const Real s519 = s237 + s251 + s505 + s516 + s517 + s518;
				const Real s520 = s191*s519;
				const Real s521 = s507 + s508 + s509 + s510 + s515 + s520;
				const Real s522 = (charge*s151*s521)/2.;
				const Real s523 = (charge*s152*s316*s431)/4.;
				const Real s524 = s522 + s523;
				const Real s525 = F[15*i+10];
				const Real s526 = N[23*i+18];
				const Real s527 = s203*s514;
				const Real s528 = s194*s519;
				const Real s529 = s527 + s528;
				const Real s530 = (charge*s151*s529)/2.;
				const Real s531 = (charge*s152*s325*s431)/4.;
				const Real s532 = s530 + s531;
				const Real s533 = F[15*i+11];
				const Real s534 = N[23*i+19];
				const Real s535 = s205*s514;
				const Real s536 = s197*s519;
				const Real s537 = s535 + s536;
				const Real s538 = (charge*s151*s537)/2.;
				const Real s539 = (charge*s152*s334*s431)/4.;
				const Real s540 = s538 + s539;
				const Real s541 = F[15*i+12];
				const Real s542 = N[23*i+20];
				const Real s543 = -(s203*s210*s225*s462);
				const Real s544 = -(s194*s223*s225*s462);
				const Real s545 = s194*s218*s464;
				const Real s546 = -2*s191*s203*s218;
				const Real s547 = s543 + s544 + s545 + s546;
				const Real s548 = s203*s547;
				const Real s549 = -(s203*s223*s225*s462);
				const Real s550 = -(s194*s215*s225*s462);
				const Real s551 = -2*s194*s201*s218;
				const Real s552 = s203*s218*s464;
				const Real s553 = s549 + s550 + s551 + s552;
				const Real s554 = s194*s553;
				const Real s555 = s548 + s554;
				const Real s556 = (charge*s151*s555)/2.;
				const Real s557 = (charge*s152*s359*s431)/4.;
				const Real s558 = s556 + s557;
				const Real s559 = F[15*i+13];
				const Real s560 = N[23*i+21];
				const Real s561 = s205*s547;
				const Real s562 = s197*s553;
				const Real s563 = s561 + s562;
				const Real s564 = (charge*s151*s563)/2.;
				const Real s565 = (charge*s152*s368*s431)/4.;
				const Real s566 = s564 + s565;
				const Real s567 = F[15*i+14];
				const Real s568 = N[23*i+22];
				const Real s569 = -(s205*s210*s225*s462);
				const Real s570 = -(s197*s223*s225*s462);
				const Real s571 = s197*s218*s464;
				const Real s572 = -2*s191*s205*s218;
				const Real s573 = s569 + s570 + s571 + s572;
				const Real s574 = s205*s573;
				const Real s575 = -(s205*s223*s225*s462);
				const Real s576 = -(s197*s215*s225*s462);
				const Real s577 = -2*s197*s201*s218;
				const Real s578 = s205*s218*s464;
				const Real s579 = s575 + s576 + s577 + s578;
				const Real s580 = s197*s579;
				const Real s581 = s574 + s580;
				const Real s582 = (charge*s151*s581)/2.;
				const Real s583 = (charge*s152*s393*s431)/4.;
				const Real s584 = s582 + s583;
				const Real s585 = 2*s10*s19;
				const Real s586 = 2*s1*s10;
				const Real s587 = -2*s2*s3*s4;
				const Real s588 = -2*s0*s4*s8;
				const Real s589 = -2*s11*s4*s6;
				const Real s590 = 2*s10*s7;
				const Real s591 = -4*s10*s2*s9;
				const Real s592 = 2*s3*s4*s9;
				const Real s593 = 2*s2*s4*s9;
				const Real s594 = 2*s10*s53;
				const Real s595 = -2*s4*s53;
				const Real s596 = -4*s0*s10*s12;
				const Real s597 = 2*s12*s4*s8;
				const Real s598 = 2*s0*s12*s4;
				const Real s599 = 2*s10*s13;
				const Real s600 = -2*s13*s4;
				const Real s601 = 2*s14*s2*s3;
				const Real s602 = -2*s14*s19;
				const Real s603 = 2*s0*s14*s8;
				const Real s604 = -2*s1*s14;
				const Real s605 = 2*s11*s14*s6;
				const Real s606 = -2*s14*s7;
				const Real s607 = -2*s14*s3*s9;
				const Real s608 = 2*s14*s2*s9;
				const Real s609 = -2*s12*s14*s8;
				const Real s610 = 2*s0*s12*s14;
				const Real s611 = 2*s11*s16*s4;
				const Real s612 = -4*s10*s16*s6;
				const Real s613 = 2*s16*s4*s6;
				const Real s614 = -2*s11*s14*s16;
				const Real s615 = 2*s14*s16*s6;
				const Real s616 = 2*s10*s17;
				const Real s617 = -2*s17*s4;
				const Real s618 = s585 + s586 + s587 + s588 + s589 + s590 + s591 + s592 + s593 + s594 + s595 + s596 + s597 + s598 + s599 + s600 + s601 + s602 + s603 + s604 + s605 + s606 + s607 + s608 + s609 + s610 + s611 + s612 + s613 + s614 + s615 + s616 + s617;
				const Real s619 = -2*s203*s210;
				const Real s620 = 2*s10;
				const Real s621 = -s4;
				const Real s622 = -s14;
				const Real s623 = s620 + s621 + s622;
				const Real s624 = -2*s207*s623;
				const Real s625 = -2*s194*s215;
				const Real s626 = s619 + s624 + s625;
				const Real s627 = -2*s10;
				const Real s628 = s14 + s4 + s627;
				const Real s629 = -(s199*s210*s225*s626);
				const Real s630 = -(s188*s223*s225*s626);
				const Real s631 = s188*s218*s628;
				const Real s632 = s345 + s629 + s630 + s631;
				const Real s633 = s199*s632;
				const Real s634 = -(s199*s223*s225*s626);
				const Real s635 = -(s188*s215*s225*s626);
				const Real s636 = s199*s218*s628;
				const Real s637 = s340 + s634 + s635 + s636;
				const Real s638 = s188*s637;
				const Real s639 = s633 + s638;
				const Real s640 = (charge*s151*s639)/2.;
				const Real s641 = (charge*s152*s264*s618)/4.;
				const Real s642 = s640 + s641;
				const Real s643 = s201*s632;
				const Real s644 = s191*s637;
				const Real s645 = s643 + s644;
				const Real s646 = (charge*s151*s645)/2.;
				const Real s647 = (charge*s152*s273*s618)/4.;
				const Real s648 = s646 + s647;
				const Real s649 = s203*s632;
				const Real s650 = s194*s637;
				const Real s651 = s238 + s239 + s240 + s241 + s649 + s650;
				const Real s652 = (charge*s151*s651)/2.;
				const Real s653 = (charge*s152*s282*s618)/4.;
				const Real s654 = s652 + s653;
				const Real s655 = s205*s632;
				const Real s656 = s197*s637;
				const Real s657 = s655 + s656;
				const Real s658 = (charge*s151*s657)/2.;
				const Real s659 = (charge*s152*s291*s618)/4.;
				const Real s660 = s658 + s659;
				const Real s661 = -(s201*s210*s225*s626);
				const Real s662 = -(s191*s223*s225*s626);
				const Real s663 = s191*s218*s628;
				const Real s664 = s551 + s661 + s662 + s663;
				const Real s665 = s201*s664;
				const Real s666 = -(s201*s223*s225*s626);
				const Real s667 = -(s191*s215*s225*s626);
				const Real s668 = s201*s218*s628;
				const Real s669 = s546 + s666 + s667 + s668;
				const Real s670 = s191*s669;
				const Real s671 = s665 + s670;
				const Real s672 = (charge*s151*s671)/2.;
				const Real s673 = (charge*s152*s316*s618)/4.;
				const Real s674 = s672 + s673;
				const Real s675 = s203*s664;
				const Real s676 = s194*s669;
				const Real s677 = s507 + s508 + s509 + s510 + s675 + s676;
				const Real s678 = (charge*s151*s677)/2.;
				const Real s679 = (charge*s152*s325*s618)/4.;
				const Real s680 = s678 + s679;
				const Real s681 = s205*s664;
				const Real s682 = s197*s669;
				const Real s683 = s681 + s682;
				const Real s684 = (charge*s151*s683)/2.;
				const Real s685 = (charge*s152*s334*s618)/4.;
				const Real s686 = s684 + s685;
				const Real s687 = -2*s194*s203*s218;
				const Real s688 = -(s203*s210*s218);
				const Real s689 = -(s194*s218*s223);
				const Real s690 = -(s203*s218*s223);
				const Real s691 = -(s194*s215*s218);
				const Real s692 = -(s203*s210*s225*s626);
				const Real s693 = -(s194*s223*s225*s626);
				const Real s694 = s194*s218*s628;
				const Real s695 = s237 + s244 + s687 + s692 + s693 + s694;
				const Real s696 = s203*s695;
				const Real s697 = -(s203*s223*s225*s626);
				const Real s698 = -(s194*s215*s225*s626);
				const Real s699 = s203*s218*s628;
				const Real s700 = s237 + s251 + s687 + s697 + s698 + s699;
				const Real s701 = s194*s700;
				const Real s702 = s688 + s689 + s690 + s691 + s696 + s701;
				const Real s703 = (charge*s151*s702)/2.;
				const Real s704 = (charge*s152*s359*s618)/4.;
				const Real s705 = s703 + s704;
				const Real s706 = s205*s695;
				const Real s707 = s197*s700;
				const Real s708 = s706 + s707;
				const Real s709 = (charge*s151*s708)/2.;
				const Real s710 = (charge*s152*s368*s618)/4.;
				const Real s711 = s709 + s710;
				const Real s712 = -(s205*s210*s225*s626);
				const Real s713 = -(s197*s223*s225*s626);
				const Real s714 = s197*s218*s628;
				const Real s715 = -2*s194*s205*s218;
				const Real s716 = s712 + s713 + s714 + s715;
				const Real s717 = s205*s716;
				const Real s718 = -(s205*s223*s225*s626);
				const Real s719 = -(s197*s215*s225*s626);
				const Real s720 = -2*s197*s203*s218;
				const Real s721 = s205*s218*s628;
				const Real s722 = s718 + s719 + s720 + s721;
				const Real s723 = s197*s722;
				const Real s724 = s717 + s723;
				const Real s725 = (charge*s151*s724)/2.;
				const Real s726 = (charge*s152*s393*s618)/4.;
				const Real s727 = s725 + s726;
				const Real s728 = 2*s11*s19;
				const Real s729 = 2*s1*s11;
				const Real s730 = 2*s11*s5;
				const Real s731 = -2*s2*s3*s6;
				const Real s732 = -2*s0*s6*s8;
				const Real s733 = -2*s10*s4*s6;
				const Real s734 = -4*s11*s2*s9;
				const Real s735 = 2*s3*s6*s9;
				const Real s736 = 2*s2*s6*s9;
				const Real s737 = 2*s11*s53;
				const Real s738 = -2*s53*s6;
				const Real s739 = -4*s0*s11*s12;
				const Real s740 = 2*s12*s6*s8;
				const Real s741 = 2*s0*s12*s6;
				const Real s742 = 2*s11*s13;
				const Real s743 = -2*s13*s6;
				const Real s744 = -4*s11*s14*s4;
				const Real s745 = 2*s10*s14*s6;
				const Real s746 = 2*s14*s4*s6;
				const Real s747 = 2*s11*s15;
				const Real s748 = -2*s15*s6;
				const Real s749 = 2*s16*s2*s3;
				const Real s750 = -2*s16*s19;
				const Real s751 = 2*s0*s16*s8;
				const Real s752 = -2*s1*s16;
				const Real s753 = 2*s10*s16*s4;
				const Real s754 = -2*s16*s5;
				const Real s755 = -2*s16*s3*s9;
				const Real s756 = 2*s16*s2*s9;
				const Real s757 = -2*s12*s16*s8;
				const Real s758 = 2*s0*s12*s16;
				const Real s759 = -2*s10*s14*s16;
				const Real s760 = 2*s14*s16*s4;
				const Real s761 = s728 + s729 + s730 + s731 + s732 + s733 + s734 + s735 + s736 + s737 + s738 + s739 + s740 + s741 + s742 + s743 + s744 + s745 + s746 + s747 + s748 + s749 + s750 + s751 + s752 + s753 + s754 + s755 + s756 + s757 + s758 + s759 + s760;
				const Real s762 = -2*s205*s210;
				const Real s763 = 2*s11;
				const Real s764 = -s6;
				const Real s765 = -s16;
				const Real s766 = s763 + s764 + s765;
				const Real s767 = -2*s207*s766;
				const Real s768 = -2*s197*s215;
				const Real s769 = s762 + s767 + s768;
				const Real s770 = -2*s11;
				const Real s771 = s16 + s6 + s770;
				const Real s772 = -(s199*s210*s225*s769);
				const Real s773 = -(s188*s223*s225*s769);
				const Real s774 = s188*s218*s771;
				const Real s775 = s379 + s772 + s773 + s774;
				const Real s776 = s199*s775;
				const Real s777 = -(s199*s223*s225*s769);
				const Real s778 = -(s188*s215*s225*s769);
				const Real s779 = s199*s218*s771;
				const Real s780 = s374 + s777 + s778 + s779;
				const Real s781 = s188*s780;
				const Real s782 = s776 + s781;
				const Real s783 = (charge*s151*s782)/2.;
				const Real s784 = (charge*s152*s264*s761)/4.;
				const Real s785 = s783 + s784;
				const Real s786 = s201*s775;
				const Real s787 = s191*s780;
				const Real s788 = s786 + s787;
				const Real s789 = (charge*s151*s788)/2.;
				const Real s790 = (charge*s152*s273*s761)/4.;
				const Real s791 = s789 + s790;
				const Real s792 = s203*s775;
				const Real s793 = s194*s780;
				const Real s794 = s792 + s793;
				const Real s795 = (charge*s151*s794)/2.;
				const Real s796 = (charge*s152*s282*s761)/4.;
				const Real s797 = s795 + s796;
				const Real s798 = s205*s775;
				const Real s799 = s197*s780;
				const Real s800 = s238 + s239 + s240 + s241 + s798 + s799;
				const Real s801 = (charge*s151*s800)/2.;
				const Real s802 = (charge*s152*s291*s761)/4.;
				const Real s803 = s801 + s802;
				const Real s804 = -(s201*s210*s225*s769);
				const Real s805 = -(s191*s223*s225*s769);
				const Real s806 = s191*s218*s771;
				const Real s807 = s577 + s804 + s805 + s806;
				const Real s808 = s201*s807;
				const Real s809 = -(s201*s223*s225*s769);
				const Real s810 = -(s191*s215*s225*s769);
				const Real s811 = s201*s218*s771;
				const Real s812 = s572 + s809 + s810 + s811;
				const Real s813 = s191*s812;
				const Real s814 = s808 + s813;
				const Real s815 = (charge*s151*s814)/2.;
				const Real s816 = (charge*s152*s316*s761)/4.;
				const Real s817 = s815 + s816;
				const Real s818 = s203*s807;
				const Real s819 = s194*s812;
				const Real s820 = s818 + s819;
				const Real s821 = (charge*s151*s820)/2.;
				const Real s822 = (charge*s152*s325*s761)/4.;
				const Real s823 = s821 + s822;
				const Real s824 = s205*s807;
				const Real s825 = s197*s812;
				const Real s826 = s507 + s508 + s509 + s510 + s824 + s825;
				const Real s827 = (charge*s151*s826)/2.;
				const Real s828 = (charge*s152*s334*s761)/4.;
				const Real s829 = s827 + s828;
				const Real s830 = -(s203*s210*s225*s769);
				const Real s831 = -(s194*s223*s225*s769);
				const Real s832 = s194*s218*s771;
				const Real s833 = s720 + s830 + s831 + s832;
				const Real s834 = s203*s833;
				const Real s835 = -(s203*s223*s225*s769);
				const Real s836 = -(s194*s215*s225*s769);
				const Real s837 = s203*s218*s771;
				const Real s838 = s715 + s835 + s836 + s837;
				const Real s839 = s194*s838;
				const Real s840 = s834 + s839;
				const Real s841 = (charge*s151*s840)/2.;
				const Real s842 = (charge*s152*s359*s761)/4.;
				const Real s843 = s841 + s842;
				const Real s844 = s205*s833;
				const Real s845 = s197*s838;
				const Real s846 = s688 + s689 + s690 + s691 + s844 + s845;
				const Real s847 = (charge*s151*s846)/2.;
				const Real s848 = (charge*s152*s368*s761)/4.;
				const Real s849 = s847 + s848;
				const Real s850 = -2*s197*s205*s218;
				const Real s851 = -(s205*s210*s218);
				const Real s852 = -(s197*s218*s223);
				const Real s853 = -(s205*s218*s223);
				const Real s854 = -(s197*s215*s218);
				const Real s855 = -(s205*s210*s225*s769);
				const Real s856 = -(s197*s223*s225*s769);
				const Real s857 = s197*s218*s771;
				const Real s858 = s237 + s244 + s850 + s855 + s856 + s857;
				const Real s859 = s205*s858;
				const Real s860 = -(s205*s223*s225*s769);
				const Real s861 = -(s197*s215*s225*s769);
				const Real s862 = s205*s218*s771;
				const Real s863 = s237 + s251 + s850 + s860 + s861 + s862;
				const Real s864 = s197*s863;
				const Real s865 = s851 + s852 + s853 + s854 + s859 + s864;
				const Real s866 = (charge*s151*s865)/2.;
				const Real s867 = (charge*s152*s393*s761)/4.;
				const Real s868 = s866 + s867;
				const Real s869 = 2*s18*s2;
				const Real s870 = 2*s2*s21;
				const Real s871 = 2*s2*s23;
				const Real s872 = -2*s0*s3*s8;
				const Real s873 = -2*s10*s3*s4;
				const Real s874 = -2*s11*s3*s6;
				const Real s875 = -2*s18*s9;
				const Real s876 = -2*s21*s9;
				const Real s877 = -2*s23*s9;
				const Real s878 = 2*s12*s3*s8;
				const Real s879 = -4*s12*s2*s8;
				const Real s880 = 2*s0*s12*s3;
				const Real s881 = 2*s12*s8*s9;
				const Real s882 = -2*s0*s12*s9;
				const Real s883 = -2*s13*s3;
				const Real s884 = 2*s13*s2;
				const Real s885 = 2*s10*s14*s3;
				const Real s886 = -4*s10*s14*s2;
				const Real s887 = 2*s14*s3*s4;
				const Real s888 = 2*s10*s14*s9;
				const Real s889 = -2*s14*s4*s9;
				const Real s890 = -2*s15*s3;
				const Real s891 = 2*s15*s2;
				const Real s892 = 2*s11*s16*s3;
				const Real s893 = -4*s11*s16*s2;
				const Real s894 = 2*s16*s3*s6;
				const Real s895 = 2*s11*s16*s9;
				const Real s896 = -2*s16*s6*s9;
				const Real s897 = -2*s17*s3;
				const Real s898 = 2*s17*s2;
				const Real s899 = s159 + s161 + s163 + s869 + s870 + s871 + s872 + s873 + s874 + s875 + s876 + s877 + s878 + s879 + s880 + s881 + s882 + s883 + s884 + s885 + s886 + s887 + s888 + s889 + s890 + s891 + s892 + s893 + s894 + s895 + s896 + s897 + s898;
				const Real s900 = -2*s199*s207;
				const Real s901 = 2*s188*s215;
				const Real s902 = s900 + s901;
				const Real s903 = -(s199*s210*s225*s902);
				const Real s904 = -(s188*s223*s225*s902);
				const Real s905 = s188*s199*s218;
				const Real s906 = s218*s223;
				const Real s907 = s903 + s904 + s905 + s906;
				const Real s908 = s199*s907;
				const Real s909 = -(s199*s223*s225*s902);
				const Real s910 = -(s188*s215*s225*s902);
				const Real s911 = -(s211*s218);
				const Real s912 = s215*s218;
				const Real s913 = s909 + s910 + s911 + s912;
				const Real s914 = s188*s913;
				const Real s915 = s260 + s261 + s908 + s914;
				const Real s916 = (charge*s151*s915)/2.;
				const Real s917 = (charge*s152*s264*s899)/4.;
				const Real s918 = s916 + s917;
				const Real s919 = s201*s907;
				const Real s920 = s191*s913;
				const Real s921 = s919 + s920;
				const Real s922 = (charge*s151*s921)/2.;
				const Real s923 = (charge*s152*s273*s899)/4.;
				const Real s924 = s922 + s923;
				const Real s925 = s203*s907;
				const Real s926 = s194*s913;
				const Real s927 = s925 + s926;
				const Real s928 = (charge*s151*s927)/2.;
				const Real s929 = (charge*s152*s282*s899)/4.;
				const Real s930 = s928 + s929;
				const Real s931 = s205*s907;
				const Real s932 = s197*s913;
				const Real s933 = s931 + s932;
				const Real s934 = (charge*s151*s933)/2.;
				const Real s935 = (charge*s152*s291*s899)/4.;
				const Real s936 = s934 + s935;
				const Real s937 = -(s201*s210*s225*s902);
				const Real s938 = -(s191*s223*s225*s902);
				const Real s939 = -(s191*s199*s218);
				const Real s940 = 2*s188*s201*s218;
				const Real s941 = s937 + s938 + s939 + s940;
				const Real s942 = s201*s941;
				const Real s943 = -(s201*s223*s225*s902);
				const Real s944 = -(s191*s215*s225*s902);
				const Real s945 = -(s199*s201*s218);
				const Real s946 = s943 + s944 + s945;
				const Real s947 = s191*s946;
				const Real s948 = s942 + s947;
				const Real s949 = (charge*s151*s948)/2.;
				const Real s950 = (charge*s152*s316*s899)/4.;
				const Real s951 = s949 + s950;
				const Real s952 = s203*s941;
				const Real s953 = s194*s946;
				const Real s954 = s952 + s953;
				const Real s955 = (charge*s151*s954)/2.;
				const Real s956 = (charge*s152*s325*s899)/4.;
				const Real s957 = s955 + s956;
				const Real s958 = s205*s941;
				const Real s959 = s197*s946;
				const Real s960 = s958 + s959;
				const Real s961 = (charge*s151*s960)/2.;
				const Real s962 = (charge*s152*s334*s899)/4.;
				const Real s963 = s961 + s962;
				const Real s964 = -(s203*s210*s225*s902);
				const Real s965 = -(s194*s223*s225*s902);
				const Real s966 = -(s194*s199*s218);
				const Real s967 = 2*s188*s203*s218;
				const Real s968 = s964 + s965 + s966 + s967;
				const Real s969 = s203*s968;
				const Real s970 = -(s203*s223*s225*s902);
				const Real s971 = -(s194*s215*s225*s902);
				const Real s972 = -(s199*s203*s218);
				const Real s973 = s970 + s971 + s972;
				const Real s974 = s194*s973;
				const Real s975 = s969 + s974;
				const Real s976 = (charge*s151*s975)/2.;
				const Real s977 = (charge*s152*s359*s899)/4.;
				const Real s978 = s976 + s977;
				const Real s979 = s205*s968;
				const Real s980 = s197*s973;
				const Real s981 = s979 + s980;
				const Real s982 = (charge*s151*s981)/2.;
				const Real s983 = (charge*s152*s368*s899)/4.;
				const Real s984 = s982 + s983;
				const Real s985 = -(s205*s210*s225*s902);
				const Real s986 = -(s197*s223*s225*s902);
				const Real s987 = -(s197*s199*s218);
				const Real s988 = 2*s188*s205*s218;
				const Real s989 = s985 + s986 + s987 + s988;
				const Real s990 = s205*s989;
				const Real s991 = -(s205*s223*s225*s902);
				const Real s992 = -(s197*s215*s225*s902);
				const Real s993 = -(s199*s205*s218);
				const Real s994 = s991 + s992 + s993;
				const Real s995 = s197*s994;
				const Real s996 = s990 + s995;
				const Real s997 = (charge*s151*s996)/2.;
				const Real s998 = (charge*s152*s393*s899)/4.;
				const Real s999 = s997 + s998;
				const Real s1000 = -2*s2*s3*s8;
				const Real s1001 = 2*s0*s26;
				const Real s1002 = 2*s0*s21;
				const Real s1003 = 2*s0*s23;
				const Real s1004 = -2*s10*s4*s8;
				const Real s1005 = -2*s11*s6*s8;
				const Real s1006 = 2*s3*s8*s9;
				const Real s1007 = 2*s2*s8*s9;
				const Real s1008 = -4*s0*s3*s9;
				const Real s1009 = -2*s53*s8;
				const Real s1010 = 2*s0*s53;
				const Real s1011 = -2*s12*s26;
				const Real s1012 = -2*s12*s21;
				const Real s1013 = -2*s12*s23;
				const Real s1014 = 2*s12*s3*s9;
				const Real s1015 = -2*s12*s2*s9;
				const Real s1016 = 2*s10*s14*s8;
				const Real s1017 = -4*s0*s10*s14;
				const Real s1018 = 2*s14*s4*s8;
				const Real s1019 = 2*s10*s12*s14;
				const Real s1020 = -2*s12*s14*s4;
				const Real s1021 = -2*s15*s8;
				const Real s1022 = 2*s0*s15;
				const Real s1023 = 2*s11*s16*s8;
				const Real s1024 = -4*s0*s11*s16;
				const Real s1025 = 2*s16*s6*s8;
				const Real s1026 = 2*s11*s12*s16;
				const Real s1027 = -2*s12*s16*s6;
				const Real s1028 = -2*s17*s8;
				const Real s1029 = 2*s0*s17;
				const Real s1030 = s1000 + s1001 + s1002 + s1003 + s1004 + s1005 + s1006 + s1007 + s1008 + s1009 + s1010 + s1011 + s1012 + s1013 + s1014 + s1015 + s1016 + s1017 + s1018 + s1019 + s1020 + s1021 + s1022 + s1023 + s1024 + s1025 + s1026 + s1027 + s1028 + s1029 + s409 + s411 + s413;
				const Real s1031 = -2*s201*s207;
				const Real s1032 = 2*s191*s215;
				const Real s1033 = s1031 + s1032;
				const Real s1034 = -(s1033*s199*s210*s225);
				const Real s1035 = -(s1033*s188*s223*s225);
				const Real s1036 = 2*s191*s199*s218;
				const Real s1037 = -(s188*s201*s218);
				const Real s1038 = s1034 + s1035 + s1036 + s1037;
				const Real s1039 = s1038*s199;
				const Real s1040 = -(s1033*s199*s223*s225);
				const Real s1041 = -(s1033*s188*s215*s225);
				const Real s1042 = s1040 + s1041 + s945;
				const Real s1043 = s1042*s188;
				const Real s1044 = s1039 + s1043;
				const Real s1045 = (charge*s1044*s151)/2.;
				const Real s1046 = (charge*s1030*s152*s264)/4.;
				const Real s1047 = s1045 + s1046;
				const Real s1048 = s1038*s201;
				const Real s1049 = s1042*s191;
				const Real s1050 = s1048 + s1049 + s260 + s261;
				const Real s1051 = (charge*s1050*s151)/2.;
				const Real s1052 = (charge*s1030*s152*s273)/4.;
				const Real s1053 = s1051 + s1052;
				const Real s1054 = s1038*s203;
				const Real s1055 = s1042*s194;
				const Real s1056 = s1054 + s1055;
				const Real s1057 = (charge*s1056*s151)/2.;
				const Real s1058 = (charge*s1030*s152*s282)/4.;
				const Real s1059 = s1057 + s1058;
				const Real s1060 = s1038*s205;
				const Real s1061 = s1042*s197;
				const Real s1062 = s1060 + s1061;
				const Real s1063 = (charge*s1062*s151)/2.;
				const Real s1064 = (charge*s1030*s152*s291)/4.;
				const Real s1065 = s1063 + s1064;
				const Real s1066 = -(s1033*s201*s210*s225);
				const Real s1067 = -(s1033*s191*s223*s225);
				const Real s1068 = s191*s201*s218;
				const Real s1069 = s1066 + s1067 + s1068 + s906;
				const Real s1070 = s1069*s201;
				const Real s1071 = -(s1033*s201*s223*s225);
				const Real s1072 = -(s1033*s191*s215*s225);
				const Real s1073 = -(s212*s218);
				const Real s1074 = s1071 + s1072 + s1073 + s912;
				const Real s1075 = s1074*s191;
				const Real s1076 = s1070 + s1075 + s312 + s313;
				const Real s1077 = (charge*s1076*s151)/2.;
				const Real s1078 = (charge*s1030*s152*s316)/4.;
				const Real s1079 = s1077 + s1078;
				const Real s1080 = s1069*s203;
				const Real s1081 = s1074*s194;
				const Real s1082 = s1080 + s1081;
				const Real s1083 = (charge*s1082*s151)/2.;
				const Real s1084 = (charge*s1030*s152*s325)/4.;
				const Real s1085 = s1083 + s1084;
				const Real s1086 = s1069*s205;
				const Real s1087 = s1074*s197;
				const Real s1088 = s1086 + s1087;
				const Real s1089 = (charge*s1088*s151)/2.;
				const Real s1090 = (charge*s1030*s152*s334)/4.;
				const Real s1091 = s1089 + s1090;
				const Real s1092 = -(s1033*s203*s210*s225);
				const Real s1093 = -(s1033*s194*s223*s225);
				const Real s1094 = -(s194*s201*s218);
				const Real s1095 = 2*s191*s203*s218;
				const Real s1096 = s1092 + s1093 + s1094 + s1095;
				const Real s1097 = s1096*s203;
				const Real s1098 = -(s1033*s203*s223*s225);
				const Real s1099 = -(s1033*s194*s215*s225);
				const Real s1100 = -(s201*s203*s218);
				const Real s1101 = s1098 + s1099 + s1100;
				const Real s1102 = s1101*s194;
				const Real s1103 = s1097 + s1102;
				const Real s1104 = (charge*s1103*s151)/2.;
				const Real s1105 = (charge*s1030*s152*s359)/4.;
				const Real s1106 = s1104 + s1105;
				const Real s1107 = s1096*s205;
				const Real s1108 = s1101*s197;
				const Real s1109 = s1107 + s1108;
				const Real s1110 = (charge*s1109*s151)/2.;
				const Real s1111 = (charge*s1030*s152*s368)/4.;
				const Real s1112 = s1110 + s1111;
				const Real s1113 = -(s1033*s205*s210*s225);
				const Real s1114 = -(s1033*s197*s223*s225);
				const Real s1115 = -(s197*s201*s218);
				const Real s1116 = 2*s191*s205*s218;
				const Real s1117 = s1113 + s1114 + s1115 + s1116;
				const Real s1118 = s1117*s205;
				const Real s1119 = -(s1033*s205*s223*s225);
				const Real s1120 = -(s1033*s197*s215*s225);
				const Real s1121 = -(s201*s205*s218);
				const Real s1122 = s1119 + s1120 + s1121;
				const Real s1123 = s1122*s197;
				const Real s1124 = s1118 + s1123;
				const Real s1125 = (charge*s1124*s151)/2.;
				const Real s1126 = (charge*s1030*s152*s393)/4.;
				const Real s1127 = s1125 + s1126;
				const Real s1128 = -2*s10*s2*s3;
				const Real s1129 = -2*s0*s10*s8;
				const Real s1130 = 2*s26*s4;
				const Real s1131 = 2*s18*s4;
				const Real s1132 = 2*s23*s4;
				const Real s1133 = -2*s10*s11*s6;
				const Real s1134 = 2*s10*s3*s9;
				const Real s1135 = 2*s10*s2*s9;
				const Real s1136 = -4*s3*s4*s9;
				const Real s1137 = -2*s10*s53;
				const Real s1138 = 2*s4*s53;
				const Real s1139 = 2*s10*s12*s8;
				const Real s1140 = 2*s0*s10*s12;
				const Real s1141 = -4*s12*s4*s8;
				const Real s1142 = -2*s10*s13;
				const Real s1143 = 2*s13*s4;
				const Real s1144 = -2*s14*s26;
				const Real s1145 = -2*s14*s18;
				const Real s1146 = -2*s14*s23;
				const Real s1147 = 2*s14*s3*s9;
				const Real s1148 = -2*s14*s2*s9;
				const Real s1149 = 2*s12*s14*s8;
				const Real s1150 = -2*s0*s12*s14;
				const Real s1151 = 2*s10*s11*s16;
				const Real s1152 = -4*s11*s16*s4;
				const Real s1153 = 2*s10*s16*s6;
				const Real s1154 = 2*s11*s14*s16;
				const Real s1155 = -2*s14*s16*s6;
				const Real s1156 = -2*s10*s17;
				const Real s1157 = 2*s17*s4;
				const Real s1158 = s1128 + s1129 + s1130 + s1131 + s1132 + s1133 + s1134 + s1135 + s1136 + s1137 + s1138 + s1139 + s1140 + s1141 + s1142 + s1143 + s1144 + s1145 + s1146 + s1147 + s1148 + s1149 + s1150 + s1151 + s1152 + s1153 + s1154 + s1155 + s1156 + s1157 + s601 + s603 + s605;
				const Real s1159 = -2*s203*s207;
				const Real s1160 = 2*s194*s215;
				const Real s1161 = s1159 + s1160;
				const Real s1162 = -(s1161*s199*s210*s225);
				const Real s1163 = -(s1161*s188*s223*s225);
				const Real s1164 = 2*s194*s199*s218;
				const Real s1165 = -(s188*s203*s218);
				const Real s1166 = s1162 + s1163 + s1164 + s1165;
				const Real s1167 = s1166*s199;
				const Real s1168 = -(s1161*s199*s223*s225);
				const Real s1169 = -(s1161*s188*s215*s225);
				const Real s1170 = s1168 + s1169 + s972;
				const Real s1171 = s1170*s188;
				const Real s1172 = s1167 + s1171;
				const Real s1173 = (charge*s1172*s151)/2.;
				const Real s1174 = (charge*s1158*s152*s264)/4.;
				const Real s1175 = s1173 + s1174;
				const Real s1176 = s1166*s201;
				const Real s1177 = s1170*s191;
				const Real s1178 = s1176 + s1177;
				const Real s1179 = (charge*s1178*s151)/2.;
				const Real s1180 = (charge*s1158*s152*s273)/4.;
				const Real s1181 = s1179 + s1180;
				const Real s1182 = s1166*s203;
				const Real s1183 = s1170*s194;
				const Real s1184 = s1182 + s1183 + s260 + s261;
				const Real s1185 = (charge*s1184*s151)/2.;
				const Real s1186 = (charge*s1158*s152*s282)/4.;
				const Real s1187 = s1185 + s1186;
				const Real s1188 = s1166*s205;
				const Real s1189 = s1170*s197;
				const Real s1190 = s1188 + s1189;
				const Real s1191 = (charge*s1190*s151)/2.;
				const Real s1192 = (charge*s1158*s152*s291)/4.;
				const Real s1193 = s1191 + s1192;
				const Real s1194 = -(s1161*s201*s210*s225);
				const Real s1195 = -(s1161*s191*s223*s225);
				const Real s1196 = 2*s194*s201*s218;
				const Real s1197 = -(s191*s203*s218);
				const Real s1198 = s1194 + s1195 + s1196 + s1197;
				const Real s1199 = s1198*s201;
				const Real s1200 = -(s1161*s201*s223*s225);
				const Real s1201 = -(s1161*s191*s215*s225);
				const Real s1202 = s1100 + s1200 + s1201;
				const Real s1203 = s1202*s191;
				const Real s1204 = s1199 + s1203;
				const Real s1205 = (charge*s1204*s151)/2.;
				const Real s1206 = (charge*s1158*s152*s316)/4.;
				const Real s1207 = s1205 + s1206;
				const Real s1208 = s1198*s203;
				const Real s1209 = s1202*s194;
				const Real s1210 = s1208 + s1209 + s312 + s313;
				const Real s1211 = (charge*s1210*s151)/2.;
				const Real s1212 = (charge*s1158*s152*s325)/4.;
				const Real s1213 = s1211 + s1212;
				const Real s1214 = s1198*s205;
				const Real s1215 = s1202*s197;
				const Real s1216 = s1214 + s1215;
				const Real s1217 = (charge*s1216*s151)/2.;
				const Real s1218 = (charge*s1158*s152*s334)/4.;
				const Real s1219 = s1217 + s1218;
				const Real s1220 = -(s1161*s203*s210*s225);
				const Real s1221 = -(s1161*s194*s223*s225);
				const Real s1222 = s194*s203*s218;
				const Real s1223 = s1220 + s1221 + s1222 + s906;
				const Real s1224 = s1223*s203;
				const Real s1225 = -(s1161*s203*s223*s225);
				const Real s1226 = -(s1161*s194*s215*s225);
				const Real s1227 = -(s213*s218);
				const Real s1228 = s1225 + s1226 + s1227 + s912;
				const Real s1229 = s1228*s194;
				const Real s1230 = s1224 + s1229 + s355 + s356;
				const Real s1231 = (charge*s1230*s151)/2.;
				const Real s1232 = (charge*s1158*s152*s359)/4.;
				const Real s1233 = s1231 + s1232;
				const Real s1234 = s1223*s205;
				const Real s1235 = s1228*s197;
				const Real s1236 = s1234 + s1235;
				const Real s1237 = (charge*s1236*s151)/2.;
				const Real s1238 = (charge*s1158*s152*s368)/4.;
				const Real s1239 = s1237 + s1238;
				const Real s1240 = -(s1161*s205*s210*s225);
				const Real s1241 = -(s1161*s197*s223*s225);
				const Real s1242 = -(s197*s203*s218);
				const Real s1243 = 2*s194*s205*s218;
				const Real s1244 = s1240 + s1241 + s1242 + s1243;
				const Real s1245 = s1244*s205;
				const Real s1246 = -(s1161*s205*s223*s225);
				const Real s1247 = -(s1161*s197*s215*s225);
				const Real s1248 = -(s203*s205*s218);
				const Real s1249 = s1246 + s1247 + s1248;
				const Real s1250 = s1249*s197;
				const Real s1251 = s1245 + s1250;
				const Real s1252 = (charge*s1251*s151)/2.;
				const Real s1253 = (charge*s1158*s152*s393)/4.;
				const Real s1254 = s1252 + s1253;
				const Real s1255 = -2*s11*s2*s3;
				const Real s1256 = -2*s0*s11*s8;
				const Real s1257 = -2*s10*s11*s4;
				const Real s1258 = 2*s26*s6;
				const Real s1259 = 2*s18*s6;
				const Real s1260 = 2*s21*s6;
				const Real s1261 = 2*s11*s3*s9;
				const Real s1262 = 2*s11*s2*s9;
				const Real s1263 = -4*s3*s6*s9;
				const Real s1264 = -2*s11*s53;
				const Real s1265 = 2*s53*s6;
				const Real s1266 = 2*s11*s12*s8;
				const Real s1267 = 2*s0*s11*s12;
				const Real s1268 = -4*s12*s6*s8;
				const Real s1269 = -2*s11*s13;
				const Real s1270 = 2*s13*s6;
				const Real s1271 = 2*s10*s11*s14;
				const Real s1272 = 2*s11*s14*s4;
				const Real s1273 = -4*s10*s14*s6;
				const Real s1274 = -2*s11*s15;
				const Real s1275 = 2*s15*s6;
				const Real s1276 = -2*s16*s26;
				const Real s1277 = -2*s16*s18;
				const Real s1278 = -2*s16*s21;
				const Real s1279 = 2*s16*s3*s9;
				const Real s1280 = -2*s16*s2*s9;
				const Real s1281 = 2*s12*s16*s8;
				const Real s1282 = -2*s0*s12*s16;
				const Real s1283 = 2*s10*s14*s16;
				const Real s1284 = -2*s14*s16*s4;
				const Real s1285 = s1255 + s1256 + s1257 + s1258 + s1259 + s1260 + s1261 + s1262 + s1263 + s1264 + s1265 + s1266 + s1267 + s1268 + s1269 + s1270 + s1271 + s1272 + s1273 + s1274 + s1275 + s1276 + s1277 + s1278 + s1279 + s1280 + s1281 + s1282 + s1283 + s1284 + s749 + s751 + s753;
				const Real s1286 = -2*s205*s207;
				const Real s1287 = 2*s197*s215;
				const Real s1288 = s1286 + s1287;
				const Real s1289 = -(s1288*s199*s210*s225);
				const Real s1290 = -(s1288*s188*s223*s225);
				const Real s1291 = 2*s197*s199*s218;
				const Real s1292 = -(s188*s205*s218);
				const Real s1293 = s1289 + s1290 + s1291 + s1292;
				const Real s1294 = s1293*s199;
				const Real s1295 = -(s1288*s199*s223*s225);
				const Real s1296 = -(s1288*s188*s215*s225);
				const Real s1297 = s1295 + s1296 + s993;
				const Real s1298 = s1297*s188;
				const Real s1299 = s1294 + s1298;
				const Real s1300 = (charge*s1299*s151)/2.;
				const Real s1301 = (charge*s1285*s152*s264)/4.;
				const Real s1302 = s1300 + s1301;
				const Real s1303 = s1293*s201;
				const Real s1304 = s1297*s191;
				const Real s1305 = s1303 + s1304;
				const Real s1306 = (charge*s1305*s151)/2.;
				const Real s1307 = (charge*s1285*s152*s273)/4.;
				const Real s1308 = s1306 + s1307;
				const Real s1309 = s1293*s203;
				const Real s1310 = s1297*s194;
				const Real s1311 = s1309 + s1310;
				const Real s1312 = (charge*s1311*s151)/2.;
				const Real s1313 = (charge*s1285*s152*s282)/4.;
				const Real s1314 = s1312 + s1313;
				const Real s1315 = s1293*s205;
				const Real s1316 = s1297*s197;
				const Real s1317 = s1315 + s1316 + s260 + s261;
				const Real s1318 = (charge*s1317*s151)/2.;
				const Real s1319 = (charge*s1285*s152*s291)/4.;
				const Real s1320 = s1318 + s1319;
				const Real s1321 = -(s1288*s201*s210*s225);
				const Real s1322 = -(s1288*s191*s223*s225);
				const Real s1323 = 2*s197*s201*s218;
				const Real s1324 = -(s191*s205*s218);
				const Real s1325 = s1321 + s1322 + s1323 + s1324;
				const Real s1326 = s1325*s201;
				const Real s1327 = -(s1288*s201*s223*s225);
				const Real s1328 = -(s1288*s191*s215*s225);
				const Real s1329 = s1121 + s1327 + s1328;
				const Real s1330 = s1329*s191;
				const Real s1331 = s1326 + s1330;
				const Real s1332 = (charge*s1331*s151)/2.;
				const Real s1333 = (charge*s1285*s152*s316)/4.;
				const Real s1334 = s1332 + s1333;
				const Real s1335 = s1325*s203;
				const Real s1336 = s1329*s194;
				const Real s1337 = s1335 + s1336;
				const Real s1338 = (charge*s1337*s151)/2.;
				const Real s1339 = (charge*s1285*s152*s325)/4.;
				const Real s1340 = s1338 + s1339;
				const Real s1341 = s1325*s205;
				const Real s1342 = s1329*s197;
				const Real s1343 = s1341 + s1342 + s312 + s313;
				const Real s1344 = (charge*s1343*s151)/2.;
				const Real s1345 = (charge*s1285*s152*s334)/4.;
				const Real s1346 = s1344 + s1345;
				const Real s1347 = -(s1288*s203*s210*s225);
				const Real s1348 = -(s1288*s194*s223*s225);
				const Real s1349 = 2*s197*s203*s218;
				const Real s1350 = -(s194*s205*s218);
				const Real s1351 = s1347 + s1348 + s1349 + s1350;
				const Real s1352 = s1351*s203;
				const Real s1353 = -(s1288*s203*s223*s225);
				const Real s1354 = -(s1288*s194*s215*s225);
				const Real s1355 = s1248 + s1353 + s1354;
				const Real s1356 = s1355*s194;
				const Real s1357 = s1352 + s1356;
				const Real s1358 = (charge*s1357*s151)/2.;
				const Real s1359 = (charge*s1285*s152*s359)/4.;
				const Real s1360 = s1358 + s1359;
				const Real s1361 = s1351*s205;
				const Real s1362 = s1355*s197;
				const Real s1363 = s1361 + s1362 + s355 + s356;
				const Real s1364 = (charge*s1363*s151)/2.;
				const Real s1365 = (charge*s1285*s152*s368)/4.;
				const Real s1366 = s1364 + s1365;
				const Real s1367 = -(s1288*s205*s210*s225);
				const Real s1368 = -(s1288*s197*s223*s225);
				const Real s1369 = s197*s205*s218;
				const Real s1370 = s1367 + s1368 + s1369 + s906;
				const Real s1371 = s1370*s205;
				const Real s1372 = -(s1288*s205*s223*s225);
				const Real s1373 = -(s1288*s197*s215*s225);
				const Real s1374 = -(s214*s218);
				const Real s1375 = s1372 + s1373 + s1374 + s912;
				const Real s1376 = s1375*s197;
				const Real s1377 = s1371 + s1376 + s389 + s390;
				const Real s1378 = (charge*s1377*s151)/2.;
				const Real s1379 = (charge*s1285*s152*s393)/4.;
				const Real s1380 = s1378 + s1379;
				const Real s1381 = -2*s18*s2;
				const Real s1382 = -2*s2*s21;
				const Real s1383 = -2*s2*s23;
				const Real s1384 = 2*s0*s3*s8;
				const Real s1385 = 2*s0*s2*s8;
				const Real s1386 = -2*s1*s3;
				const Real s1387 = 2*s10*s3*s4;
				const Real s1388 = 2*s10*s2*s4;
				const Real s1389 = -2*s3*s5;
				const Real s1390 = 2*s11*s3*s6;
				const Real s1391 = 2*s11*s2*s6;
				const Real s1392 = -2*s3*s7;
				const Real s1393 = 2*s18*s9;
				const Real s1394 = 2*s21*s9;
				const Real s1395 = 2*s23*s9;
				const Real s1396 = -4*s0*s8*s9;
				const Real s1397 = 2*s1*s9;
				const Real s1398 = -4*s10*s4*s9;
				const Real s1399 = 2*s5*s9;
				const Real s1400 = -4*s11*s6*s9;
				const Real s1401 = 2*s7*s9;
				const Real s1402 = -2*s12*s3*s8;
				const Real s1403 = -2*s0*s12*s2;
				const Real s1404 = -2*s10*s14*s3;
				const Real s1405 = -2*s14*s2*s4;
				const Real s1406 = -2*s11*s16*s3;
				const Real s1407 = -2*s16*s2*s6;
				const Real s1408 = s1381 + s1382 + s1383 + s1384 + s1385 + s1386 + s1387 + s1388 + s1389 + s1390 + s1391 + s1392 + s1393 + s1394 + s1395 + s1396 + s1397 + s1398 + s1399 + s1400 + s1401 + s1402 + s1403 + s1404 + s1405 + s1406 + s1407 + s165 + s172 + s179 + s880 + s887 + s894;
				const Real s1409 = 2*s199*s210;
				const Real s1410 = -2*s188*s207;
				const Real s1411 = s1409 + s1410;
				const Real s1412 = -(s1411*s199*s210*s225);
				const Real s1413 = -(s1411*s188*s223*s225);
				const Real s1414 = -(s189*s218);
				const Real s1415 = s210*s218;
				const Real s1416 = s1412 + s1413 + s1414 + s1415;
				const Real s1417 = s1416*s199;
				const Real s1418 = -(s1411*s199*s223*s225);
				const Real s1419 = -(s1411*s188*s215*s225);
				const Real s1420 = s1418 + s1419 + s905 + s906;
				const Real s1421 = s1420*s188;
				const Real s1422 = s1417 + s1421 + s256 + s257;
				const Real s1423 = (charge*s1422*s151)/2.;
				const Real s1424 = (charge*s1408*s152*s264)/4.;
				const Real s1425 = s1423 + s1424;
				const Real s1426 = s1416*s201;
				const Real s1427 = s1420*s191;
				const Real s1428 = s1426 + s1427;
				const Real s1429 = (charge*s1428*s151)/2.;
				const Real s1430 = (charge*s1408*s152*s273)/4.;
				const Real s1431 = s1429 + s1430;
				const Real s1432 = s1416*s203;
				const Real s1433 = s1420*s194;
				const Real s1434 = s1432 + s1433;
				const Real s1435 = (charge*s1434*s151)/2.;
				const Real s1436 = (charge*s1408*s152*s282)/4.;
				const Real s1437 = s1435 + s1436;
				const Real s1438 = s1416*s205;
				const Real s1439 = s1420*s197;
				const Real s1440 = s1438 + s1439;
				const Real s1441 = (charge*s1440*s151)/2.;
				const Real s1442 = (charge*s1408*s152*s291)/4.;
				const Real s1443 = s1441 + s1442;
				const Real s1444 = -(s1411*s201*s210*s225);
				const Real s1445 = -(s1411*s191*s223*s225);
				const Real s1446 = -(s188*s191*s218);
				const Real s1447 = s1444 + s1445 + s1446;
				const Real s1448 = s1447*s201;
				const Real s1449 = -(s1411*s201*s223*s225);
				const Real s1450 = -(s1411*s191*s215*s225);
				const Real s1451 = s1036 + s1037 + s1449 + s1450;
				const Real s1452 = s1451*s191;
				const Real s1453 = s1448 + s1452;
				const Real s1454 = (charge*s1453*s151)/2.;
				const Real s1455 = (charge*s1408*s152*s316)/4.;
				const Real s1456 = s1454 + s1455;
				const Real s1457 = s1447*s203;
				const Real s1458 = s1451*s194;
				const Real s1459 = s1457 + s1458;
				const Real s1460 = (charge*s1459*s151)/2.;
				const Real s1461 = (charge*s1408*s152*s325)/4.;
				const Real s1462 = s1460 + s1461;
				const Real s1463 = s1447*s205;
				const Real s1464 = s1451*s197;
				const Real s1465 = s1463 + s1464;
				const Real s1466 = (charge*s1465*s151)/2.;
				const Real s1467 = (charge*s1408*s152*s334)/4.;
				const Real s1468 = s1466 + s1467;
				const Real s1469 = -(s1411*s203*s210*s225);
				const Real s1470 = -(s1411*s194*s223*s225);
				const Real s1471 = -(s188*s194*s218);
				const Real s1472 = s1469 + s1470 + s1471;
				const Real s1473 = s1472*s203;
				const Real s1474 = -(s1411*s203*s223*s225);
				const Real s1475 = -(s1411*s194*s215*s225);
				const Real s1476 = s1164 + s1165 + s1474 + s1475;
				const Real s1477 = s1476*s194;
				const Real s1478 = s1473 + s1477;
				const Real s1479 = (charge*s1478*s151)/2.;
				const Real s1480 = (charge*s1408*s152*s359)/4.;
				const Real s1481 = s1479 + s1480;
				const Real s1482 = s1472*s205;
				const Real s1483 = s1476*s197;
				const Real s1484 = s1482 + s1483;
				const Real s1485 = (charge*s1484*s151)/2.;
				const Real s1486 = (charge*s1408*s152*s368)/4.;
				const Real s1487 = s1485 + s1486;
				const Real s1488 = -(s1411*s205*s210*s225);
				const Real s1489 = -(s1411*s197*s223*s225);
				const Real s1490 = -(s188*s197*s218);
				const Real s1491 = s1488 + s1489 + s1490;
				const Real s1492 = s1491*s205;
				const Real s1493 = -(s1411*s205*s223*s225);
				const Real s1494 = -(s1411*s197*s215*s225);
				const Real s1495 = s1291 + s1292 + s1493 + s1494;
				const Real s1496 = s1495*s197;
				const Real s1497 = s1492 + s1496;
				const Real s1498 = (charge*s1497*s151)/2.;
				const Real s1499 = (charge*s1408*s152*s393)/4.;
				const Real s1500 = s1498 + s1499;
				const Real s1501 = 2*s2*s3*s8;
				const Real s1502 = -2*s19*s8;
				const Real s1503 = -2*s0*s26;
				const Real s1504 = -2*s0*s21;
				const Real s1505 = -2*s0*s23;
				const Real s1506 = 2*s0*s2*s3;
				const Real s1507 = 2*s10*s4*s8;
				const Real s1508 = 2*s0*s10*s4;
				const Real s1509 = -2*s5*s8;
				const Real s1510 = 2*s11*s6*s8;
				const Real s1511 = 2*s0*s11*s6;
				const Real s1512 = -2*s7*s8;
				const Real s1513 = -2*s3*s8*s9;
				const Real s1514 = -2*s0*s2*s9;
				const Real s1515 = 2*s12*s26;
				const Real s1516 = 2*s12*s21;
				const Real s1517 = 2*s12*s23;
				const Real s1518 = -4*s12*s2*s3;
				const Real s1519 = 2*s12*s19;
				const Real s1520 = -4*s10*s12*s4;
				const Real s1521 = 2*s12*s5;
				const Real s1522 = -4*s11*s12*s6;
				const Real s1523 = 2*s12*s7;
				const Real s1524 = -2*s10*s14*s8;
				const Real s1525 = -2*s0*s14*s4;
				const Real s1526 = -2*s11*s16*s8;
				const Real s1527 = -2*s0*s16*s6;
				const Real s1528 = s1007 + s1018 + s1025 + s1501 + s1502 + s1503 + s1504 + s1505 + s1506 + s1507 + s1508 + s1509 + s1510 + s1511 + s1512 + s1513 + s1514 + s1515 + s1516 + s1517 + s1518 + s1519 + s1520 + s1521 + s1522 + s1523 + s1524 + s1525 + s1526 + s1527 + s405 + s417 + s424;
				const Real s1529 = 2*s201*s210;
				const Real s1530 = -2*s191*s207;
				const Real s1531 = s1529 + s1530;
				const Real s1532 = -(s1531*s199*s210*s225);
				const Real s1533 = -(s1531*s188*s223*s225);
				const Real s1534 = s1446 + s1532 + s1533;
				const Real s1535 = s1534*s199;
				const Real s1536 = -(s1531*s199*s223*s225);
				const Real s1537 = -(s1531*s188*s215*s225);
				const Real s1538 = s1536 + s1537 + s939 + s940;
				const Real s1539 = s1538*s188;
				const Real s1540 = s1535 + s1539;
				const Real s1541 = (charge*s151*s1540)/2.;
				const Real s1542 = (charge*s152*s1528*s264)/4.;
				const Real s1543 = s1541 + s1542;
				const Real s1544 = s1534*s201;
				const Real s1545 = s1538*s191;
				const Real s1546 = s1544 + s1545 + s256 + s257;
				const Real s1547 = (charge*s151*s1546)/2.;
				const Real s1548 = (charge*s152*s1528*s273)/4.;
				const Real s1549 = s1547 + s1548;
				const Real s1550 = s1534*s203;
				const Real s1551 = s1538*s194;
				const Real s1552 = s1550 + s1551;
				const Real s1553 = (charge*s151*s1552)/2.;
				const Real s1554 = (charge*s152*s1528*s282)/4.;
				const Real s1555 = s1553 + s1554;
				const Real s1556 = s1534*s205;
				const Real s1557 = s1538*s197;
				const Real s1558 = s1556 + s1557;
				const Real s1559 = (charge*s151*s1558)/2.;
				const Real s1560 = (charge*s152*s1528*s291)/4.;
				const Real s1561 = s1559 + s1560;
				const Real s1562 = -(s1531*s201*s210*s225);
				const Real s1563 = -(s1531*s191*s223*s225);
				const Real s1564 = -(s192*s218);
				const Real s1565 = s1415 + s1562 + s1563 + s1564;
				const Real s1566 = s1565*s201;
				const Real s1567 = -(s1531*s201*s223*s225);
				const Real s1568 = -(s1531*s191*s215*s225);
				const Real s1569 = s1068 + s1567 + s1568 + s906;
				const Real s1570 = s1569*s191;
				const Real s1571 = s1566 + s1570 + s308 + s309;
				const Real s1572 = (charge*s151*s1571)/2.;
				const Real s1573 = (charge*s152*s1528*s316)/4.;
				const Real s1574 = s1572 + s1573;
				const Real s1575 = s1565*s203;
				const Real s1576 = s1569*s194;
				const Real s1577 = s1575 + s1576;
				const Real s1578 = (charge*s151*s1577)/2.;
				const Real s1579 = (charge*s152*s1528*s325)/4.;
				const Real s1580 = s1578 + s1579;
				const Real s1581 = s1565*s205;
				const Real s1582 = s1569*s197;
				const Real s1583 = s1581 + s1582;
				const Real s1584 = (charge*s151*s1583)/2.;
				const Real s1585 = (charge*s152*s1528*s334)/4.;
				const Real s1586 = s1584 + s1585;
				const Real s1587 = -(s1531*s203*s210*s225);
				const Real s1588 = -(s1531*s194*s223*s225);
				const Real s1589 = -(s191*s194*s218);
				const Real s1590 = s1587 + s1588 + s1589;
				const Real s1591 = s1590*s203;
				const Real s1592 = -(s1531*s203*s223*s225);
				const Real s1593 = -(s1531*s194*s215*s225);
				const Real s1594 = s1196 + s1197 + s1592 + s1593;
				const Real s1595 = s1594*s194;
				const Real s1596 = s1591 + s1595;
				const Real s1597 = (charge*s151*s1596)/2.;
				const Real s1598 = (charge*s152*s1528*s359)/4.;
				const Real s1599 = s1597 + s1598;
				const Real s1600 = s1590*s205;
				const Real s1601 = s1594*s197;
				const Real s1602 = s1600 + s1601;
				const Real s1603 = (charge*s151*s1602)/2.;
				const Real s1604 = (charge*s152*s1528*s368)/4.;
				const Real s1605 = s1603 + s1604;
				const Real s1606 = -(s1531*s205*s210*s225);
				const Real s1607 = -(s1531*s197*s223*s225);
				const Real s1608 = -(s191*s197*s218);
				const Real s1609 = s1606 + s1607 + s1608;
				const Real s1610 = s1609*s205;
				const Real s1611 = -(s1531*s205*s223*s225);
				const Real s1612 = -(s1531*s197*s215*s225);
				const Real s1613 = s1323 + s1324 + s1611 + s1612;
				const Real s1614 = s1613*s197;
				const Real s1615 = s1610 + s1614;
				const Real s1616 = (charge*s151*s1615)/2.;
				const Real s1617 = (charge*s152*s1528*s393)/4.;
				const Real s1618 = s1616 + s1617;
				const Real s1619 = 2*s10*s2*s3;
				const Real s1620 = -2*s10*s19;
				const Real s1621 = 2*s0*s10*s8;
				const Real s1622 = -2*s1*s10;
				const Real s1623 = -2*s26*s4;
				const Real s1624 = -2*s18*s4;
				const Real s1625 = -2*s23*s4;
				const Real s1626 = 2*s2*s3*s4;
				const Real s1627 = 2*s0*s4*s8;
				const Real s1628 = 2*s10*s11*s6;
				const Real s1629 = 2*s11*s4*s6;
				const Real s1630 = -2*s10*s7;
				const Real s1631 = -2*s10*s3*s9;
				const Real s1632 = -2*s2*s4*s9;
				const Real s1633 = -2*s10*s12*s8;
				const Real s1634 = -2*s0*s12*s4;
				const Real s1635 = 2*s14*s26;
				const Real s1636 = 2*s14*s18;
				const Real s1637 = 2*s14*s23;
				const Real s1638 = -4*s14*s2*s3;
				const Real s1639 = 2*s14*s19;
				const Real s1640 = -4*s0*s14*s8;
				const Real s1641 = 2*s1*s14;
				const Real s1642 = -4*s11*s14*s6;
				const Real s1643 = 2*s14*s7;
				const Real s1644 = -2*s10*s11*s16;
				const Real s1645 = -2*s16*s4*s6;
				const Real s1646 = s1135 + s1140 + s1153 + s1619 + s1620 + s1621 + s1622 + s1623 + s1624 + s1625 + s1626 + s1627 + s1628 + s1629 + s1630 + s1631 + s1632 + s1633 + s1634 + s1635 + s1636 + s1637 + s1638 + s1639 + s1640 + s1641 + s1642 + s1643 + s1644 + s1645 + s592 + s597 + s611;
				const Real s1647 = 2*s203*s210;
				const Real s1648 = -2*s194*s207;
				const Real s1649 = s1647 + s1648;
				const Real s1650 = -(s1649*s199*s210*s225);
				const Real s1651 = -(s1649*s188*s223*s225);
				const Real s1652 = s1471 + s1650 + s1651;
				const Real s1653 = s1652*s199;
				const Real s1654 = -(s1649*s199*s223*s225);
				const Real s1655 = -(s1649*s188*s215*s225);
				const Real s1656 = s1654 + s1655 + s966 + s967;
				const Real s1657 = s1656*s188;
				const Real s1658 = s1653 + s1657;
				const Real s1659 = (charge*s151*s1658)/2.;
				const Real s1660 = (charge*s152*s1646*s264)/4.;
				const Real s1661 = s1659 + s1660;
				const Real s1662 = s1652*s201;
				const Real s1663 = s1656*s191;
				const Real s1664 = s1662 + s1663;
				const Real s1665 = (charge*s151*s1664)/2.;
				const Real s1666 = (charge*s152*s1646*s273)/4.;
				const Real s1667 = s1665 + s1666;
				const Real s1668 = s1652*s203;
				const Real s1669 = s1656*s194;
				const Real s1670 = s1668 + s1669 + s256 + s257;
				const Real s1671 = (charge*s151*s1670)/2.;
				const Real s1672 = (charge*s152*s1646*s282)/4.;
				const Real s1673 = s1671 + s1672;
				const Real s1674 = s1652*s205;
				const Real s1675 = s1656*s197;
				const Real s1676 = s1674 + s1675;
				const Real s1677 = (charge*s151*s1676)/2.;
				const Real s1678 = (charge*s152*s1646*s291)/4.;
				const Real s1679 = s1677 + s1678;
				const Real s1680 = -(s1649*s201*s210*s225);
				const Real s1681 = -(s1649*s191*s223*s225);
				const Real s1682 = s1589 + s1680 + s1681;
				const Real s1683 = s1682*s201;
				const Real s1684 = -(s1649*s201*s223*s225);
				const Real s1685 = -(s1649*s191*s215*s225);
				const Real s1686 = s1094 + s1095 + s1684 + s1685;
				const Real s1687 = s1686*s191;
				const Real s1688 = s1683 + s1687;
				const Real s1689 = (charge*s151*s1688)/2.;
				const Real s1690 = (charge*s152*s1646*s316)/4.;
				const Real s1691 = s1689 + s1690;
				const Real s1692 = s1682*s203;
				const Real s1693 = s1686*s194;
				const Real s1694 = s1692 + s1693 + s308 + s309;
				const Real s1695 = (charge*s151*s1694)/2.;
				const Real s1696 = (charge*s152*s1646*s325)/4.;
				const Real s1697 = s1695 + s1696;
				const Real s1698 = s1682*s205;
				const Real s1699 = s1686*s197;
				const Real s1700 = s1698 + s1699;
				const Real s1701 = (charge*s151*s1700)/2.;
				const Real s1702 = (charge*s152*s1646*s334)/4.;
				const Real s1703 = s1701 + s1702;
				const Real s1704 = -(s1649*s203*s210*s225);
				const Real s1705 = -(s1649*s194*s223*s225);
				const Real s1706 = -(s195*s218);
				const Real s1707 = s1415 + s1704 + s1705 + s1706;
				const Real s1708 = s1707*s203;
				const Real s1709 = -(s1649*s203*s223*s225);
				const Real s1710 = -(s1649*s194*s215*s225);
				const Real s1711 = s1222 + s1709 + s1710 + s906;
				const Real s1712 = s1711*s194;
				const Real s1713 = s1708 + s1712 + s351 + s352;
				const Real s1714 = (charge*s151*s1713)/2.;
				const Real s1715 = (charge*s152*s1646*s359)/4.;
				const Real s1716 = s1714 + s1715;
				const Real s1717 = s1707*s205;
				const Real s1718 = s1711*s197;
				const Real s1719 = s1717 + s1718;
				const Real s1720 = (charge*s151*s1719)/2.;
				const Real s1721 = (charge*s152*s1646*s368)/4.;
				const Real s1722 = s1720 + s1721;
				const Real s1723 = -(s1649*s205*s210*s225);
				const Real s1724 = -(s1649*s197*s223*s225);
				const Real s1725 = -(s194*s197*s218);
				const Real s1726 = s1723 + s1724 + s1725;
				const Real s1727 = s1726*s205;
				const Real s1728 = -(s1649*s205*s223*s225);
				const Real s1729 = -(s1649*s197*s215*s225);
				const Real s1730 = s1349 + s1350 + s1728 + s1729;
				const Real s1731 = s1730*s197;
				const Real s1732 = s1727 + s1731;
				const Real s1733 = (charge*s151*s1732)/2.;
				const Real s1734 = (charge*s152*s1646*s393)/4.;
				const Real s1735 = s1733 + s1734;
				const Real s1736 = 2*s11*s2*s3;
				const Real s1737 = -2*s11*s19;
				const Real s1738 = 2*s0*s11*s8;
				const Real s1739 = -2*s1*s11;
				const Real s1740 = 2*s10*s11*s4;
				const Real s1741 = -2*s11*s5;
				const Real s1742 = -2*s26*s6;
				const Real s1743 = -2*s18*s6;
				const Real s1744 = -2*s21*s6;
				const Real s1745 = 2*s2*s3*s6;
				const Real s1746 = 2*s0*s6*s8;
				const Real s1747 = 2*s10*s4*s6;
				const Real s1748 = -2*s11*s3*s9;
				const Real s1749 = -2*s2*s6*s9;
				const Real s1750 = -2*s11*s12*s8;
				const Real s1751 = -2*s0*s12*s6;
				const Real s1752 = -2*s10*s11*s14;
				const Real s1753 = -2*s14*s4*s6;
				const Real s1754 = 2*s16*s26;
				const Real s1755 = 2*s16*s18;
				const Real s1756 = 2*s16*s21;
				const Real s1757 = -4*s16*s2*s3;
				const Real s1758 = 2*s16*s19;
				const Real s1759 = -4*s0*s16*s8;
				const Real s1760 = 2*s1*s16;
				const Real s1761 = -4*s10*s16*s4;
				const Real s1762 = 2*s16*s5;
				const Real s1763 = s1262 + s1267 + s1272 + s1736 + s1737 + s1738 + s1739 + s1740 + s1741 + s1742 + s1743 + s1744 + s1745 + s1746 + s1747 + s1748 + s1749 + s1750 + s1751 + s1752 + s1753 + s1754 + s1755 + s1756 + s1757 + s1758 + s1759 + s1760 + s1761 + s1762 + s735 + s740 + s745;
				const Real s1764 = 2*s205*s210;
				const Real s1765 = -2*s197*s207;
				const Real s1766 = s1764 + s1765;
				const Real s1767 = -(s1766*s199*s210*s225);
				const Real s1768 = -(s1766*s188*s223*s225);
				const Real s1769 = s1490 + s1767 + s1768;
				const Real s1770 = s1769*s199;
				const Real s1771 = -(s1766*s199*s223*s225);
				const Real s1772 = -(s1766*s188*s215*s225);
				const Real s1773 = s1771 + s1772 + s987 + s988;
				const Real s1774 = s1773*s188;
				const Real s1775 = s1770 + s1774;
				const Real s1776 = (charge*s151*s1775)/2.;
				const Real s1777 = (charge*s152*s1763*s264)/4.;
				const Real s1778 = s1776 + s1777;
				const Real s1779 = s1769*s201;
				const Real s1780 = s1773*s191;
				const Real s1781 = s1779 + s1780;
				const Real s1782 = (charge*s151*s1781)/2.;
				const Real s1783 = (charge*s152*s1763*s273)/4.;
				const Real s1784 = s1782 + s1783;
				const Real s1785 = s1769*s203;
				const Real s1786 = s1773*s194;
				const Real s1787 = s1785 + s1786;
				const Real s1788 = (charge*s151*s1787)/2.;
				const Real s1789 = (charge*s152*s1763*s282)/4.;
				const Real s1790 = s1788 + s1789;
				const Real s1791 = s1769*s205;
				const Real s1792 = s1773*s197;
				const Real s1793 = s1791 + s1792 + s256 + s257;
				const Real s1794 = (charge*s151*s1793)/2.;
				const Real s1795 = (charge*s152*s1763*s291)/4.;
				const Real s1796 = s1794 + s1795;
				const Real s1797 = -(s1766*s201*s210*s225);
				const Real s1798 = -(s1766*s191*s223*s225);
				const Real s1799 = s1608 + s1797 + s1798;
				const Real s1800 = s1799*s201;
				const Real s1801 = -(s1766*s201*s223*s225);
				const Real s1802 = -(s1766*s191*s215*s225);
				const Real s1803 = s1115 + s1116 + s1801 + s1802;
				const Real s1804 = s1803*s191;
				const Real s1805 = s1800 + s1804;
				const Real s1806 = (charge*s151*s1805)/2.;
				const Real s1807 = (charge*s152*s1763*s316)/4.;
				const Real s1808 = s1806 + s1807;
				const Real s1809 = s1799*s203;
				const Real s1810 = s1803*s194;
				const Real s1811 = s1809 + s1810;
				const Real s1812 = (charge*s151*s1811)/2.;
				const Real s1813 = (charge*s152*s1763*s325)/4.;
				const Real s1814 = s1812 + s1813;
				const Real s1815 = s1799*s205;
				const Real s1816 = s1803*s197;
				const Real s1817 = s1815 + s1816 + s308 + s309;
				const Real s1818 = (charge*s151*s1817)/2.;
				const Real s1819 = (charge*s152*s1763*s334)/4.;
				const Real s1820 = s1818 + s1819;
				const Real s1821 = -(s1766*s203*s210*s225);
				const Real s1822 = -(s1766*s194*s223*s225);
				const Real s1823 = s1725 + s1821 + s1822;
				const Real s1824 = s1823*s203;
				const Real s1825 = -(s1766*s203*s223*s225);
				const Real s1826 = -(s1766*s194*s215*s225);
				const Real s1827 = s1242 + s1243 + s1825 + s1826;
				const Real s1828 = s1827*s194;
				const Real s1829 = s1824 + s1828;
				const Real s1830 = (charge*s151*s1829)/2.;
				const Real s1831 = (charge*s152*s1763*s359)/4.;
				const Real s1832 = s1830 + s1831;
				const Real s1833 = s1823*s205;
				const Real s1834 = s1827*s197;
				const Real s1835 = s1833 + s1834 + s351 + s352;
				const Real s1836 = (charge*s151*s1835)/2.;
				const Real s1837 = (charge*s152*s1763*s368)/4.;
				const Real s1838 = s1836 + s1837;
				const Real s1839 = -(s1766*s205*s210*s225);
				const Real s1840 = -(s1766*s197*s223*s225);
				const Real s1841 = -(s198*s218);
				const Real s1842 = s1415 + s1839 + s1840 + s1841;
				const Real s1843 = s1842*s205;
				const Real s1844 = -(s1766*s205*s223*s225);
				const Real s1845 = -(s1766*s197*s215*s225);
				const Real s1846 = s1369 + s1844 + s1845 + s906;
				const Real s1847 = s1846*s197;
				const Real s1848 = s1843 + s1847 + s385 + s386;
				const Real s1849 = (charge*s151*s1848)/2.;
				const Real s1850 = (charge*s152*s1763*s393)/4.;
				const Real s1851 = s1849 + s1850;
				B[12*i+0] = charge*((charge*s152*s186*s396)/4. + (charge*s152*s186*s397)/4. + (charge*s10*s152*s186*s433)/4. + (charge*s11*s152*s186*s434)/4. + (charge*s152*s186*s2*s435)/4. + (charge*s0*s152*s186*s436)/4. + (charge*s152*s186*s4*s437)/4. + (charge*s12*s152*s186*s442)/4. + (charge*s14*s152*s186*s443)/4. + (charge*s152*s186*s444*s445)/12. + (charge*s152*s16*s186*s446)/4. + (charge*s152*s186*s447*s448)/12. + (charge*s152*s186*s449*s450)/12. + s440*((charge*s152*s186*s441)/12. + s451) + s432*((charge*s152*s186*s3)/4. + s453) + s266*s454 + s266*s465 + s275*s480 + s275*s481 + s284*s488 + s284*s489 + s293*s496 + s293*s497 + s318*s504 + s318*s506 + s327*s525 + s327*s526 + s336*s533 + s336*s534 + s361*s541 + s361*s542 + s370*s559 + s370*s560 + s395*s567 + s395*s568 + (charge*s152*s186*s438*s6)/4. + (charge*s152*s186*s452*s8)/4. + (charge*s152*s186*s439*s9)/4.);
				B[12*i+1] = charge*((charge*s152*s396*s431)/4. + (charge*s152*s397*s431)/4. + (charge*s152*s3*s431*s432)/4. + (charge*s10*s152*s431*s433)/4. + (charge*s11*s152*s431*s434)/4. + (charge*s152*s2*s431*s435)/4. + (charge*s0*s152*s431*s436)/4. + (charge*s152*s4*s431*s437)/4. + (charge*s152*s431*s440*s441)/12. + (charge*s12*s152*s431*s442)/4. + (charge*s14*s152*s431*s443)/4. + (charge*s152*s431*s444*s445)/12. + (charge*s152*s16*s431*s446)/4. + (charge*s152*s431*s447*s448)/12. + s449*((charge*s152*s431*s450)/12. + s451) + s454*s479 + s465*s479 + s480*s487 + s481*s487 + s488*s495 + s489*s495 + s496*s503 + s497*s503 + s504*s524 + s506*s524 + s525*s532 + s526*s532 + s533*s540 + s534*s540 + s541*s558 + s542*s558 + s559*s566 + s560*s566 + s567*s584 + s568*s584 + (charge*s152*s431*s438*s6)/4. + s452*(s453 + (charge*s152*s431*s8)/4.) + (charge*s152*s431*s439*s9)/4.);
				B[12*i+2] = charge*((charge*s152*s396*s618)/4. + (charge*s152*s397*s618)/4. + (charge*s152*s3*s432*s618)/4. + (charge*s11*s152*s434*s618)/4. + (charge*s152*s2*s435*s618)/4. + (charge*s0*s152*s436*s618)/4. + (charge*s152*s4*s437*s618)/4. + (charge*s152*s440*s441*s618)/12. + (charge*s12*s152*s442*s618)/4. + (charge*s14*s152*s443*s618)/4. + (charge*s152*s16*s446*s618)/4. + (charge*s152*s447*s448*s618)/12. + (charge*s152*s449*s450*s618)/12. + (charge*s152*s438*s6*s618)/4. + s433*(s453 + (charge*s10*s152*s618)/4.) + s444*(s451 + (charge*s152*s445*s618)/12.) + s454*s642 + s465*s642 + s480*s648 + s481*s648 + s488*s654 + s489*s654 + s496*s660 + s497*s660 + s504*s674 + s506*s674 + s525*s680 + s526*s680 + s533*s686 + s534*s686 + s541*s705 + s542*s705 + s559*s711 + s560*s711 + s567*s727 + s568*s727 + (charge*s152*s452*s618*s8)/4. + (charge*s152*s439*s618*s9)/4.);
				B[12*i+3] = charge*((charge*s152*s396*s761)/4. + (charge*s152*s397*s761)/4. + (charge*s152*s3*s432*s761)/4. + (charge*s10*s152*s433*s761)/4. + (charge*s152*s2*s435*s761)/4. + (charge*s0*s152*s436*s761)/4. + (charge*s152*s4*s437*s761)/4. + (charge*s152*s440*s441*s761)/12. + (charge*s12*s152*s442*s761)/4. + (charge*s14*s152*s443*s761)/4. + (charge*s152*s444*s445*s761)/12. + (charge*s152*s16*s446*s761)/4. + (charge*s152*s449*s450*s761)/12. + (charge*s152*s438*s6*s761)/4. + s434*(s453 + (charge*s11*s152*s761)/4.) + s447*(s451 + (charge*s152*s448*s761)/12.) + s454*s785 + s465*s785 + s480*s791 + s481*s791 + s488*s797 + s489*s797 + (charge*s152*s452*s761*s8)/4. + s496*s803 + s497*s803 + s504*s817 + s506*s817 + s525*s823 + s526*s823 + s533*s829 + s534*s829 + s541*s843 + s542*s843 + s559*s849 + s560*s849 + s567*s868 + s568*s868 + (charge*s152*s439*s761*s9)/4.);
				B[12*i+4] = charge*((charge*s152*s396*s899)/4. + (charge*s152*s397*s899)/4. + (charge*s152*s3*s432*s899)/4. + (charge*s10*s152*s433*s899)/4. + (charge*s11*s152*s434*s899)/4. + (charge*s0*s152*s436*s899)/4. + (charge*s152*s4*s437*s899)/4. + (charge*s12*s152*s442*s899)/4. + (charge*s14*s152*s443*s899)/4. + (charge*s152*s444*s445*s899)/12. + (charge*s152*s16*s446*s899)/4. + (charge*s152*s447*s448*s899)/12. + (charge*s152*s449*s450*s899)/12. + (charge*s152*s438*s6*s899)/4. + (charge*s152*s452*s8*s899)/4. + s435*(s453 + (charge*s152*s2*s899)/4.) + s440*(s451 + (charge*s152*s441*s899)/12.) + (charge*s152*s439*s899*s9)/4. + s454*s918 + s465*s918 + s480*s924 + s481*s924 + s488*s930 + s489*s930 + s496*s936 + s497*s936 + s504*s951 + s506*s951 + s525*s957 + s526*s957 + s533*s963 + s534*s963 + s541*s978 + s542*s978 + s559*s984 + s560*s984 + s567*s999 + s568*s999);
				B[12*i+5] = charge*((charge*s1030*s152*s396)/4. + (charge*s1030*s152*s397)/4. + (charge*s1030*s152*s3*s432)/4. + (charge*s10*s1030*s152*s433)/4. + (charge*s1030*s11*s152*s434)/4. + (charge*s1030*s152*s2*s435)/4. + (charge*s1030*s152*s4*s437)/4. + (charge*s1030*s152*s440*s441)/12. + (charge*s1030*s12*s152*s442)/4. + (charge*s1030*s14*s152*s443)/4. + (charge*s1030*s152*s444*s445)/12. + (charge*s1030*s152*s16*s446)/4. + (charge*s1030*s152*s447*s448)/12. + s449*((charge*s1030*s152*s450)/12. + s451) + s436*((charge*s0*s1030*s152)/4. + s453) + s1047*s454 + s1047*s465 + s1053*s480 + s1053*s481 + s1059*s488 + s1059*s489 + s1065*s496 + s1065*s497 + s1079*s504 + s1079*s506 + s1085*s525 + s1085*s526 + s1091*s533 + s1091*s534 + s1106*s541 + s1106*s542 + s1112*s559 + s1112*s560 + s1127*s567 + s1127*s568 + (charge*s1030*s152*s438*s6)/4. + (charge*s1030*s152*s452*s8)/4. + (charge*s1030*s152*s439*s9)/4.);
				B[12*i+6] = charge*((charge*s1158*s152*s396)/4. + (charge*s1158*s152*s397)/4. + (charge*s1158*s152*s3*s432)/4. + (charge*s10*s1158*s152*s433)/4. + (charge*s11*s1158*s152*s434)/4. + (charge*s1158*s152*s2*s435)/4. + (charge*s0*s1158*s152*s436)/4. + (charge*s1158*s152*s440*s441)/12. + (charge*s1158*s12*s152*s442)/4. + (charge*s1158*s14*s152*s443)/4. + (charge*s1158*s152*s16*s446)/4. + (charge*s1158*s152*s447*s448)/12. + (charge*s1158*s152*s449*s450)/12. + s444*((charge*s1158*s152*s445)/12. + s451) + s437*((charge*s1158*s152*s4)/4. + s453) + s1175*s454 + s1175*s465 + s1181*s480 + s1181*s481 + s1187*s488 + s1187*s489 + s1193*s496 + s1193*s497 + s1207*s504 + s1207*s506 + s1213*s525 + s1213*s526 + s1219*s533 + s1219*s534 + s1233*s541 + s1233*s542 + s1239*s559 + s1239*s560 + s1254*s567 + s1254*s568 + (charge*s1158*s152*s438*s6)/4. + (charge*s1158*s152*s452*s8)/4. + (charge*s1158*s152*s439*s9)/4.);
				B[12*i+7] = charge*((charge*s1285*s152*s396)/4. + (charge*s1285*s152*s397)/4. + (charge*s1285*s152*s3*s432)/4. + (charge*s10*s1285*s152*s433)/4. + (charge*s11*s1285*s152*s434)/4. + (charge*s1285*s152*s2*s435)/4. + (charge*s0*s1285*s152*s436)/4. + (charge*s1285*s152*s4*s437)/4. + (charge*s1285*s152*s440*s441)/12. + (charge*s12*s1285*s152*s442)/4. + (charge*s1285*s14*s152*s443)/4. + (charge*s1285*s152*s444*s445)/12. + (charge*s1285*s152*s16*s446)/4. + (charge*s1285*s152*s449*s450)/12. + s447*((charge*s1285*s152*s448)/12. + s451) + s1302*s454 + s1302*s465 + s1308*s480 + s1308*s481 + s1314*s488 + s1314*s489 + s1320*s496 + s1320*s497 + s1334*s504 + s1334*s506 + s1340*s525 + s1340*s526 + s1346*s533 + s1346*s534 + s1360*s541 + s1360*s542 + s1366*s559 + s1366*s560 + s1380*s567 + s1380*s568 + s438*(s453 + (charge*s1285*s152*s6)/4.) + (charge*s1285*s152*s452*s8)/4. + (charge*s1285*s152*s439*s9)/4.);
				B[12*i+8] = charge*((charge*s1408*s152*s396)/4. + (charge*s1408*s152*s397)/4. + (charge*s1408*s152*s3*s432)/4. + (charge*s10*s1408*s152*s433)/4. + (charge*s11*s1408*s152*s434)/4. + (charge*s1408*s152*s2*s435)/4. + (charge*s0*s1408*s152*s436)/4. + (charge*s1408*s152*s4*s437)/4. + (charge*s12*s1408*s152*s442)/4. + (charge*s14*s1408*s152*s443)/4. + (charge*s1408*s152*s444*s445)/12. + (charge*s1408*s152*s16*s446)/4. + (charge*s1408*s152*s447*s448)/12. + (charge*s1408*s152*s449*s450)/12. + s440*((charge*s1408*s152*s441)/12. + s451) + s1425*s454 + s1425*s465 + s1431*s480 + s1431*s481 + s1437*s488 + s1437*s489 + s1443*s496 + s1443*s497 + s1456*s504 + s1456*s506 + s1462*s525 + s1462*s526 + s1468*s533 + s1468*s534 + s1481*s541 + s1481*s542 + s1487*s559 + s1487*s560 + s1500*s567 + s1500*s568 + (charge*s1408*s152*s438*s6)/4. + (charge*s1408*s152*s452*s8)/4. + s439*(s453 + (charge*s1408*s152*s9)/4.));
				B[12*i+9] = charge*((charge*s152*s1528*s396)/4. + (charge*s152*s1528*s397)/4. + (charge*s152*s1528*s3*s432)/4. + (charge*s10*s152*s1528*s433)/4. + (charge*s11*s152*s1528*s434)/4. + (charge*s152*s1528*s2*s435)/4. + (charge*s0*s152*s1528*s436)/4. + (charge*s152*s1528*s4*s437)/4. + (charge*s152*s1528*s440*s441)/12. + (charge*s14*s152*s1528*s443)/4. + (charge*s152*s1528*s444*s445)/12. + (charge*s152*s1528*s16*s446)/4. + (charge*s152*s1528*s447*s448)/12. + s449*((charge*s152*s1528*s450)/12. + s451) + s442*((charge*s12*s152*s1528)/4. + s453) + s1543*s454 + s1543*s465 + s1549*s480 + s1549*s481 + s1555*s488 + s1555*s489 + s1561*s496 + s1561*s497 + s1574*s504 + s1574*s506 + s1580*s525 + s1580*s526 + s1586*s533 + s1586*s534 + s1599*s541 + s1599*s542 + s1605*s559 + s1605*s560 + s1618*s567 + s1618*s568 + (charge*s152*s1528*s438*s6)/4. + (charge*s152*s1528*s452*s8)/4. + (charge*s152*s1528*s439*s9)/4.);
				B[12*i+10] = charge*((charge*s152*s1646*s396)/4. + (charge*s152*s1646*s397)/4. + (charge*s152*s1646*s3*s432)/4. + (charge*s10*s152*s1646*s433)/4. + (charge*s11*s152*s1646*s434)/4. + (charge*s152*s1646*s2*s435)/4. + (charge*s0*s152*s1646*s436)/4. + (charge*s152*s1646*s4*s437)/4. + (charge*s152*s1646*s440*s441)/12. + (charge*s12*s152*s1646*s442)/4. + (charge*s152*s16*s1646*s446)/4. + (charge*s152*s1646*s447*s448)/12. + (charge*s152*s1646*s449*s450)/12. + s444*((charge*s152*s1646*s445)/12. + s451) + s443*((charge*s14*s152*s1646)/4. + s453) + s1661*s454 + s1661*s465 + s1667*s480 + s1667*s481 + s1673*s488 + s1673*s489 + s1679*s496 + s1679*s497 + s1691*s504 + s1691*s506 + s1697*s525 + s1697*s526 + s1703*s533 + s1703*s534 + s1716*s541 + s1716*s542 + s1722*s559 + s1722*s560 + s1735*s567 + s1735*s568 + (charge*s152*s1646*s438*s6)/4. + (charge*s152*s1646*s452*s8)/4. + (charge*s152*s1646*s439*s9)/4.);
				B[12*i+11] = charge*((charge*s152*s1763*s396)/4. + (charge*s152*s1763*s397)/4. + (charge*s152*s1763*s3*s432)/4. + (charge*s10*s152*s1763*s433)/4. + (charge*s11*s152*s1763*s434)/4. + (charge*s152*s1763*s2*s435)/4. + (charge*s0*s152*s1763*s436)/4. + (charge*s152*s1763*s4*s437)/4. + (charge*s152*s1763*s440*s441)/12. + (charge*s12*s152*s1763*s442)/4. + (charge*s14*s152*s1763*s443)/4. + (charge*s152*s1763*s444*s445)/12. + (charge*s152*s1763*s449*s450)/12. + s447*((charge*s152*s1763*s448)/12. + s451) + s446*((charge*s152*s16*s1763)/4. + s453) + s1778*s454 + s1778*s465 + s1784*s480 + s1784*s481 + s1790*s488 + s1790*s489 + s1796*s496 + s1796*s497 + s1808*s504 + s1808*s506 + s1814*s525 + s1814*s526 + s1820*s533 + s1820*s534 + s1832*s541 + s1832*s542 + s1838*s559 + s1838*s560 + s1851*s567 + s1851*s568 + (charge*s152*s1763*s438*s6)/4. + (charge*s152*s1763*s452*s8)/4. + (charge*s152*s1763*s439*s9)/4.);
			},
			simplices.Dim(0),
			ThreadCount()
		);
    }

	}; // SimplicialMeshDetails<2,4,Real,Int>

} // namespace Repulsor