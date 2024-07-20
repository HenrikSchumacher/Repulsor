(* ::Package:: *)

ClearAll[ComputeNearFarData];


ComputeNearFarData[n_,m_]:=Module[{name},
name="ComputeNearFarData";
StringJoin[
"
	void ",name,"(
		const Tensor2<Real,Int> & V_coords,
		const Tensor2<Int ,Int> & simplices,
			  Tensor2<Real,Int> & P_near,
			  Tensor2<Real,Int> & P_far
	) const
    {
        ptic(ClassName()+\"::",name,"\");
        eprint(ClassName()+\"::",name," not implemented. Doing nothing.\");
        ptoc(ClassName()+\"::",name,"\");
    }
"]
];


ComputeNearFarData[0,m_]:=Module[{name,DomDim,AmbDim,HullSize,neardim,fardim,numdigits,tostr,n,s},
n=0;
s=ToString;
name="ComputeNearFarData";
DomDim=s[n];
AmbDim=s[m];
HullSize=s[m(n+1)];
neardim=1+(n+1)m+m (m+1)/2;
fardim=1+m+m (m+1)/2;
numdigits=Max[Length[IntegerDigits[m(n+1),10]],Length[IntegerDigits[Max[neardim,fardim],10]]];
tostr[i_]:=ConstantArray[" ",numdigits-Length[IntegerDigits[i,10]]]<>s[i];

StringJoin[
"
	void ",name,"( 
		const Tensor2<Real,Int> & V_coords,
		const Tensor2<Int ,Int> & simplices,
			  Tensor2<Real,Int> & P_near,
			  Tensor2<Real,Int> & P_far
	) const
    {
        ptic(ClassName()+\"::",name,"\");
        
        //Int size       = "<>s[n+1]<>";
        //Int amb_dim    = "<>AmbDim<>";
        //Int dom_dim    = "<>DomDim<>";

        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
		ParallelDo(
			[&]( const Int thread )
			{
				cptr<Real> V_coords__      = V_coords.data();	
				cptr<Int>  simplices__     = simplices.data();
	
				Real hull    [",s[n+1],"][",AmbDim,"];
	
				Int simplex  [",s[n+1],"];
				
				const Int i_begin = job_ptr[thread];
				const Int i_end   = job_ptr[thread+1];
	
	            for( Int i = i_begin; i < i_end; ++i )
	            {
					mptr<Real> near = P_near.data(i);                    
					mptr<Real> far  = P_far.data(i);   
            
",
Table[line[5,"simplex[",s[j],"] = simplices__[",s[n+1],"*i +",s[j],"]"],{j,0,n+1-1}],
"\n",
Table[line[5,"near[",s[1+m j+k],"] = hull[",s[j],"][",s[k],"] = V_coords__[",AmbDim,"*simplex[",s[j],"]+",s[k],"]"],{j,0,n+1-1},{k,0,m-1}],
"\n",
Table[line[5,"far["<>s[k+1]<>"] = ",floatcast@s[1./(n+1),CForm]," * ( ",Riffle[Table[{"hull[",s[j],"][",s[k],"]"},{j,0,n+1-1}]," + "]," )"],{k,0,m-1}],
"
					near[",tostr[0],"] = far[",tostr[0],"] = static_cast<Real>(1);
",
Module[{c=m},
Table[line[5,"near[",tostr[++c+n m],"] = far[",tostr[c],"] = ",floatcast@s[Boole[i==j]]],{i,0,m-1},{j,i,m-1}]],
"
	            } // for( Int i = i_begin; i < i_end; ++i )

			},
			ThreadCount()
		);

        ptoc(ClassName()+\"::",name,"\");
    }
"]
];


ComputeNearFarData[1,m_]:=Module[{name,DomDim,AmbDim,HullSize,neardim,fardim,numdigits,tostr,n,s},
n=1;
s=ToString;
name="ComputeNearFarData";
DomDim=s[n];
AmbDim=s[m];
HullSize=s[m(n+1)];
neardim=1+(n+1)m+m (m+1)/2;
fardim=1+m+m (m+1)/2;
numdigits=Max[Length[IntegerDigits[m(n+1),10]],Length[IntegerDigits[Max[neardim,fardim],10]]];
tostr[i_]:=ConstantArray[" ",numdigits-Length[IntegerDigits[i,10]]]<>s[i];

StringJoin[
"
	void ",name,"( 
		const Tensor2<Real,Int> & V_coords,
		const Tensor2<Int ,Int> & simplices,
			  Tensor2<Real,Int> & P_near,
			  Tensor2<Real,Int> & P_far
	) const
    {
        ptic(ClassName()+\"::",name,"\");
        
        //Int size       = "<>s[n+1]<>";
        //Int amb_dim    = "<>AmbDim<>";
        //Int dom_dim    = "<>DomDim<>";

        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
		ParallelDo(
			[&]( const Int thread )
			{
				cptr<Real> V_coords__      = V_coords.data();	
				cptr<Int > simplices__     = simplices.data();
	
				Real hull    [",s[n+1],"][",AmbDim,"];
				Real df      [",AmbDim,"][",DomDim,"];
				Real dfdagger[",DomDim,"][",AmbDim,"];
				Real g       [",DomDim,"][",DomDim,"];
				Real ginv    [",DomDim,"][",DomDim,"];
	
				Int simplex  [",s[n+1],"];
				
				const Int i_begin = job_ptr[thread];
				const Int i_end   = job_ptr[thread+1];
	
	            for( Int i = i_begin; i < i_end; ++i )
	            {
					mptr<Real> near = P_near.data(i);                    
					mptr<Real> far  = P_far.data(i);   
            
",
Table[line[5,"simplex[",s[j],"] = simplices__[",s[n+1],"*i +",s[j],"]"],{j,0,n+1-1}],
"\n",
Table[line[5,"near[",s[1+m j+k],"] = hull[",s[j],"][",s[k],"] = V_coords__[",AmbDim,"*simplex[",s[j],"]+",s[k],"]"],{j,0,n+1-1},{k,0,m-1}],
"\n",
Table[line[5,"far["<>s[k+1]<>"] = ",floatcast@s[1./(n+1),CForm]," * ( ",Riffle[Table[{"hull[",s[j],"][",s[k],"]"},{j,0,n+1-1}]," + "]," )"],{k,0,m-1}],
"\n",
Table[line[5,"df["<>s[k],"][",s[j]<>"] = hull[",s[j+1],"][",s[k],"] - hull[0][",s[k],"]"],{k,0,m-1},{j,0,n-1}],
"\n",
Table[line[5,"g[",s[i],"][",s[j],"] = ",Riffle[Table[{"df[",s[k],"][",s[i]<>"] * df[",s[k],"][",s[j],"]"},{k,0,m-1}]," + "]],{i,0,n-1},{j,0,n-1}],
"

					near[0] = far[0] = std::sqrt( std::fabs(g[0][0]) );
	
	                ginv[0][0] = static_cast<Real>(1)/g[0][0];
	                
	                //  dfdagger = g^{-1} * df^T (",DomDim," x ",AmbDim," matrix)
",
  Table[line[5,"dfdagger[",s[j],"][",s[k],"] = ",Riffle[Table[{"ginv[",s[Min[i,j]],"][",s[Max[i,j]],"] * df[",s[k],"][",s[i],"]"},{i,0,n-1}]," + "]],{j,0,n-1},{k,0,m-1}],"            
",
Module[{c=m},
Table[line[5,"near[",tostr[++c+n m],"] = far[",tostr[c],"]  = ",If[i==j,floatcast["1"],"  "],Table[{" - df[",s[i],"][",s[k],"] * dfdagger[",s[k],"][",s[j],"]"},{k,0,n-1}]],{i,0,m-1},{j,i,m-1}]],
"
	            } // for( Int i = i_begin; i < i_end; ++i )
			},
			ThreadCount()
		);

        ptoc(ClassName()+\"::",name,"\");
    }
"]
];


ComputeNearFarData[2,m_]:=Module[{name,DomDim,AmbDim,HullSize,neardim,fardim,numdigits,tostr,n,s},
n=2;
s=ToString;
name="ComputeNearFarData";
DomDim=s[n];
AmbDim=s[m];
HullSize=s[m(n+1)];
neardim=1+(n+1)m+m (m+1)/2;
fardim=1+m+m (m+1)/2;
numdigits=Max[Length[IntegerDigits[m(n+1),10]],Length[IntegerDigits[Max[neardim,fardim],10]]];
tostr[i_]:=ConstantArray[" ",numdigits-Length[IntegerDigits[i,10]]]<>s[i];

StringJoin[
"
	void ",name,"( 
		const Tensor2<Real,Int> & V_coords,
		const Tensor2<Int ,Int> & simplices,
			  Tensor2<Real,Int> & P_near,
			  Tensor2<Real,Int> & P_far
	) const
    {
        ptic(ClassName()+\"::",name,"\");
        
        //Int size       = "<>s[n+1]<>";
        //Int amb_dim    = "<>AmbDim<>";
        //Int dom_dim    = "<>DomDim<>";

        const JobPointers<Int> job_ptr ( simplices.Dimension(0), ThreadCount() );
        
		ParallelDo(
			[&]( const Int thread )
			{
				cptr<Real> V_coords__      = V_coords.data();	
				cptr<Int > simplices__     = simplices.data();
	
				Real hull    [",s[n+1],"][",AmbDim,"];
				Real df      [",AmbDim,"][",DomDim,"];
				Real dfdagger[",DomDim,"][",AmbDim,"];
				Real g       [",DomDim,"][",DomDim,"];
				Real ginv    [",DomDim,"][",DomDim,"];
	
				Int simplex  [",s[n+1],"];
				
				const Int i_begin = job_ptr[thread];
				const Int i_end   = job_ptr[thread+1];
	
	            for( Int i = i_begin; i < i_end; ++i )
	            {
					mptr<Real> near = P_near.data(i);                    
					mptr<Real> far  = P_far.data(i);   
            
",
Table[line[5,"simplex[",s[j],"] = simplices__[",s[n+1],"*i +",s[j],"]"],{j,0,n+1-1}],
"\n",
Table[line[5,"near[",s[1+m j+k],"] = hull[",s[j],"][",s[k],"] = V_coords__[",AmbDim,"*simplex[",s[j],"]+",s[k],"]"],{j,0,n+1-1},{k,0,m-1}],
"\n",
Table[line[5,"far["<>s[k+1]<>"] = ",floatcast@s[1./(n+1),CForm]," * ( ",Riffle[Table[{"hull[",s[j],"][",s[k],"]"},{j,0,n+1-1}]," + "]," )"],{k,0,m-1}],
"\n",
Table[line[5,"df["<>s[k],"][",s[j]<>"] = hull[",s[j+1],"][",s[k],"] - hull[0][",s[k],"]"],{k,0,m-1},{j,0,n-1}],
"\n",
Table[line[5,"g[",s[i],"][",s[j],"] = ",Riffle[Table[{"df[",s[k],"][",s[i]<>"] * df[",s[k],"][",s[j],"]"},{k,0,m-1}]," + "]],{i,0,n-1},{j,0,n-1}],
"
	                Real det = g[0][0] * g[1][1] - g[0][1] * g[0][1];
	
					near[0] = far[0] = std::sqrt( std::fabs(det) ) * ",floatcast["0.5"],";
	
	                Real invdet = ",floatcast["1"],"/det;
	                ginv[0][0] =  g[1][1] * invdet;
	                ginv[0][1] = -g[0][1] * invdet;
	                ginv[1][1] =  g[0][0] * invdet;
	                
	                //  dfdagger = g^{-1} * df^T (",DomDim," x ",AmbDim," matrix)
",
 Table[line[5,"dfdagger[",s[j],"][",s[k],"] = ",Riffle[Table[{"ginv[",s[Min[i,j]],"][",s[Max[i,j]],"] * df[",s[k],"][",s[i],"]"},{i,0,n-1}]," + "]],{j,0,n-1},{k,0,m-1}],"            
",
Module[{c=m},
	Table[line[5,"near[",tostr[++c+n m],"] = far[",tostr[c],"]  = ",If[i==j,floatcast["1"],"  "],Table[{" - df[",s[i],"][",s[k],"] * dfdagger[",s[k],"][",s[j],"]"},{k,0,n-1}]],{i,0,m-1},{j,i,m-1}]
],
"
	            } // for( Int i = i_begin; i < i_end; ++i )
			},
			ThreadCount()
		);

        ptoc(ClassName()+\"::",name,"\");
    }
"]
];
