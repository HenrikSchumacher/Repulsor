(* ::Package:: *)

ClearAll[ComputeNearFarDataOps];


ComputeNearFarDataOps[n_,m_]:=Module[{name},
name="ComputeNearFarDataOps";
StringJoin[
"
	void ",name,"( 
		const Tensor2<Real,Int> & V_coords,
        const Tensor2<Int ,Int> & simplices,
		      Tensor2<Real,Int> & P_coords,
		      Tensor3<Real,Int> & P_hull_coords,
		      Tensor2<Real,Int> & P_near,
		      Tensor2<Real,Int> & P_far,
		SparseMatrix_T & DiffOp,
		SparseMatrix_T & AvOp 
	) const
    {
        ptic(ClassName()+\"::",name,"\");
        eprint(ClassName()+\"::",name," not implemented. Doing nothing.\");
        ptoc(ClassName()+\"::",name,"\");
    }
"]
];


ComputeNearFarDataOps[0,m_]:=Module[{name,DomDim,AmbDim,HullSize,neardim,fardim,numdigits,tostr,n,s},
n=0;
s=ToString;
name="ComputeNearFarDataOps";
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
		      Tensor2<Real,Int> & P_coords,
		      Tensor3<Real,Int> & P_hull_coords,
		      Tensor2<Real,Int> & P_near,
		      Tensor2<Real,Int> & P_far,
		SparseMatrix_T & DiffOp,
		SparseMatrix_T & AvOp 
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
				mptr<LInt> AvOp_outer = AvOp.Outer().data();
				mptr< Int> AvOp_inner = AvOp.Inner().data();
				AvOp.Value().Fill(static_cast<Real>(1));
	
				mptr<LInt> DiffOp_outer = DiffOp.Outer().data();
				mptr< Int> DiffOp_inner = DiffOp.Inner().data();
				DiffOp.Value().SetZero();
	
				cptr<Real> V_coords__      = V_coords.data();
				
				cptr<Int>  simplices__     = simplices.data();
				mptr<Real> P_hull_coords__ = P_hull_coords.data();
			    mptr<Real> P_coords__      = P_coords.data();
	
				Int simplex        [",s[n+1],"];
				Int sorted_simplex [",s[n+1],"];
	
				const Int i_begin = job_ptr[thread];
				const Int i_end   = job_ptr[thread+1];
	
	            for( Int i = i_begin; i < i_end; ++i )
	            {
	
					mptr<Real> near = P_near.data(i);                    
					mptr<Real> far  = P_far.data(i);

",
Table[line[5,"simplex[",s[j],"] = sorted_simplex[",s[j],"] = simplices__[",s[n+1],"*i +",s[j],"]"],{j,0,n+1-1}],
"                  
	                // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
	                Sort( sorted_simplex, sorted_simplex + ",s[n+1]," );
	
					AvOp_outer[i+1] = (i+1) * ",s[n+1],";  
                      
",
Table[line[5,"AvOp_inner[",s[n+1],"*i+",s[j],"] = sorted_simplex[",s[j],"]"],{j,0,n+1-1}],
"\n",
Table[line[5,"DiffOp_outer[",AmbDim,"*i+",s[k],"] = (",AmbDim," * i + ",s[k],") * ",s[n+1]],{k,0,m-1}],
"\n",
Table[line[5,"DiffOp_inner[(i*",AmbDim,"+",s[k],")*",s[n+1],"+",s[j],"] = sorted_simplex[",s[j],"]"],{j,0,n+1-1},{k,0,m-1}],
"\n",
Table[line[5,"near[",s[1+m j+k],"] = P_hull_coords__[",HullSize,"*i+",s[m j+k],"] = V_coords__[",AmbDim,"*simplex[",s[j],"]+",s[k],"]"],{j,0,n+1-1},{k,0,m-1}],
"\n",
Table[line[5,"far["<>s[k+1]<>"] = P_coords__[",AmbDim,"*i+",s[k],"] = ",s[1./(n+1),CForm]," * ( ",Riffle[Table[{"P_hull_coords__[",HullSize,"*i+",s[m j+k],"]"},{j,0,n+1-1}]," + "]," )"],{k,0,m-1}],
"
	                near[",tostr[0],"] = far[",tostr[0],"] = ",floatcast[s[1]],";
",
Module[{c=m},
	Table[line[5,"near[",tostr[++c+n m],"] = far[",tostr[c],"] = ",floatcast[s[Boole[i==j]]]],{i,0,m-1},{j,i,m-1}]
],
"
	            }
			},
			ThreadCount()
		);

        ptoc(ClassName()+\"::",name,"\");
    }
"
]
];


ComputeNearFarDataOps[1,m_]:=Module[{name,DomDim,AmbDim,HullSize,neardim,fardim,numdigits,tostr,n,s},
n=1;
s=ToString;
name="ComputeNearFarDataOps";
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
		      Tensor2<Real,Int> & P_coords,
		      Tensor3<Real,Int> & P_hull_coords,
		      Tensor2<Real,Int> & P_near,
		      Tensor2<Real,Int> & P_far,
		SparseMatrix_T & DiffOp,
		SparseMatrix_T & AvOp 
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
				mptr<LInt> AvOp_outer = AvOp.Outer().data();
				mptr< Int> AvOp_inner = AvOp.Inner().data();
				mptr<Real> AvOp_value = AvOp.Values().data();
	
				mptr<LInt> DiffOp_outer = DiffOp.Outer().data();
				mptr< Int> DiffOp_inner = DiffOp.Inner().data();
				mptr<Real> DiffOp_value = DiffOp.Value().data();
	
				cptr<Real> V_coords__      = V_coords.data();
				
				cptr<Int>  simplices__     = simplices.data();
			    mptr<Real> P_hull_coords__ = P_hull_coords.data();
				mptr<Real> P_coords__      = P_coords.data();
	
				Real df       [",AmbDim,"][",DomDim,"];
				Real dfdagger [",DomDim,"][",AmbDim,"];
				Real g        [",DomDim,"][",DomDim,"];
				Real ginv     [",DomDim,"][",DomDim,"];
	
				Int simplex        [",s[n+1],"];
				Int sorted_simplex [",s[n+1],"];
	
				const Int i_begin = job_ptr[thread];
				const Int i_end   = job_ptr[thread+1];
	
	            for( Int i = i_begin; i < i_end; ++i )
	            {
	
					mptr<Real> near = P_near.data(i);                    
					mptr<Real> far  = P_far.data(i);

",
Table[line[5,"simplex[",s[j],"] = sorted_simplex[",s[j],"] = simplices__[",s[n+1],"*i +",s[j],"]"],{j,0,n+1-1}],
"                  
	                // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
	                Sort( sorted_simplex, sorted_simplex + ",s[n+1]," );
	
					AvOp_outer[i+1] = (i+1) * ",s[n+1],";  
	                      
",
Table[line[5,"AvOp_inner[",s[n+1],"*i+",s[j],"] = sorted_simplex[",s[j],"]"],{j,0,n+1-1}],
"\n",
Table[line[5,"AvOp_value[",s[n+1],"*i+",s[j],"] = ",s[1./(n+1),CForm]],{j,0,n+1-1}],
"\n",
Table[line[5,"DiffOp_outer[",AmbDim,"*i+",s[k],"] = (",AmbDim," * i + ",s[k],") * ",s[n+1]],{k,0,m-1}],
"\n",
Table[line[5,"DiffOp_inner[(i*",AmbDim,"+",s[k],")*",s[n+1],"+",s[j],"] = sorted_simplex[",s[j],"]"],{j,0,n+1-1},{k,0,m-1}],
"\n",
Table[line[5,"near[",s[1+m j+k],"] = P_hull_coords__[",HullSize,"*i+",s[m j+k],"] = V_coords__[",AmbDim,"*simplex[",s[j],"]+",s[k],"]"],{j,0,n+1-1},{k,0,m-1}],
"\n",
Table[line[5,"far["<>s[k+1]<>"] = P_coords__[",AmbDim,"*i+",s[k],"] = ",s[1./(n+1),CForm]," * ( ",Riffle[Table[{"P_hull_coords__[",HullSize,"*i+",s[m j+k],"]"},{j,0,n+1-1}]," + "]," )"],{k,0,m-1}],
"\n",
Table[line[5,"df["<>s[k],"][",s[j]<>"] = V_coords__[",AmbDim,"*sorted_simplex[",s[j+1],"]+",s[k],"] - V_coords__[",AmbDim,"*sorted_simplex[0]+",s[k],"]"],{k,0,m-1},{j,0,n-1}],
"\n",
Table[line[5,"g[",s[i],"][",s[j],"] = ",Riffle[Table[{"df[",s[k],"][",s[i]<>"] * df[",s[k],"][",s[j],"]"},{k,0,m-1}]," + "]],{i,0,n-1},{j,0,n-1}],
"
	                near[0] = far[0] = std::sqrt( std::fabs(g[0][0]) );
	
	                ginv[0][0] =  ",floatcast["1"],"/g[0][0];
	                
	                //  dfdagger = g^{-1} * df^T (",DomDim," x ",AmbDim," matrix)
",
Table[line[5,"dfdagger[",s[j],"][",s[k],"] = ",Riffle[Table[{"ginv[",s[Min[i,j]],"][",s[Max[i,j]],"] * df[",s[k],"][",s[i],"]"},{i,0,n-1}]," + "]],{j,0,n-1},{k,0,m-1}],
"\n",
Module[{c=m},
	Table[line[5,"near[",tostr[++c+n m],"] = far[",tostr[c],"] = ",If[i==j,floatcast["1.0"],"  "],Table[{" - df[",s[i],"][",s[k],"] * dfdagger[",s[k],"][",s[j],"]"},{k,0,n-1}]],{i,0,m-1},{j,i,m-1}]
],
"
	                // derivative operator  ("<>AmbDim<>" x "<>s[n+1]<>" matrix)
	
	                Real * Df = &DiffOp_value[ ",s[m(n+1)]," * i ];

",Table[{line[5,"Df["<>tostr[k * (n+1)]<>"] =",Table[{" - dfdagger[",s[j],"][",s[k],"]"},{j,0,n-1}]],Table[line[4,"Df["<>tostr[k * (n+1)+j+1]<>"] =   dfdagger[",s[j],"][",s[k],"]"],{j,0,n-1}]}
,
{k,0,m-1}],
"
	            }
			},
			ThreadCount()
		);

        ptoc(ClassName()+\"::",name,"\");
    }
"
]
];


ComputeNearFarDataOps[2,m_]:=Module[{name,DomDim,AmbDim,HullSize,neardim,fardim,numdigits,tostr,n,s},
n=2;
s=ToString;
name="ComputeNearFarDataOps";
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
		      Tensor2<Real,Int> & P_coords,
		      Tensor3<Real,Int> & P_hull_coords,
		      Tensor2<Real,Int> & P_near,
		      Tensor2<Real,Int> & P_far,
		SparseMatrix_T & DiffOp,
		SparseMatrix_T & AvOp 
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
				mptr<LInt> AvOp_outer = AvOp.Outer().data();
				mptr< Int> AvOp_inner = AvOp.Inner().data();
				mptr<Real> AvOp_value = AvOp.Values().data();
	
				mptr<LInt> DiffOp_outer = DiffOp.Outer().data();
				mptr< Int> DiffOp_inner = DiffOp.Inner().data();
				mptr<Real> DiffOp_value = DiffOp.Values().data();
	
				cptr<Real> V_coords__      = V_coords.data();
				
				cptr< Int> simplices__     = simplices.data();
			    mptr<Real> P_hull_coords__ = P_hull_coords.data();
				mptr<Real> P_coords__      = P_coords.data();
	
				Real df       [",AmbDim,"][",DomDim,"];
				Real dfdagger [",DomDim,"][",AmbDim,"];
				Real g        [",DomDim,"][",DomDim,"];
				Real ginv     [",DomDim,"][",DomDim,"];
	
				Int simplex        [",s[n+1],"];
	            Int sorted_simplex [",s[n+1],"];
	
				const Int i_begin = job_ptr[thread];
				const Int i_end   = job_ptr[thread+1];
	
	            for( Int i = i_begin; i < i_end; ++i )
	            {
	
					mptr<Real> near = P_near.data(i);                    
					mptr<Real> far  = P_far.data(i);

",
Table[line[5,"simplex[",s[j],"] = sorted_simplex[",s[j],"] = simplices__[",s[n+1],"*i +",s[j],"]"],{j,0,n+1-1}],
"                  
	                // sorting simplex so that we do not have to sort the sparse arrays to achieve CSR format later
	                Sort( sorted_simplex, sorted_simplex + ",s[n+1]," );
","
					AvOp_outer[i+1] = (i+1) * ",s[n+1],";                      
",
Table[line[5,"AvOp_inner[",s[n+1],"*i+",s[j],"] = sorted_simplex[",s[j],"]"],{j,0,n+1-1}],"\n",
Table[line[5,"AvOp_value[",s[n+1],"*i+",s[j],"] = ",s[1./(n+1),CForm]],{j,0,n+1-1}],"\n",
Table[line[5,"DiffOp_outer[",AmbDim,"*i+",s[k],"] = (",AmbDim," * i + ",s[k],") * ",s[n+1]]
,{k,0,m-1}],
"\n",
Table[line[5,"DiffOp_inner[(i * ",AmbDim," + ",s[k],") * ",s[n+1]," + ",s[j]," ] = sorted_simplex[",s[j],"]"],{k,0,m-1},{j,0,n+1-1}],
"\n",
Table[
line[5,"near[",s[1+m j+k],"] = P_hull_coords__[",HullSize,"*i+",s[m j+k],"] = V_coords__[",AmbDim,"*simplex[",s[j],"]+",s[k],"]"],{j,0,n+1-1},{k,0,m-1}],
"\n",
Table[
line[5,"far["<>s[k+1]<>"] = P_coords__[",AmbDim,"*i+",s[k],"] = ",s[1./(n+1),CForm]," * ( ",Riffle[Table[{"P_hull_coords__[",HullSize,"*i+",s[m j+k],"]"},{j,0,n+1-1}]," + "]," )"],{k,0,m-1}],
"\n",
Table[line[5,"df["<>s[k],"][",s[j]<>"] = V_coords__[",AmbDim,"*sorted_simplex[",s[j+1],"]+",s[k],"] - V_coords__[",AmbDim,"*sorted_simplex[0]+",s[k],"]"],{k,0,m-1},{j,0,n-1}],
"\n",
Table[line[5,"g[",s[i],"][",s[j],"] = ",Riffle[Table[{"df[",s[k],"][",s[i]<>"] * df[",s[k],"][",s[j],"]"},{k,0,m-1}]," + "]],{i,0,n-1},{j,0,n-1}],
"
	                Real det = g[0][0] * g[1][1] - g[0][1] * g[0][1];
	
	                near[0] = far[0] = std::sqrt( std::fabs(det) ) * ",floatcast["0.5"],";
	
	                Real invdet = ",floatcast["1.0"],"/det;
	                ginv[0][0] =  g[1][1] * invdet;
	                ginv[0][1] = -g[0][1] * invdet;
	                ginv[1][1] =  g[0][0] * invdet;
	                
	                //  dfdagger = g^{-1} * df^T ("<>DomDim<>" x "<>AmbDim<>" matrix)
",  
Table[line[5,"dfdagger[",s[j],"][",s[k],"] = ",Riffle[Table[{"ginv[",s[Min[i,j]],"][",s[Max[i,j]],"] * df[",s[k],"][",s[i],"]"},{i,0,n-1}]," + "]],{j,0,n-1},{k,0,m-1}],"            
",
Module[{c=m},
	Table[line[5,"near[",tostr[++c+n m],"] = far[",tostr[c],"]  = ",If[i==j,floatcast["1.0"],"  "],Table[{" - df[",s[i],"][",s[k],"] * dfdagger[",s[k],"][",s[j],"]"},{k,0,n-1}]],{i,0,m-1},{j,i,m-1}]
],
"
	                // derivative operator  ("<>AmbDim<>" x "<>s[n+1]<>" matrix)
	
	                mptr<Real> Df = &DiffOp_value[ ",s[m(n+1)]," * i ];

",
StringJoin[
	Table[{
		line[5,"Df["<>tostr[k * (n+1)]<>"] =",Table[{" - dfdagger[",s[j],"][",s[k],"]"},{j,0,n-1}]],
		Table[line[5,"Df["<>tostr[k * (n+1)+j+1]<>"] =   dfdagger[",s[j],"][",s[k],"]"],{j,0,n-1}]
		},{k,0,m-1}]],
"
	            }
			},
			ThreadCount()
		);

        ptoc(ClassName()+\"::",name,"\");
    }
"]
];
