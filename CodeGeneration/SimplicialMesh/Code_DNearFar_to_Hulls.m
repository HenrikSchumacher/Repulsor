(* ::Package:: *)

ClearAll[DNearFarToHulls];


DNearFarToHulls[n_,m_]:=Module[{name},
name="DNearFarToHulls";
StringJoin[
"
	void ",name,"( 
		cref<Tensor2<Real,Int>> V_coords, 
		cref<Tensor2<Int ,Int>> simplices, 
		cref<Tensor2<Real,Int>> P_D_near, 
		cref<Tensor2<Real,Int>> P_D_far, 
		cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
		mref<Tensor3<Real,Int>> buffer
	) const
    {
        ptic(ClassName()+\"::",name,"\");
        eprint(ClassName()+\"::",name," not implemented. Returning 0.\");
		
		buffer.Fill(static_cast<Real>(0));

        ptoc(ClassName()+\"::",name,"\");
    }
"]
];


DNearFarToHulls[0,m_/;m>=0]:=Module[{name,DomDim,AmbDim,HullSize,numdigits,tostr,n,s,
neardim,fardim,NearDim,FarDim,PP,P,pr,\[Tau],vol,dXXnear,dXXfar,A,PDnear,PDfar,vcpneardata,vcpfardata,code0,scratchsymbols,scratchrules,$,code1,codestringAssign
},
n=0;
s=ToString;
name="DNearFarToHulls";
DomDim=ToString[n];
AmbDim=ToString[m];
HullSize=ToString[m(n+1)];
neardim=1+(n+1)m+m (m+1)/2;
fardim=1+m+m (m+1)/2;

NearDim=ToString[neardim];
FarDim=ToString[fardim];

numdigits=Max[Length[IntegerDigits[m(n+1),10]],Length[IntegerDigits[fardim,10]]+Length[IntegerDigits[neardim,10]]];
tostr[i_]:=ConstantArray[" ",numdigits-Length[IntegerDigits[i,10]]]<>ToString[i];

PP=Table[Indexed[P,{i,j}],{i,1,n+1},{j,1,m}];
pr=IdentityMatrix[m];
vol=1;

dXXnear=Table[Indexed[PDnear,i],{i,1,neardim}];
dXXfar=Table[Indexed[PDfar,i],{i,1,fardim}];
vcpneardata=Join[{vol},vol Flatten[PP],vol Flatten[Table[pr[[i,j]],{i,1,m},{j,i,m}]]];
vcpfardata=Join[{vol},vol Mean[PP],vol Flatten[Table[pr[[i,j]],{i,1,m},{j,i,m}]]];
A= "charge" (dXXnear . D[vcpneardata,{Flatten[PP],1}] + dXXfar . D[vcpfardata,{Flatten[PP],1}]);

code0=Part[Experimental`OptimizeExpression[A,"OptimizationLevel"->2,
"OptimizationSymbol"->$
]/.CompoundExpression->List/.Set-> CAssign,
-1]/.Thread[
	Flatten[PP]->Flatten[Table["X["<>AmbDim<>"*S["<>s[n+1]<>"*i+"<>s[j]<>"]+"<>s[k]<>"]",{j,0,n},{k,0,m-1}]]
]/.Thread[
	Flatten[dXXnear]->Flatten[Table["N["<>NearDim<>"*i+"<>s[j]<>"]",{j,0,neardim-1}]]
]/.Thread[
	Flatten[dXXfar]->Flatten[Table["F["<>FarDim<>"*i+"<>s[j]<>"]",{j,0,fardim-1}]]
]/.Power[x_,2]:>HoldForm[x x];

scratchsymbols=Cases[code0,CAssign[symbol_,rhs_]:>symbol,\[Infinity]];
scratchrules=MapIndexed[#1->"s"<>ToString[#2[[1]]-1]&,scratchsymbols];


code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];

code1=Flatten@MapIndexed[CAssign["B["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1,m],{2}];

codestringAssign=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->4],"Sqrt("->"std::sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];

StringJoin["
	void ",name,"( 
		cref<Tensor2<Real,Int>> V_coords, 
		cref<Tensor2<Int ,Int>> simplices, 
		cref<Tensor2<Real,Int>> P_D_near,
		cref<Tensor2<Real,Int>> P_D_far, 
		cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
		mref<Tensor3<Real,Int>> buffer
	) const
    {
        ptic(ClassName()+\"::"<>name<>"\");

        if( P_D_near.Dimension(1) != "<>NearDim<>" )
        {
            eprint(\"in "<>name<>": P_D_near.Dimension(1) != "<>NearDim<>". Aborting\");
        }

		//cptr<Real> X = V_coords.data();
		//cptr<Int>  S = simplices.data();
		cptr<Real> N = P_D_near.data();
		cptr<Real> F = P_D_far.data();
		mptr<Real> B = buffer.data();
        
		ParallelDo(
			[=]( const Int i )
			{
				Real charge = 0;
				for( Int k = 0; k < SIZE; ++k )
				{
					charge += V_charges[simplices(i,k)];
				}
				charge *= nth;
"<>codestringAssign<>"
			},
			simplices.Dimension(0),
			ThreadCount()
		);

        ptoc(ClassName()+\"::"<>name<>"\");
        
    }
"]
]


DNearFarToHulls[1,m_/;m>=0]:=Module[{name,DomDim,AmbDim,HullSize,numdigits,tostr,n,s,
neardim,fardim,NearDim,FarDim,PP,P,pr,\[Tau],vol,dXXnear,dXXfar,A,PDnear,PDfar,vcpneardata,vcpfardata,code0,scratchsymbols,scratchrules,$,code1,codestringAssign
},
n=1;
s=ToString;
name="DNearFarToHulls";
DomDim=ToString[n];
AmbDim=ToString[m];
HullSize=ToString[m(n+1)];
neardim=1+(n+1)m+m (m+1)/2;
fardim=1+m+m (m+1)/2;
NearDim=ToString[neardim];
FarDim=ToString[fardim];

numdigits=Max[Length[IntegerDigits[m(n+1),10]],Length[IntegerDigits[fardim,10]]+Length[IntegerDigits[neardim,10]]];
tostr[i_]:=ConstantArray[" ",numdigits-Length[IntegerDigits[i,10]]]<>ToString[i];

PP=Table[Indexed[P,{i,j}],{i,1,n+1},{j,1,m}];

\[Tau]=PP[[2]]-PP[[1]];
pr=IdentityMatrix[m]-KroneckerProduct[\[Tau],\[Tau]]/Dot[\[Tau] . \[Tau]];
vol=Sqrt[# . #]&[\[Tau]];

dXXnear=Table[Indexed[PDnear,i],{i,1,neardim}];
dXXfar=Table[Indexed[PDfar,i],{i,1,fardim}];
vcpneardata=Join[{vol},vol Flatten[PP],vol Flatten[Table[pr[[i,j]],{i,1,m},{j,i,m}]]];
vcpfardata=Join[{vol},vol Mean[PP],vol Flatten[Table[pr[[i,j]],{i,1,m},{j,i,m}]]];
A= "charge" (dXXnear . D[vcpneardata,{Flatten[PP],1}] + dXXfar . D[vcpfardata,{Flatten[PP],1}]);

code0=Part[Experimental`OptimizeExpression[A,"OptimizationLevel"->2,
"OptimizationSymbol"->$
]/.CompoundExpression->List/.Set-> CAssign,
-1,-1]/.Thread[
	Flatten[PP]->Flatten[Table["X["<>AmbDim<>"*S["<>s[n+1]<>"*i+"<>s[j]<>"]+"<>s[k]<>"]",{j,0,n},{k,0,m-1}]]
]/.Thread[
	Flatten[dXXnear]->Flatten[Table["N["<>NearDim<>"*i+"<>s[j]<>"]",{j,0,neardim-1}]]
]/.Thread[
	Flatten[dXXfar]->Flatten[Table["F["<>FarDim<>"*i+"<>s[j]<>"]",{j,0,fardim-1}]]
]/.Power[x_,2]:>HoldForm[x x];

scratchsymbols=Cases[code0,CAssign[symbol_,rhs_]:>symbol,\[Infinity]];
scratchrules=MapIndexed[#1->"s"<>ToString[#2[[1]]-1]&,scratchsymbols];


code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];

code1[[-1]]=Flatten@MapIndexed[CAssign["B["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1[[-1]],m],{2}];

codestringAssign=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->4],"Sqrt("->"std::sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];

StringJoin["
	void ",name,"( 
		cref<Tensor2<Real,Int>> V_coords, 
		cref<Tensor2<Int ,Int>> simplices, 
		cref<Tensor2<Real,Int>> P_D_near, 
		cref<Tensor2<Real,Int>> P_D_far, 
		cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
		mref<Tensor3<Real,Int>> buffer
	) const
    {
        ptic(ClassName()+\"::"<>name<>"\");

        if( P_D_near.Dimension(1) != "<>NearDim<>" )
        {
            eprint(\"in "<>name<>": P_D_near.Dimension(1) != "<>NearDim<>". Aborting\");
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
					charge += V_charges[simplices(i,k)];
				}
				charge *= nth;
"<>codestringAssign<>"
			},
			simplices.Dimension(0),
			ThreadCount()
		);

        ptoc(ClassName()+\"::"<>name<>"\");
        
    }
"]
];


DNearFarToHulls[n_/;n>=2,m_/;m>=0]:=Module[{name,DomDim,AmbDim,HullSize,numdigits,tostr,neardim,fardim,NearDim,FarDim,\[Nu],g,df,s,
PP,P,c,pr,vol,dXXnear,dXXfar,A,PDnear,PDfar,vcpneardata,vcpfardata,code0,scratchsymbols,scratchrules,$,code1,codestringAssign
},
s=ToString;
name="DNearFarToHulls";
DomDim=s[n];
AmbDim=s[m];
HullSize=s[m(n+1)];
neardim=1+(n+1)m+m (m+1)/2;
fardim=1+m+m (m+1)/2;
NearDim=ToString[neardim];
FarDim=ToString[fardim];

numdigits=Max[Length[IntegerDigits[m(n+1),10]],Length[IntegerDigits[fardim,10]]+Length[IntegerDigits[neardim,10]]];
tostr[i_]:=ConstantArray[" ",numdigits-Length[IntegerDigits[i,10]]]<>s[i];

PP=Table[Indexed[P,{i,j}],{i,1,n+1},{j,1,m}];

dXXnear=Table[Indexed[PDnear,i],{i,1,neardim}];
dXXfar=Table[Indexed[PDfar,i],{i,1,fardim}];
	
If[m==3,
(
	\[Nu]=Cross[PP[[2]]-PP[[1]],PP[[3]]-PP[[1]]];
	pr=KroneckerProduct[\[Nu],\[Nu]]/Dot[\[Nu] . \[Nu]];
	vol=1/2 Sqrt[ # . #]&[\[Nu]];
)
,
(
	df=Transpose[Table[PP[[j+1]]-PP[[1]],{j,1,n}]];
	g=df\[Transpose] . df;
	vol= "charge" Sqrt[Det[g]]/n!;
	pr=df . Inverse[g] . Transpose[df];
)
];

vcpneardata=Join[{vol},vol Flatten[PP],vol Flatten[Table[pr[[i,j]],{i,1,m},{j,i,m}]]];
vcpfardata=Join[{vol},vol Mean[PP],vol Flatten[Table[pr[[i,j]],{i,1,m},{j,i,m}]]];
A= "charge" (dXXnear . D[vcpneardata,{Flatten[PP],1}] + dXXfar . D[vcpfardata,{Flatten[PP],1}]);
code0=Part[
Experimental`OptimizeExpression[A,
"OptimizationLevel"->2,
"OptimizationSymbol"->$
]/.CompoundExpression->List/.Set-> CAssign,
-1,-1]/.Thread[
Flatten[PP]->
Flatten[Table["X["<>AmbDim<>"*S["<>s[n+1]<>"*i+"<>s[j]<>"]+"<>s[k]<>"]", {j,0,n},{k,0,m-1}]]
]/.Thread[
	Flatten[dXXnear]->Flatten[Table["N["<>NearDim<>"*i+"<>s[j]<>"]",{j,0,neardim-1}]]
]/.Thread[
	Flatten[dXXfar]->Flatten[Table["F["<>FarDim<>"*i+"<>s[j]<>"]",{j,0,fardim-1}]]
]/.Power[x_,2]:>HoldForm[x x];

scratchsymbols=Cases[code0,CAssign[symbol_,rhs_]:>symbol,\[Infinity]];
scratchrules=MapIndexed[#1->"s"<>ToString[#2[[1]]-1]&,scratchsymbols];

code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];

code1[[-1]]=Flatten@MapIndexed[CAssign["B["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1[[-1]],m],{2}];

codestringAssign=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->4],"Sqrt("->"std::sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];

StringJoin["
    void ",name,"( 
		cref<Tensor2<Real,Int>> V_coords, 
		cref<Tensor2<Int ,Int>> simplices, 
		cref<Tensor2<Real,Int>> P_D_near, 
		cref<Tensor2<Real,Int>> P_D_far, 
		cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
		mref<Tensor3<Real,Int>> buffer
	) const
    {
        ptic(ClassName()+\"::"<>name<>"\");

        if( P_D_near.Dimension(1) != "<>NearDim<>" )
        {
            eprint(\"in "<>name<>": P_D_near.Dimension(1) != "<>NearDim<>". Aborting\");
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
					charge += V_charges[simplices(i,k)];
				}
				charge *= nth;

"<>codestringAssign<>"
			},
			simplices.Dimension(0),
			ThreadCount()
		);

        ptoc(ClassName()+\"::"<>name<>"\");
        
    }
"]
];
