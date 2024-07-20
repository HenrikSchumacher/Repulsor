(* ::Package:: *)

ClearAll[DNearToHulls];


DNearToHulls[n_,m_]:=Module[{name},
name="DNearToHulls";
StringJoin[
"
	template<AddTo_T addtoQ> 
	void ",name,"( 
		cref<Tensor2<Real,Int>> V_coords,
        cref<Tensor2<Int ,Int>> simplices,
        cref<Tensor2<Real,Int>> P_D_near,
        cref<Tensor1<Real,Int>> V_charges,
        // cppcheck-suppress [constParameter]
        mref<Tensor3<Real,Int>> buffer
	) const
    {
        ptic(ClassName()+\"::",name,"\");
        eprint(ClassName()+\"::",name," not implemented. Returning 0.\");
		
		if constexpr ( addtoQ != AddTo )
		{
			buffer.Fill(static_cast<Real>(0));
		}
        ptoc(ClassName()+\"::",name,"\");
    }
"]
];


DNearToHulls[0,m_/;m>=0]:=Module[{name,DomDim,AmbDim,HullSize,numdigits,tostr,n,s,
datadim,NearDim,PP,P,pr,\[Tau],vol,dXX,A,PDdata,vcpdata,code0,scratchsymbols,scratchrules,$,code1,codestringAssign,codestringAddTo
},
n=0;
s=ToString;
name="DNearToHulls";
DomDim=ToString[n];
AmbDim=ToString[m];
HullSize=ToString[m(n+1)];
datadim=1+(n+1)m+m (m+1)/2;
NearDim=ToString[datadim];

numdigits=Max[Length[IntegerDigits[m(n+1),10]],Length[IntegerDigits[datadim,10]]];
tostr[i_]:=ConstantArray[" ",numdigits-Length[IntegerDigits[i,10]]]<>ToString[i];

PP=Table[Indexed[P,{i,j}],{i,1,n+1},{j,1,m}];
pr=IdentityMatrix[m];
vol=1;

dXX=Table[Indexed[PDdata,i],{i,1,datadim}];
vcpdata=Join[{vol},vol Flatten[PP],vol Flatten[Table[pr[[i,j]],{i,1,m},{j,i,m}]]];
A= "charge" (dXX . D[vcpdata,{Flatten[PP],1}]);

code0=Part[Experimental`OptimizeExpression[A,"OptimizationLevel"->2,
"OptimizationSymbol"->$
]/.CompoundExpression->List/.Set-> CAssign,
-1]/.Thread[
	Flatten[PP]->Flatten[Table["V_coords__["<>AmbDim<>"*simplices__["<>s[n+1]<>"*i+"<>s[j]<>"]+"<>s[k]<>"]",{j,0,n},{k,0,m-1}]]
]/.Thread[
	Flatten[dXX]->Flatten[Table["P_D_near__["<>NearDim<>"*i+"<>s[j]<>"]",{j,0,datadim-1}]]
]/.Power[x_,2]:>HoldForm[x x];

scratchsymbols=Cases[code0,CAssign[symbol_,rhs_]:>symbol,\[Infinity]];
scratchrules=MapIndexed[#1->"s"<>ToString[#2[[1]]-1]&,scratchsymbols];


code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];

code1=Flatten@MapIndexed[CAssign["buffer__["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1,m],{2}];

codestringAssign=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->5],"Sqrt("->"std::sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];

code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];

code1=Flatten@MapIndexed[CAddTo["buffer__["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1,m],{2}];

codestringAddTo=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->5],"Sqrt("->"std::sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];


StringJoin["
	template<AddTo_T addtoQ> 
	void ",name,"( 
        cref<Tensor2<Real,Int>> V_coords,
        cref<Tensor2<Int ,Int>> simplices,
        cref<Tensor2<Real,Int>> P_D_near,
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

		//cptr<Real> V_coords__  = V_coords.data();
		//cptr<Int>  simplices__ = simplices.data();
		cptr<Real> P_D_near__  = P_D_near.data();
		mptr<Real> buffer__    = buffer.data();
        
        if constexpr ( addtoQ == AddTo )
		{
			ParallelDo(
				[=]( const Int i )
				{
					Real charge = 0;
					for( Int k = 0; k < SIZE; ++k )
					{
						charge += V_charges[simplices(i,k)];
					}
					charge *= nth;

"<>codestringAddTo<>"
				},
				simplices.Dimension(0),
				ThreadCount()
			);
		}
		else
		{
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
		}

        ptoc(ClassName()+\"::"<>name<>"\");
        
    }
"]
]


DNearToHulls[1,m_/;m>=0]:=Module[{name,DomDim,AmbDim,HullSize,numdigits,tostr,n,s,
datadim,NearDim,PP,P,pr,\[Tau],vol,dXX,A,PDdata,vcpdata,code0,scratchsymbols,scratchrules,$,code1,codestringAssign,codestringAddTo
},
n=1;
s=ToString;
name="DNearToHulls";
DomDim=ToString[n];
AmbDim=ToString[m];
HullSize=ToString[m(n+1)];
datadim=1+(n+1)m+m (m+1)/2;
NearDim=ToString[datadim];

numdigits=Max[Length[IntegerDigits[m(n+1),10]],Length[IntegerDigits[datadim,10]]];
tostr[i_]:=ConstantArray[" ",numdigits-Length[IntegerDigits[i,10]]]<>ToString[i];

PP=Table[Indexed[P,{i,j}],{i,1,n+1},{j,1,m}];

\[Tau]=PP[[2]]-PP[[1]];
pr=IdentityMatrix[m]-KroneckerProduct[\[Tau],\[Tau]]/Dot[\[Tau] . \[Tau]];
vol=Sqrt[# . #]&[\[Tau]];

dXX=Table[Indexed[PDdata,i],{i,1,datadim}];
vcpdata=Join[{vol},vol Flatten[PP],vol Flatten[Table[pr[[i,j]],{i,1,m},{j,i,m}]]];
A= "charge" (dXX . D[vcpdata,{Flatten[PP],1}]);

code0=Part[Experimental`OptimizeExpression[A,"OptimizationLevel"->2,
"OptimizationSymbol"->$
]/.CompoundExpression->List/.Set-> CAssign,
-1,-1]/.Thread[
	Flatten[PP]->Flatten[Table["V_coords__["<>AmbDim<>"*simplices__["<>s[n+1]<>"*i+"<>s[j]<>"]+"<>s[k]<>"]",{j,0,n},{k,0,m-1}]]
]/.Thread[
	Flatten[dXX]->Flatten[Table["P_D_near__["<>NearDim<>"*i+"<>s[j]<>"]",{j,0,datadim-1}]]
]/.Power[x_,2]:>HoldForm[x x];

scratchsymbols=Cases[code0,CAssign[symbol_,rhs_]:>symbol,\[Infinity]];
scratchrules=MapIndexed[#1->"s"<>ToString[#2[[1]]-1]&,scratchsymbols];


code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];

code1[[-1]]=Flatten@MapIndexed[CAssign["buffer__["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1[[-1]],m],{2}];

codestringAssign=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->5],"Sqrt("->"std::sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];

code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];

code1[[-1]]=Flatten@MapIndexed[CAddTo["buffer__["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1[[-1]],m],{2}];

codestringAddTo=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->5],"Sqrt("->"std::sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];


StringJoin["
	template<AddTo_T addtoQ> 
	void ",name,"( 
        cref<Tensor2<Real,Int>> V_coords,
        cref<Tensor2<Int ,Int>> simplices,
        cref<Tensor2<Real,Int>> P_D_near,
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

		cptr<Real> V_coords__  = V_coords.data();
		cptr<Int>  simplices__ = simplices.data();
		cptr<Real> P_D_near__  = P_D_near.data();
		mptr<Real> buffer__    = buffer.data();
        
        if constexpr ( addtoQ == AddTo )
		{
			ParallelDo(
				[=]( const Int i )
				{
					Real charge = 0;
					for( Int k = 0; k < SIZE; ++k )
					{
						charge += V_charges[simplices(i,k)];
					}
					charge *= nth;

"<>codestringAddTo<>"
				},
				simplices.Dimension(0),
				ThreadCount()
			);
		}
		else
		{
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
		}

        ptoc(ClassName()+\"::"<>name<>"\");
        
    }
"]
];


DNearToHulls[n_/;n>=2,m_/;m>=0]:=Module[{name,DomDim,AmbDim,HullSize,numdigits,tostr,datadim,NearDim,\[Nu],g,df,s,
PP,P,c,pr,vol,dXX,A,PDdata,vcpdata,code0,scratchsymbols,scratchrules,$,code1,codestringAssign,codestringAddTo
},
s=ToString;
name="DNearToHulls";
DomDim=s[n];
AmbDim=s[m];
HullSize=s[m(n+1)];
datadim=1+(n+1)m+m (m+1)/2;
NearDim=s[datadim];
numdigits=Max[Length[IntegerDigits[m(n+1),10]],Length[IntegerDigits[datadim,10]]];
tostr[i_]:=ConstantArray[" ",numdigits-Length[IntegerDigits[i,10]]]<>s[i];

PP=Table[Indexed[P,{i,j}],{i,1,n+1},{j,1,m}];
If[m==3,
(
	\[Nu]=Cross[PP[[2]]-PP[[1]],PP[[3]]-PP[[1]]];
	pr=KroneckerProduct[\[Nu],\[Nu]]/Dot[\[Nu] . \[Nu]];
	vol=1/2 Sqrt[ # . #]&[\[Nu]];

	dXX=Table[Indexed[PDdata,i],{i,1,datadim}];
	vcpdata=Join[{vol},vol Flatten[PP],vol Flatten[Table[pr[[i,j]],{i,1,m},{j,i,m}]]];
)
,
(
	df=Transpose[Table[PP[[j+1]]-PP[[1]],{j,1,n}]];
	g=df\[Transpose] . df;
	vol= "charge" Sqrt[Det[g]]/n!;
	pr=Flatten[Table[(df . Inverse[g] . Transpose[df])[[i,j]],{i,1,m},{j,i,m}]];

	dXX=Table[Indexed[PDdata,i],{i,1,datadim}];
	vcpdata=Join[{vol},vol Flatten[PP],vol pr];
)
];

A= "charge" (dXX . D[vcpdata,{Flatten[PP],1}]);

code0=Part[
Experimental`OptimizeExpression[A,
"OptimizationLevel"->2,
"OptimizationSymbol"->$
]/.CompoundExpression->List/.Set-> CAssign,
-1,-1]/.Thread[
Flatten[PP]->
Flatten[Table["V_coords__["<>AmbDim<>"*simplices__["<>s[n+1]<>"*i+"<>s[j]<>"]+"<>s[k]<>"]",
{j,0,n},{k,0,m-1}]]
]/.Thread[
Flatten[dXX]->Flatten[Table["P_D_near__["<>NearDim<>"*i+"<>s[j]<>"]",{j,0,datadim-1}]]
]/.Power[x_,2]:>HoldForm[x x];

scratchsymbols=Cases[code0,CAssign[symbol_,rhs_]:>symbol,\[Infinity]];
scratchrules=MapIndexed[#1->"s"<>ToString[#2[[1]]-1]&,scratchsymbols];

code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];

code1[[-1]]=Flatten@MapIndexed[CAssign["buffer__["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1[[-1]],m],{2}];

codestringAssign=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->5],"Sqrt("->"std::sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];

code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];

code1[[-1]]=Flatten@MapIndexed[CAddTo["buffer__["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1[[-1]],m],{2}];

codestringAddTo=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->5],"Sqrt("->"std::sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];

StringJoin["
	template<AddTo_T addtoQ> 
    void ",name,"( 
        cref<Tensor2<Real,Int>> V_coords,
        cref<Tensor2<Int ,Int>> simplices,
        cref<Tensor2<Real,Int>> P_D_near,
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
        
		cptr<Real> V_coords__  = V_coords.data();
        cptr<Int>  simplices__ = simplices.data();
		cptr<Real> P_D_near__  = P_D_near.data();
        mptr<Real> buffer__    = buffer.data();

		if constexpr ( addtoQ == AddTo )
		{
			ParallelDo(
				[=]( const Int i )
				{
					Real charge = 0;
					for( Int k = 0; k < SIZE; ++k )
					{
						charge += V_charges[simplices(i,k)];
					}
					charge *= nth;

"<>codestringAddTo<>"
				},
				simplices.Dimension(0),
				ThreadCount()
			);
		}
		else
		{
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
		}

        ptoc(ClassName()+\"::"<>name<>"\");
        
    }
"]
];
