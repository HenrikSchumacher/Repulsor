(* ::Package:: *)

ClearAll[DFarToHulls];


DFarToHulls[n_,m_]:=Module[{name},
name="DFarToHulls";
StringJoin[
"
	void ",name,"( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_far,
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+\"::",name,"\");
                eprint(ClassName()+\"::",name," not implemented. Returning 0.\");
		
		if(!addTo)
		{
			buffer.Fill(static_cast<Real>(0));
		}
        ptoc(ClassName()+\"::",name,"\");
    }
"]
];


DFarToHulls[0,m_/;m>0]:=Module[{name,DomDim,AmbDim,HullSize,numdigits,tostr,n,s,
datadim,FarDim,PP,P,pr,\[Tau],vol,dXX,A,PDdata,vcpdata,code0,scratchsymbols,scratchrules,$,code1,codestringAssign,codestringAddTo
},
n=0;
s=ToString;
name="DFarToHulls";
DomDim=ToString[n];
AmbDim=ToString[m];
HullSize=ToString[m(n+1)];
datadim=1+m+m (m+1)/2;
FarDim=ToString[datadim];

numdigits=Max[Length[IntegerDigits[m(n+1),10]],Length[IntegerDigits[datadim,10]]];
tostr[i_]:=ConstantArray[" ",numdigits-Length[IntegerDigits[i,10]]]<>ToString[i];


PP=Table[Indexed[P,{i,j}],{i,1,n+1},{j,1,m}];
pr=IdentityMatrix[m];
vol=1;

dXX=Table[Indexed[PDdata,i],{i,1,datadim}];
vcpdata=Join[{vol},vol Mean[PP],vol Flatten[Table[pr[[i,j]],{i,1,m},{j,i,m}]]];
A=dXX . D[vcpdata,{Flatten[PP],1}];

code0=Part[
Experimental`OptimizeExpression[A,"OptimizationLevel"->2,"OptimizationSymbol"->$],
-1
]/.CompoundExpression->List/.Set-> CAssign/.Thread[
	Flatten[PP]->Flatten[Table["V_coords__["<>AmbDim<>"*simplices__["<>s[n+1]<>"*i+"<>s[j]<>"]+"<>s[k]<>"]",{j,0,n},{k,0,m-1}]]
]/.Thread[
	Flatten[dXX]->Flatten[Table["P_D_far__["<>FarDim<>"*i+"<>s[j]<>"]",{j,0,datadim-1}]]
]/.Power[x_,2]:>HoldForm[x x];

scratchsymbols=Cases[code0,CAssign[symbol_,rhs_]:>symbol,\[Infinity]];
scratchrules=MapIndexed[#1->"s"<>ToString[#2[[1]]-1]&,scratchsymbols];


code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];


code1=Flatten@MapIndexed[CAssign["buffer__["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1,m],{2}];

codestringAssign=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->5],"Sqrt("->"sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];

code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];

code1=Flatten@MapIndexed[CAddTo["buffer__["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1,m],{2}];

codestringAddTo=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->5],"Sqrt("->"sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];

StringJoin["
	void ",name,"( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_far, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+\"::"<>name<>"\");

        if( P_D_far.Dimension(1) != "<>FarDim<>" )
        {
            eprint(\"in "<>name<>": P_D_far.Dimension(1) != "<>FarDim<>". Aborting\");
        }

		//ptr<Real> V_coords__  = V_coords.data();
		//ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_far__   = P_D_far.data();
		mut<Real> buffer__    = buffer.data();
        
        if( addTo )
		{
			ParallelDo(
				[=]( const Int i )
				{
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


DFarToHulls[1,m_/;m>0]:=Module[{name,DomDim,AmbDim,HullSize,numdigits,tostr,n,s,
datadim,FarDim,PP,P,pr,\[Tau],vol,dXX,A,PDdata,vcpdata,code0,scratchsymbols,scratchrules,$,code1,codestringAssign,codestringAddTo
},
n=1;
s=ToString;
name="DFarToHulls";
DomDim=ToString[n];
AmbDim=ToString[m];
HullSize=ToString[m(n+1)];
datadim=1+m+m (m+1)/2;
FarDim=ToString[datadim];

numdigits=Max[Length[IntegerDigits[m(n+1),10]],Length[IntegerDigits[datadim,10]]];
tostr[i_]:=ConstantArray[" ",numdigits-Length[IntegerDigits[i,10]]]<>ToString[i];


PP=Table[Indexed[P,{i,j}],{i,1,n+1},{j,1,m}];
\[Tau]=PP[[2]]-PP[[1]];
pr=IdentityMatrix[m]-KroneckerProduct[\[Tau],\[Tau]]/Dot[\[Tau] . \[Tau]];
vol=Sqrt[# . #]&[\[Tau]];

dXX=Table[Indexed[PDdata,i],{i,1,datadim}];
vcpdata=Join[{vol},vol Mean[PP],vol Flatten[Table[pr[[i,j]],{i,1,m},{j,i,m}]]];
A=dXX . D[vcpdata,{Flatten[PP],1}];

code0=Part[Experimental`OptimizeExpression[A,"OptimizationLevel"->2,
"OptimizationSymbol"->$
]/.CompoundExpression->List/.Set-> CAssign,
-1,-1]/.Thread[
	Flatten[PP]->Flatten[Table["V_coords__["<>AmbDim<>"*simplices__["<>s[n+1]<>"*i+"<>s[j]<>"]+"<>s[k]<>"]",{j,0,n},{k,0,m-1}]]
]/.Thread[
	Flatten[dXX]->Flatten[Table["P_D_far__["<>FarDim<>"*i+"<>s[j]<>"]",{j,0,datadim-1}]]
]/.Power[x_,2]:>HoldForm[x x];

scratchsymbols=Cases[code0,CAssign[symbol_,rhs_]:>symbol,\[Infinity]];
scratchrules=MapIndexed[#1->"s"<>ToString[#2[[1]]-1]&,scratchsymbols];


code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];

code1[[-1]]=Flatten@MapIndexed[CAssign["buffer__["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1[[-1]],m],{2}];

codestringAssign=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->5],"Sqrt("->"sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];

code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];

code1[[-1]]=Flatten@MapIndexed[CAddTo["buffer__["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1[[-1]],m],{2}];

codestringAddTo=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->5],"Sqrt("->"sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];

StringJoin["
	void ",name,"( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_far, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+\"::"<>name<>"\");

        if( P_D_far.Dimension(1) != "<>FarDim<>" )
        {
            eprint(\"in "<>name<>": P_D_far.Dimension(1) != "<>FarDim<>". Aborting\");
        }

		ptr<Real> V_coords__  = V_coords.data();
		ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_far__   = P_D_far.data();
		mut<Real> buffer__    = buffer.data();
        
        if( addTo )
		{
			ParallelDo(
				[=]( const Int i )
				{
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


DFarToHulls[n_/;n>=2,m_/;m>0]:=Module[{name,DomDim,AmbDim,HullSize,numdigits,tostr,datadim,FarDim,\[Nu],g,df,s,
PP,P,c,pr,vol,dXX,A,PDdata,vcpdata,code0,scratchsymbols,scratchrules,$,code1,codestringAssign,codestringAddTo
},
s=ToString;
name="DFarToHulls";
DomDim=s[n];
AmbDim=s[m];
HullSize=s[m(n+1)];
datadim=1+m+m (m+1)/2;
FarDim=s[datadim];
numdigits=Max[Length[IntegerDigits[m(n+1),10]],Length[IntegerDigits[datadim,10]]];
tostr[i_]:=ConstantArray[" ",numdigits-Length[IntegerDigits[i,10]]]<>s[i];

PP=Table[Indexed[P,{i,j}],{i,1,n+1},{j,1,m}];

If[m==3,
(
	\[Nu]=Cross[PP[[2]]-PP[[1]],PP[[3]]-PP[[1]]];
	pr=KroneckerProduct[\[Nu],\[Nu]]/Dot[\[Nu] . \[Nu]];
	vol=1/2Sqrt[ # . #]&[\[Nu]];

	dXX=Table[Indexed[PDdata,i],{i,1,datadim}];
	vcpdata=Join[{vol},vol Mean[PP],vol Flatten[Table[pr[[i,j]],{i,1,m},{j,i,m}]]];
)
,
(
	df=Transpose[Table[PP[[j+1]]-PP[[1]],{j,1,n}]];
	g=df\[Transpose] . df;
	vol=Sqrt[Det[g]]/n!;
	pr=Flatten[Table[(df . Inverse[g] . Transpose[df])[[i,j]],{i,1,m},{j,i,m}]];

	dXX=Table[Indexed[PDdata,i],{i,1,datadim}];
	vcpdata=Join[{vol},vol Mean[PP],vol pr];
)
];

A=dXX . D[vcpdata,{Flatten[PP],1}];

toC[a_String]:=a;
toC[a_]:=ToString[a,CForm];
arule={
CAddTo[a_,b_]:>CAddTo[toC[a],toC[b]],
CAssign[a_,b_]:>CAssign[toC[a],toC[b]],
CDeclareAssign["const Real",a_,b_]:>CDeclareAssign["const Real",toC[a],toC[b]]
};

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
Flatten[dXX]->Flatten[Table["P_D_far__["<>FarDim<>"*i+"<>s[j]<>"]",{j,0,datadim-1}]]
]/.Power[x_,2]:>HoldForm[x x];

scratchsymbols=Cases[code0,CAssign[symbol_,rhs_]:>symbol,\[Infinity]];
scratchrules=MapIndexed[#1->"s"<>ToString[#2[[1]]-1]&,scratchsymbols];

code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];

code1[[-1]]=Flatten@MapIndexed[CAssign["buffer__["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1[[-1]],m],{2}];

codestringAssign=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->5],"Sqrt("->"sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];

code1=code0/.scratchrules/.CAssign[symbol_,rhs_]:>CDeclareAssign["const Real",symbol,myCForm[rhs]];

code1[[-1]]=Flatten@MapIndexed[CAddTo["buffer__["<>HullSize<>"*i+"<>s[m (#2[[1]]-1)+(#2[[2]]-1)]<>"]",myCForm[#1]]&,Partition[code1[[-1]],m],{2}];

codestringAddTo=StringCases[StringReplace[GenerateCode[CBlock[code1],Indent->5],"Sqrt("->"sqrt("],Longest["{\n"~~x___~~"\n}"]:>x][[1]];


StringJoin["
    void ",name,"( 
		const Tensor2<Real,Int> & V_coords, 
		const Tensor2<Int ,Int> & simplices, 
		const Tensor2<Real,Int> & P_D_far, 
        // cppcheck-suppress [constParameter]
		      Tensor3<Real,Int> & buffer, 
		bool addTo 
	) const
    {
        ptic(ClassName()+\"::"<>name<>"\");

        if( P_D_far.Dimension(1) != "<>FarDim<>" )
        {
            eprint(\"in "<>name<>": P_D_far.Dimension(1) != "<>FarDim<>". Aborting\");
        }
        
		ptr<Real> V_coords__  = V_coords.data();
        ptr<Int>  simplices__ = simplices.data();
		ptr<Real> P_D_far__   = P_D_far.data();
        mut<Real> buffer__    = buffer.data();

		if( addTo )
		{
			ParallelDo(
				[=]( const Int i )
				{
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
