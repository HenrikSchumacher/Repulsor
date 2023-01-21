(* ::Package:: *)

ClearAll[line];

$MaxDomainDimension=2;
$MaxAmbientDimension=3;

line[k_,s___]:=StringJoin[ConstantArray["\t",k],s,";\n"];

CDeclareAssign::usage = "CDeclareAssign[type, var, value] represents 'type var = value;'.";

SymbolicC`Private`IsCExpression[ _CDeclareAssign ] := True;

CDeclareAssign::usage="";
GenerateCode[CDeclareAssign[typeArg_, idArg_, rhs_], opts : OptionsPattern[]] :=
    Module[{type, id},
      type = Flatten[{typeArg}];
      id = Flatten[{idArg}];
      type = Riffle[ Map[ GenerateCode[#, opts] &, type], " "];
      id = Riffle[ Map[ GenerateCode[#, opts] &, id], ", "];
      GenerateCode[CAssign[type <> " " <> id, rhs], opts]
    ];

CAddTo::usage = "CAddTo[var, value] represents 'var += value;'.";

SymbolicC`Private`IsCExpression[ _CAddTo ] := True;

GenerateCode[CAddTo[lhs_, rhs_], opts : OptionsPattern[]] :=  GenerateCode[lhs,opts]<>" += "<>GenerateCode[rhs,opts];


toC[a_String]:=a;
toC[a_]:=ToString[a,CForm];

arule={
	CAddTo[a_,b_]:>CAddTo[toC[a],toC[b]],
	CAssign[a_,b_]:>CAssign[toC[a],toC[b]],
	CDeclareAssign["const Real",a_,b_]:>CDeclareAssign["const Real",toC[a],toC[b]]
};

myCForm=StringReplace[ToString[#,CForm],"\""->""]&;

floatcast[x_]:="static_cast<Real>("<>x<>")";
