(* ::Package:: *)

(* ::Input::Initialization:: *)
lightspeed = 299792;
MeVm3tocm2g=4.57821356*10^(-6);
massSp4units=100;(*MeV*)
massLSunits=16700;(*MeV*)
MeVfm=197.3;
PATHTOCODE= Directory[];
vfunc[vmin_,vmax_,steps_,i_]:= vmin*(vmax/vmin)^((i-1)/(steps-1))//N;
numsteps = 200;
vmin =25;
vmax =2800;
sigmav[a_,re_,\[CapitalDelta]v_,mass_]:=\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
RowBox[{"Re", "[", 
FractionBox[
RowBox[{"4", " ", 
SuperscriptBox["a", "2"], " ", "\[Pi]"}], 
SuperscriptBox[
RowBox[{"Abs", "[", 
RowBox[{"1", "-", 
FractionBox[
RowBox[{"a", " ", "re", " ", 
SuperscriptBox["mass", "2"], " ", 
SuperscriptBox["\[CapitalDelta]v", "2"]}], 
RowBox[{"8", 
RowBox[{"(", 
RowBox[{"1", "-", 
RowBox[{
RowBox[{"\[CapitalDelta]v", "^", "2"}], "/", "4"}]}], ")"}]}]], "+", 
FractionBox[
RowBox[{"I", " ", "a", " ", "mass", " ", "\[CapitalDelta]v"}], 
RowBox[{"2", 
SqrtBox[
RowBox[{"1", "-", 
RowBox[{
RowBox[{"\[CapitalDelta]v", "^", "2"}], "/", "4"}]}]]}]]}], "]"}], "2"]], "]"}], 
RowBox[{"0", "<", "\[CapitalDelta]v", "<", "1"}]},
{"0", 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)
MBdistvmean[\[CapitalDelta]v_,vmean_]:=(32(\[CapitalDelta]v^2) (E^(-((4 \[CapitalDelta]v^2)/(\[Pi] vmean^2)))) )/(\[Pi]^2 vmean^3)
intfuncMB[a_,re_,mass_,\[CapitalDelta]v_,vmean_]:=\[CapitalDelta]v*MBdistvmean[\[CapitalDelta]v,vmean]*sigmav[a,re,\[CapitalDelta]v,mass]
IntedfuncMBvmean[a_,re_,mass_,vmean_]:=NIntegrate[intfuncMB[a,re,mass,\[CapitalDelta]v,vmean],{\[CapitalDelta]v,0,1},MinRecursion->18,MaxRecursion->24,AccuracyGoal->9]
Sp4data =  Import[PATHTOCODE<>"/output/tables/Sp(4)_data.csv","CSV"];
aSp4units=Sp4data[[All,3]]/massSp4units;(*MeV^-1*)
reSp4units=Sp4data[[All,4]]/massSp4units;(*MeV^-1*)
aofmpfpichiPT[mpifpi2_]:=mpifpi2/32
achilow = aofmpfpichiPT[4.6^2]/massSp4units;
achihigh = aofmpfpichiPT[6.0^2]/massSp4units;
aLS=22.2/MeVfm;(*MeV^-1*)
reLS=-2.59*10^(-3)/MeVfm;(*MeV^-1*)(*MeV^-1*)
Clear[Exporttable]
Exporttable = List[];
AppendTo[Exporttable,Table[vfunc[vmin,vmax,numsteps,i],{i,1,numsteps}]];
AppendTo[Exporttable,Table[IntedfuncMBvmean[achilow,0,massSp4units,vfunc[vmin,vmax,numsteps,i]/lightspeed]*lightspeed/MeVm3tocm2g/massSp4units,{i,1,numsteps}]];
AppendTo[Exporttable,Table[IntedfuncMBvmean[achihigh,0,massSp4units,vfunc[vmin,vmax,numsteps,i]/lightspeed]*lightspeed/MeVm3tocm2g/massSp4units,{i,1,numsteps}]];
For[i=1,i<Length[aSp4units]+1,i++,(a=Table[IntedfuncMBvmean[aSp4units[[i]],reSp4units[[i]],massSp4units,vfunc[vmin,vmax,numsteps,j]/lightspeed]*lightspeed/MeVm3tocm2g/massSp4units,{j,1,numsteps}];
AppendTo[Exporttable,a])];
AppendTo[Exporttable,Table[IntedfuncMBvmean[aLS,reLS,massLSunits,vfunc[vmin,vmax,numsteps,i]/lightspeed]*lightspeed/MeVm3tocm2g/massLSunits,{i,1,numsteps}]];
Export[PATHTOCODE<>"/output/tables/sigma_v_data.csv",Transpose[Exporttable]];
