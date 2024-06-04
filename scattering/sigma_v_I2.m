(* ::Package:: *)

(* ::Input::Initialization:: *)
lightspeed = 299792;
MeVm3tocm2g=4.57821356*10^(-6);

(**PATHTOCODE ="/home/XXX"**)
PATHTOCODE =Directory[]
(**numsteps = 200;**)
numsteps = 200;
vmin =18090;
vmax = lightspeed*0.99;
massSp4units=100;(*MeV*)
MJdistgamma[gamma_, theta_]:= (E^(-(gamma/theta)) gamma^2*Sqrt[1-1/gamma^2])/(theta BesselK[2,1/theta])
vmeanthetaexact[theta_]:=2 (theta (1+theta)(E^(-1/theta)) )/BesselK[2,1/theta]
vmeanthetaexp[theta_]:=Sqrt[theta]*Sqrt[8/Pi]-Sqrt[theta]^3*7/Sqrt[8*Pi]+Sqrt[theta]^5*105/(32*Sqrt[2*Pi])-Sqrt[theta]^7*525/(256*Sqrt[2*Pi])-Sqrt[theta]^9*9765/(8192*Sqrt[2*Pi])
vmeantheta[theta_]:=If[theta<0.00241148,vmeanthetaexp[theta],vmeanthetaexact[theta]]
thetavmean[vmean_]:=If[vmean<6*10^(-6),0,theta/.FindRoot[vmeantheta[theta]-vmean, {theta,0.1}]]
MJdistv[v_,vmean_]:=MJdistgamma[gammaofv[v],thetavmean[vmean]]

gammaofv[v_]:=1/Sqrt[1-v^2]
vofgamma[gamma_]:=Sqrt[1-1/gamma^2]
Pofv[v_,mass_]:=Piecewise[{{v*mass*gammaofv[v],0<v<1}},Indeterminate]
vofP[P_,mass_]:=P/Sqrt[mass^2+P^2]
sigmav[a_,re_,v_,mass_]:=\!\(\*
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
RowBox[{"a", " ", 
SuperscriptBox["mass", "2"], " ", "re", " ", 
SuperscriptBox["v", "2"]}], 
RowBox[{"2", " ", 
RowBox[{"(", 
RowBox[{"1", "-", 
SuperscriptBox["v", "2"]}], ")"}]}]], "+", 
FractionBox[
RowBox[{"I", " ", "a", " ", "mass", " ", "v"}], 
SqrtBox[
RowBox[{"1", "-", 
SuperscriptBox["v", "2"]}]]]}], "]"}], "2"]], "]"}], 
RowBox[{"0", "<", "v", "<", "1"}]},
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
intfuncMJ[a_,re_,mass_,gamma_,theta_]:=vofgamma[gamma]*MJdistgamma[gamma,theta]*sigmav[a,re,vofgamma[gamma],mass]
IntedfuncMJtheta[a_,re_,mass_,theta_]:=NIntegrate[intfuncMJ[a,re,mass,gamma,theta],{gamma,1,5},MinRecursion->12,MaxRecursion->16,AccuracyGoal->5]
IntedfuncMJvmean[a_,re_,mass_,vmean_]:=IntedfuncMJtheta[a,re,mass,thetavmean[vmean]]
Sp4data =  Import[PATHTOCODE<>"/output/Sp(4)_data.csv","Table"];
aSp4units=Sp4data[[All,3]]/massSp4units;(*MeV^-1*)
reSp4units=Sp4data[[All,4]]/massSp4units;(*MeV^-1*)
aofmpfpichiPT[mpifpi2_]:=mpifpi2/32
(**achilow = aofmpfpichiPT[Min[Sp4data[[All,6]]]^2]/massSp4units
achihigh = aofmpfpichiPT[Max[Sp4data[[All,6]]]^2]/massSp4units**)
achilow = aofmpfpichiPT[4.6^2]/massSp4units
achihigh = aofmpfpichiPT[6.0^2]/massSp4units
vfunc[vmin_,vmax_,steps_,i_]:= vmin*(vmax/vmin)^((i-1)/(steps-1))//N;
Clear[Exporttable]
Exporttable = List[];
AppendTo[Exporttable,Table[vfunc[vmin,vmax,numsteps,i],{i,1,numsteps}]];
AppendTo[Exporttable,Table[IntedfuncMJvmean[achilow,0,massSp4units,vfunc[vmin,vmax,numsteps,i]/lightspeed]*lightspeed/MeVm3tocm2g/massSp4units,{i,1,numsteps}]];
AppendTo[Exporttable,Table[IntedfuncMJvmean[achihigh,0,massSp4units,vfunc[vmin,vmax,numsteps,i]/lightspeed]*lightspeed/MeVm3tocm2g/massSp4units,{i,1,numsteps}]];
For[i=1,i<Length[aSp4units]+1,i++,(a=Table[IntedfuncMJvmean[aSp4units[[i]],reSp4units[[i]],massSp4units,vfunc[vmin,vmax,numsteps,j]/lightspeed]*lightspeed/MeVm3tocm2g/massSp4units,{j,1,numsteps}];
AppendTo[Exporttable,a])];
Export[PATHTOCODE<>"/output/sigma_v_data.csv",Transpose[Exporttable]];
Quit[]

