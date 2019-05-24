(* ::Package:: *)

(* ::Section:: *)
(*Velocity Gradient \[Del] *)


(*None Type Arithmetic Rules *)
Unprotect[None,Power];
None/:Subtract[None,None]:=None;
None/:Subtract[None,_Real|_Integer]:=None;
None/:Subtract[_Real|_Integer,None]:=None;
None/:Plus[None,-None]:=None;
None/:Plus[None,_Real|_Integer]:=None;
None/:Times[Rational[_,_],None]:=None;
None/:Times[_Integer|_Real,None]:=None;
None/:Times[None,None]:=None;
Power/:Times[_,Power[None,_]]:=None;
Protect[None,Power];


Unevaluated[x]/:grad[Unevaluated[x],array_,div_]:=Block[{arr,nx,ny,diffL,diffR,diff},
{nx,ny}={#[[1]],#[[2]]}&[Dimensions[array]];
arr=PadLeft[PadRight[array,{nx,ny+1},None],{nx,ny+2},None];
diffL=arr[[All,2;;ny+1]]-arr[[All,1;;ny]];
diffR=arr[[All,3;;ny+2]]-arr[[All,2;;ny+1]];
diff=(diffL+diffR)/2;
diff/div
];


Unevaluated[y]/:grad[Unevaluated[y],array_,div_]:=Block[{arr,nx,ny,diffL,diffR,diff},
{nx,ny}={#[[1]],#[[2]]}&[Dimensions[array]];
arr=PadLeft[PadRight[array,{nx+1,ny},None],{nx+2,ny},None];
diffL=arr[[2;;nx+1,All]]-arr[[1;;nx,All]];
diffR=arr[[3;;nx+2,All]]-arr[[2;;nx+1,All]];
diff=(diffL+diffR)/2;
diff/div
];


velocityGrad[dgrad_,vectorfield_]:=Module[{pts,vel,xmin,xmax,ymin,ymax,array,grid,
ux,vx,uy,vy,DVxx,DVyx,DVxy,DVyy,DVww,trace},
{pts,vel}=Transpose[vectorfield];
{ymin,ymax}=MinMax@pts[[All,1]];
{xmin,xmax}=MinMax@pts[[All,2]];
array = meshgrid[Range[ymin,ymax,dgrad],Range[xmin,xmax,dgrad]];
grid = Replace[Replace[array,Dispatch[Thread[pts->vel]],{2}],{__Integer}->{None,None},{2}];
ux = grad[Unevaluated[x],grid[[All,All,2]],dgrad];
vx = grad[Unevaluated[x],grid[[All,All,1]],dgrad];
uy = grad[Unevaluated[y],grid[[All,All,2]],dgrad];
vy = grad[Unevaluated[y],grid[[All,All,1]],dgrad];
DVxx = (ux);
DVyx = DVxy =(uy+vx)/2;
DVyy = (vy);
DVww = (-uy+vx); (*rotation*)
trace = (DVxx+DVyy)/2;
Print@Thread[{"xx","yy","yx","xy","ww","trace"}->
(MatrixPlot[#,ColorFunction->"Rainbow",ColorRules->{0->Black}]&/@({DVxx,DVyy,DVyx,DVxy,DVww,trace}/. None->0))];
{DVxx,DVyy,DVyx,DVxy,DVww,trace,array}
];


(* ::Section:: *)
(*Strain Rate Measures*)


SetAttributes[makeVecs,HoldAll];
makeVecs[vals_,vecs_,s_,rQ_]:=With[{n=Length[vals]},
If[TrueQ[rQ],
Do[Switch[Sign[Im[vals[[k]]]],
0, Null,
1, vecs[[k]]+=s I vecs[[k+1]],
-1, vecs[[k]]=Conjugate[vecs[[k-1]]]
],{k,n}]]];

Options[getEigensystem]={Mode->Automatic};
getEigensystem[mat_?SquareMatrixQ,opts:OptionsPattern[]]:=Module[{m=mat,chk,ei,ev,lm,lv,n,rm,rQ,rv},
n=Length[mat];
Switch[OptionValue[Mode],
Right|Automatic,
{lm,rm}={"N","V"},
Left,{lm,rm}={"V","N"},
All,{lm,rm}={"V","V"},
_,{lm,rm}={"N","V"}];
LinearAlgebra`LAPACK`GEEV[lm,rm,m,ev,ei,lv,rv,chk,rQ];
If[!TrueQ[chk],
Message[getEigensystem::eivec0];
Return[$Failed,Module]
];
If[rQ,ev+=I ei];
Switch[OptionValue[Mode],Right|Automatic,rv=ConjugateTranspose[rv];
makeVecs[ev,rv,1,rQ];{ev,rv},Left,lv=ConjugateTranspose[lv];
makeVecs[ev,lv,-1,rQ];{ev,lv},All,{lv,rv}=ConjugateTranspose/@{lv,rv};
makeVecs[ev,rv,1,rQ];makeVecs[ev,lv,-1,rQ];
{ev,rv,lv},_,rv=ConjugateTranspose[rv];
makeVecs[ev,rv,1,rQ];{ev,rv}]
];


plotbars[pt_,\[Phi]_,eigval_,scale_,inds_,imgdim_]:=Module[{Amp,x0,y0,xmaj,ymaj},
Amp=scale Abs[eigval];
{y0,x0}=pt;
{xmaj,ymaj}=Amp Through[{Cos,Sin}[\[Phi]]];
Line[Abs[{0,imgdim}-{##}]&@@@Thread[{xmaj {-1,1}+x0,ymaj {-1,1}+y0}]]
];


SRMeasures[img_,flow_,DVxx_,DVyy_,DVxy_,array_,scale_:500]:=Module[{graphicsPrimitive={},pos,DV,eVa,eVe,
vp1,val1,ind1,\[Phi],rvect,speed,anglevelmean,rotationTrans,rvectTurned,scalar,minDv,maxDv,Tracee,ptsTransfer,err,transFnPts,pts,
transVec,meanflowRot,imgdim=Last@ImageDimensions[img],pt,dir,\[Theta],vel},

pos = SparseArray[(DVxx*DVyy*DVxy)/. None -> 0]["NonzeroPositions"];
{pts,vel} = {flow[[All,1]],flow[[All,-1]]};
ptsTransfer = {#[[2]],imgdim-#[[1]]}&/@pts;
{err,transFnPts}=FindGeometricTransform[ptsTransfer,pts];
{pt,dir}=Mean/@{N@pts,vel};
\[Theta]=(ArcTan[#[[2]]/#[[1]]]&[dir]);
Do[
DV={{(DVxx[[Sequence@@i]]-DVyy[[Sequence@@i]])/2, DVxy[[Sequence@@i]]},
{DVxy[[Sequence@@i]], -(DVxx[[Sequence@@i]]-DVyy[[Sequence@@i]])/2}};
If[FreeQ[DV,None],
{eVa,eVe}={DiagonalMatrix[Reverse@#[[1]]],Reverse[-#[[2]]]}&[Re@getEigensystem[DV]];
(*{eVa,eVe}={DiagonalMatrix[#[[1]]],#[[2]]}&@Eigensystem[DV];*)
vp1 = Diagonal[eVa];
{val1,ind1}={Sort@vp1,Ordering@vp1};
\[Phi] = ArcTan[eVe[[2,ind1[[2]] ]]/eVe[[1,ind1[[2]] ]]];
rvect ={eVe[[1,ind1[[2]] ]],eVe[[2,ind1[[2]] ]]};

rotationTrans=RotationTransform[\[Theta]];
rvectTurned=rotationTrans[rvect];
scalar=Abs[First@rvectTurned];

{minDv,maxDv} = MinMax[eVa];
Tracee =(DVxx[[Sequence@@i]]+DVyy[[Sequence@@i]])/2;

Block[{prim,transpt=transFnPts@array[[Sequence@@i]]},
 If[Abs[minDv]*Abs[maxDv]>0,
 prim = GeometricTransformation[Circle[transpt,
 {Round[scale Abs[Tracee]],Round[scale Abs[Tracee]]}],
 RotationTransform[\[Phi],transpt]];

 If[Tracee>0,
 (*isotropic strain-rate*)
 AppendTo[graphicsPrimitive, Prepend[{prim},XYZColor[1,0,0,0.5]]],
 AppendTo[graphicsPrimitive, Prepend[{prim},XYZColor[0,0,1,0.5]]]
];
(*anisotropic strain-rate*)
 AppendTo[graphicsPrimitive,{ColorData["CandyColors"][scalar],Thickness[0.005],
 plotbars[array[[Sequence@@i]],\[Phi],maxDv,scale,i,imgdim]}]
  ]
 ]
],{i,pos}];
transVec=TransformationFunction[ReplacePart[Chop@transFnPts[[1]],{2,3}-> 0]];
meanflowRot=Thread[{ptsTransfer,transVec[#]&/@vel}];
{graphicsPrimitive,meanflowRot}
];


(* ::Section:: *)
(*Plotting F(x)s*)


(*plot the vector field along with the strain rates *)
plotStreamField[graphicsPrimitive_,flowfield_,image_,vecscale_:{0.075,0.30},imgres_:500,
imgsize_:500,arrowstretch_:200,arrowheadsize_:0.035]:=
Module[{pt,dir},
{pt,dir}=Mean/@{N@flowfield[[All,1]],flowfield[[All,2]]};
Image[Rasterize[Show[
ImageAdjust@image,
With[{reg=ConvexHullMesh@flowfield[[All,1]]},
ListStreamPlot[flowfield,StreamColorFunction->"Rainbow",StreamPoints->Fine,
RegionFunction->Function[{x,y,xu,yv},RegionMember[reg,{x,y}]],StreamScale->vecscale]
], 
Graphics[{GrayLevel[0.],Thick,Arrowheads[{{arrowheadsize,1,{Graphics[Line[{{-1, Rational[1, 2]}, {0, 0}, {-1, Rational[-1, 2]}, {-1, Rational[1, 2]}}], ImageSize -> {27.60000000000103, Automatic}],1}}}],Arrow[{pt,(pt+arrowstretch*dir)}]}],
ImageSize->Large],"Image",ImageResolution->imgres],
ImageSize->imgsize]
];


(*plot the vector field along with the strain rates *)
plotSRField[graphicsPrimitive_,flowfield_,image_,
imgres_:500,imgsize_:500,arrowstretch_:200,arrowheadsize_:0.035]:=
Module[{pt,dir,graphicsPrimitiveC,cases,posprims,pickprims,pickprimsOrig,plt},
{pt,dir}=Mean/@{N@flowfield[[All,1]],flowfield[[All,2]]};
Image[
Rasterize[Show[
ImageAdjust@image,
Graphics[If[!FreeQ[#,0,\[Infinity]],{Thick,#},{Thin,#}]&/@graphicsPrimitive],
Graphics[{GrayLevel[0],Thick,Arrowheads[{{arrowheadsize,1,{Graphics[Line[{{-1, Rational[1, 2]}, {0, 0}, {-1, Rational[-1, 2]}, {-1, Rational[1, 2]}}], ImageSize -> {27.60000000000103, Automatic}],1}}}],Arrow[{pt,(pt+arrowstretch dir)}]}],
ImageSize->Large],"Image",ImageResolution->imgres],
ImageSize->imgsize]
];
