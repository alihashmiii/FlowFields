(* ::Package:: *)

(* ::Section:: *)
(*Build Image Pyramid*)


conv2[A_?NumericQ,B_]:=conv2[{{A}},B];
conv2[A_?VectorQ,B_]:=conv2[{A},B];
conv2[A_,B_]:=ListConvolve[A,B,{1,-1},0];


(* 
'level' is the number of level of the output pyramid; 'blurRadius' \[Rule] convolve pyramid images by a blur radius 
As output we get pyramids at level i (as in [Bouguet 2002]) and gradient Images.
*)
makePyramid[image_,level_Integer,blurRadius_]:=Module[{filter,F,imdata,ind0x,ind0y,ind1x,ind1y,
ind2x,ind2y,simg,pyramid,pyramidgradX,pyramidgradY,gradX,gradY,sx,sy},
pyramid=pyramidgradX=pyramidgradY={};
If[blurRadius>0,
filter=1;
F=conv2[{1,2,1},List/@{1,2,1}]/16.;
filter=Nest[conv2[#,F]&,filter,blurRadius]
];
Do[
If[i==1,
imdata=image,
imdata=pyramid[[i-1]];
{sx,sy}=Dimensions[imdata];
{ind0x,ind0y } = Range[1,#,2]&/@{sx-1,sy-1};
{ind1x,ind1y} = Range[2,#,2]&/@{sx,sy};
{ind2x,ind2y} = Range[3,#,2]&/@{sx,sy};

If[Length[ind2x]<Length[ind1x],AppendTo[ind2x,Last@ind2x]];
If[Length[ind2y]<Length[ind1y],AppendTo[ind2y,Last@ind2y]];

simg=0.25*imdata[[ind1x,ind1y]];
simg=Fold[#1+0.125imdata[[#2[[1]],#2[[2]]]]&,simg,{{ind0x,ind1y},{ind1x,ind0y},{ind2x,ind1y},
{ind1x,ind2y}}];
imdata=Fold[#1+0.0625imdata[[#2[[1]],#2[[2]]]]&,simg,{{ind0x,ind0y},{ind2x,ind0y},{ind2x,ind2y},
{ind0x,ind2y}}];
];

If[filter!= {},imdata=conv2[filter,imdata][[2;;-2,2;;-2]] ];
gradX = conv2[{{1,0,-1}},imdata][[All,2;;-2]];
gradY=Thread[conv2[{{1,0,-1}},Transpose@imdata]][[2;;-2]];
AppendTo[pyramid,imdata];
AppendTo[pyramidgradX,gradX];
AppendTo[pyramidgradY,gradY],{i,level}];
{pyramid,pyramidgradX,pyramidgradY}
];


(* ::Section:: *)
(*Pyramidal LK*)


meshgrid[int1_Integer,int2_Integer]:=Transpose[
{ConstantArray[Range[int1],int2],Transpose@ConstantArray[Range[int2],int1]},
{3,2,1}];
meshgrid[list1_List,list2_List]:=Transpose[
{ConstantArray[list1,Length@list2],Transpose@ConstantArray[list2,Length@list1]},
{3,2,1}];


interp2D[meshgrid_,array_]:=Block[{interp},
Off[InterpolatingFunction::dmval];
interp=Interpolation[Flatten[Transpose[{meshgrid,array},{3,2,1}],1],InterpolationOrder->1];
On[InterpolatingFunction::dmval];
interp
];


pyrLK::description="pyrLK tracks 'features' from images in pyr1 to pyr2. features contains the pts to track. winsize is the
size of the window to estimate local gradient. We iterate at most maxIter times for each pyramid or stop if convergence
criteria is met i.e. < some limit (in pixels)";
pyrLK[{pyr1_,pyr1x_,pyr1y_},{pyr2_,pyr2x_,pyr2y_},mask_,features_,winsize_:5,maxIter_:2,threshold_:0.5]:=Block[{pyrNum,
windowSize=winsize,thresh=threshold,maxiter=maxIter,pts0=features,pts,sp,ds,warn,img1,img2,Px,Py,count,r,ldspl,Sdsp,feature,
featureNa,featureNb,A,a,B,b,dsp,II,revdsp,ind,ind2,dsval,dim,interpFuncX,interpFuncY,mesh,ptsToForceBack,meshg,nearestFunc,
interpImg1,interpImg2,f},
Off[InterpolatingFunction::dmval];
(*mesh=ImageMesh[mask,Method\[Rule]"Exact"];*)
pyrNum=Length[pyr1];
If[Length[windowSize]<pyrNum,
windowSize=ConstantArray[windowSize,pyrNum-Length[windowSize]]];
If[Length[thresh]<pyrNum,
thresh=ConstantArray[thresh,pyrNum-Length[threshold]]];
If[Length[maxiter]<pyrNum,
maxiter=ConstantArray[maxiter,pyrNum-Length[threshold]]];
(*initialze some var: sp \[Rule] (displacement/speed) of features & ds \[Rule] radius of neighbourhood *)
If[Last@Dimensions[pts0]>= 4,
 sp = pts0[[All,3;;4]]/(2^pyrNum);
 pts0 = pts0[[All,1;;2]],
 sp = ConstantArray[0,{First[Dimensions@pts0],2}]
];
ds = Floor[windowSize/2];
windowSize=2 ds + 1;
warn = ConstantArray[0,{First[Dimensions@pts0],2}];
Do[
sp=2*sp;
pts=pts0/(2^(i-1));
{img1,img2}={pyr1[[i]],pyr2[[i]]};
{Px,Py}={pyr1x[[i]],pyr1y[[i]]};
(*neighbourhood of a feature*)
r=Tuples[Range[-ds[[i]],ds[[i]]],2];
dim=Dimensions@Px;
meshg=meshgrid[Sequence@@dim];
{interpFuncX,interpFuncY}={interp2D[meshg,Px],interp2D[meshg,Py]};
{interpImg1,interpImg2}={interp2D[meshg,img1],interp2D[meshg,img2]};
Do[
(*initialize variable use for current feature*)
count=0; (* iteration count *)
ldspl=thresh[[i]] + 1; (*norm of last iteration displacement *)
Sdsp={0,0}; (* \[CapitalSigma] of displacement evaluated at the level *)
feature=pts[[k]];
featureNa=(feature+#)&/@r;
a=Thread[{interpFuncX@@@featureNa,interpFuncY@@@featureNa}];
Check[
A = Inverse[a\[Transpose].a];
II=interpImg1@@@featureNa,
warn[[k,1]]+=1;
Continue[]
];
(* iterate until it converges (or fails to) *)
With[{maxitercount=maxiter[[i]],threshParam=thresh[[i]]},
f=(feature+#)&/@r;
While[(count<maxitercount)&&(ldspl>threshParam)&&Norm[Sdsp]<5,
featureNb=With[{speed=sp[[k]]},speed+#&/@f];
b=II -(interpImg2@@@featureNb);
B=(a\[Transpose]).b;
dsp=A.B;
ldspl=Norm[dsp];
revdsp=Reverse[dsp];
sp[[k]]+=revdsp;
Sdsp+=revdsp;
count+=1;
];
],{k,Length@pts}];
(*clamp estimated speed to stay in image *)
(*lower clamp*)
dsval=ds[[i]];
ind=Position[sp,x_/;x<(-pts+dsval),{2}];
Scan[(sp[[Sequence@@#]]=-pts[[Sequence@@#]]+dsval)&,ind];
ind=ind[[All,1]];
(*higher clamp*)
sp = sp+pts;
ind2=Position[sp[[All,1]],x_/;x>(First[Dimensions@img1]-dsval),{2}];
ind2=Thread[{Flatten@ind2,1}];
sp=ReplacePart[sp,ind2->First[Dimensions@img1]-dsval];
ind2=ind2[[All,1]];
ind=Union[ind,ind2];
ind2=Position[sp[[All,2]],x_/;x>(Last[Dimensions@img1]-dsval),{2}];
ind2=Thread[{ind2,2}];
sp=ReplacePart[sp,ind2->Last[Dimensions@img1]-dsval];
sp = sp-pts;
ind2=ind2[[All,1]];
ind=Union[ind,ind2];
Scan[(warn[[Sequence@@#]]+=1)&,Thread[{ind,2}]],{i,pyrNum,1,-1}];
On[InterpolatingFunction::dmval];
{sp,warn}
];


(* ::Section:: *)
(*Main *)


Options[KLTracker]= {"dgrid" -> 8, "threshold" -> 0.99, "windowFilter" -> 30, "pyramidNum" -> 2, "winSize" -> {20,20},
"maxIterations" -> {10,10}, "kthreshold" -> {0.05,0.05}, "blurRadius" -> 1};


KLTracker[images_,imagedata_,masks_,istart_,iend_,dt_,OptionsPattern[]]:= Block[{index,time,nX,nY,img,
y,x,mask,gridpts,boxmeanvals,pyr1,pyr1X,pyr1Y,pyr2,pyr2X,pyr2Y,ptsToTrack,sol,warn,lenImages = Length@images,
winSize=OptionValue["winSize"],dgrid=OptionValue["dgrid"],windowFilter=OptionValue["windowFilter"],
threshold=OptionValue["threshold"],pyramidNum=OptionValue["pyramidNum"],blurRadius=OptionValue["blurRadius"],
kthreshold=OptionValue["kthreshold"],maxIterations=OptionValue["maxIterations"]},

time=Array[None &, lenImages];
{nX,nY}=ImageDimensions[First@images];
First@Last@Reap@Do[
index=i-istart+1;
Print[index];
time[[index]]=i;

If[i==istart,
img=images[[i]];
x=Range[winSize[[i]],nX-winSize[[i]],dgrid];
y=Range[winSize[[i]],nY-winSize[[i]],dgrid];
gridpts=Flatten[meshgrid[x,y],1];
];

img=images[[i]]; mask=masks[[i]];

With[{filtsize=Round[windowFilter/2]},
boxmeanvals=ParallelTable[
Mean@ImageTake[mask,{pt[[2]]-filtsize,pt[[2]]+filtsize},
{pt[[1]]-filtsize,pt[[1]]+filtsize}],
{pt,gridpts}]
];
(* pts to track *)
ptsToTrack=Reverse[Pick[gridpts,Boole[Thread[boxmeanvals>threshold]],1],2];
(* create image pyramid for the first image *)
{pyr1,pyr1X,pyr1Y}=makePyramid[imagedata[[i]],2,1];
{pyr2,pyr2X,pyr2Y}=makePyramid[imagedata[[i+dt]],2,1];
{sol,warn}=pyrLK[{pyr1,pyr1X,pyr1Y},{pyr2,pyr2X,pyr2Y},masks[[i]],ptsToTrack,winSize,maxIterations,kthreshold];
warn=Total/@warn;
sol=sol/dt;
Sow[{ptsToTrack,sol}],{i,istart,iend-dt}]
];


(* ::Section:: *)
(*Misc F(x)*)


meanFlowField[flowfield_]:=Normal@GroupBy[Flatten[Thread/@flowfield,1],First->Last,Mean]/.Rule-> List;


movingAvgFlowField[flowfield_,block_:5,trans_:1]:=BlockMap[
 Normal@GroupBy[Flatten[Thread/@#,1],First->Last,Mean]/.Rule ->List &,
flowfield,block,trans];


(* ::Section:: *)
(*Plot F(x)*)


(*plot the vector field*)
plotVectorField[vectorfield_,image_,imageres_:250,imgsize_:Medium,vecscale_:{0.10,0.25}]:= 
With[{reg=ConvexHullMesh@vectorfield[[All,1]]},
Image[
ImageRotate[
Rasterize[Show[
ImageRotate[ImageAdjust@image,Pi/2],
ListVectorPlot[vectorfield, VectorColorFunction->"Rainbow", VectorPoints->Fine,
RegionFunction->Function[{x,y,xu,yv},RegionMember[reg,{x,y}]],VectorScale->vecscale],
ImageSize->Large],"Image",ImageResolution->imageres],
-Pi/2],
ImageSize->imgsize]
];


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
plotStreamField[flowfield_,image_,vecscale_:{0.075,0.30},imgres_:500,
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
