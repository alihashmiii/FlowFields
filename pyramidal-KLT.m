(* ::Package:: *)

(* ::Section::Closed:: *)
(*Build Image Pyramid*)


conv2[A_?NumericQ,B_]:=conv2[{{A}},B];
conv2[A_?VectorQ,B_]:=conv2[{A},B];
conv2[A_,B_]:=ListConvolve[A,B,{1,-1},0];


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


(* ::Section::Closed:: *)
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


pyrLK::description="pyrLK tracks'features' from images in pyr1 to pyr2. pt0 contains the features to track. winsize is the size
of the window to estimate local gradient. iterate at most maxIter times per pyramid level or stop if convergence is less than 
the stop threshold (in pixels). winsize, maxiter and threshold can be scalars or vectors with a different value for each pyramid
level.

Output:
speed -> estimated speed of features (in lower pyramid image)
failure -> tracking failure
error -> final difference of pixel color.

'failure' is an array containing 2 counts of failures for each particle. The first is the number of times LK has failed due to a
weak gradient of color intensity in the area around the features. The second is the number of times it has been tracked out of the
image and forced back inside, or lost due to algorithmic failure (this should not happen anymore, otherwise warnings are displayed)";
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
 sp = pts0[[All,3;;4]] /(2^pyrNum);
pts0=pts0[[All,1;;2]],
sp=ConstantArray[0,{First[Dimensions@pts0],2}]
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
sp=sp + pts;
ind2=Position[sp[[All,1]],x_/;x>(First[Dimensions@img1]-dsval),{2}];
ind2=Thread[{Flatten@ind2,1}];
sp=ReplacePart[sp,ind2->First[Dimensions@img1]-dsval];
ind2=ind2[[All,1]];
ind=Union[ind,ind2];
ind2=Position[sp[[All,2]],x_/;x>(Last[Dimensions@img1]-dsval),{2}];
ind2=Thread[{ind2,2}];
sp=ReplacePart[sp,ind2->Last[Dimensions@img1]-dsval];
sp=sp - pts;
ind2=ind2[[All,1]];
ind=Union[ind,ind2];
Scan[(warn[[Sequence@@#]]+=1)&,Thread[{ind,2}]],{i,pyrNum,1,-1}];
On[InterpolatingFunction::dmval];
{sp,warn}
];


(* ::Section::Closed:: *)
(*Main *)


Options[KLTracker]= {"dgrid" -> 8, "threshold" -> 0.99, "windowFilter" -> 30, "pyramidNum" -> 2, "winSize" -> {20,20},
"maxIterations" -> {10,10}, "kthreshold" -> {0.05,0.05}, "blurRadius" -> 1};


KLTracker[images_,imagedata_,images_,masks_,istart_,iend_,dt_,OptionsPattern[]]:= Block[{index,time,nX,nY,img,
y,x,xeul,yeul,indexpoints,mask,boxmeanvals,pyr1,pyr1X,pyr1Y,pyr2,pyr2X,pyr2Y,ptsToTrack,sol,warn,lenImages = Length@images,
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
With[{p=Range[winSize[[i]],nX-winSize[[i]],dgrid]},
x=Array[p&,Length[p]]
];
With[{p=Range[winSize[[i]],nY-winSize[[i]],dgrid]},
y=Transpose@Array[p&,Length[p]]
];
With[{dim=Times@@Dimensions[x]},
xeul=ArrayReshape[x,dim];
yeul=ArrayReshape[y,dim];
indexpoints=ConstantArray[None,{dim,lenImages}]
];
];

img=images[[i]]; mask=masks[[i]];
With[{filtsize=Round[windowFilter/2]},
boxmeanvals=ParallelTable[
Mean@ImageTake[mask,{xeul[[j]]-filtsize,xeul[[j]]+filtsize},
{yeul[[j]]-filtsize,yeul[[j]]+filtsize}],
{j,Length@xeul}]
];
indexpoints[[All,i]] = Boole[Thread[boxmeanvals>threshold]];
(* pts to track *)
ptsToTrack = Thread[{Extract[xeul,#],Extract[yeul,#]}]&[Position[indexpoints[[All,i]],1]];
(* create image pyramid for the first image *)
{pyr1,pyr1X,pyr1Y} = makePyramid[imagedata[[i]],pyramidNum,blurRadius];
{pyr2,pyr2X,pyr2Y} = makePyramid[imagedata[[i+dt]],pyramidNum,blurRadius];
{sol,warn} = pyrLK[{pyr1,pyr1X,pyr1Y},{pyr2,pyr2X,pyr2Y},masks[[i]],ptsToTrack,winSize,maxIterations,kthreshold];
warn=Total/@warn;
sol=sol/dt;
Sow[{ptsToTrack,sol}],{i,istart,iend-dt}]
];


(* ::Section::Closed:: *)
(*Misc F(x)*)


meanFlowField[flowfield_]:=Normal@GroupBy[Flatten[Thread/@flowfield,1],First->Last,Mean]/.Rule-> List;


movingAvgFlowField[flowfield_,block_:5,trans_:1]:=BlockMap[
 Normal@GroupBy[Flatten[Thread/@#,1],First->Last,Mean]/.Rule ->List &,
flowfield,block,trans];


(* ::Section::Closed:: *)
(*Velocity Gradient \[Del] *)


Unevaluated[x]/:grad[Unevaluated[x],array_,div_]:=Block[{arr,nx,ny,diffL,diffR,diff},
{nx,ny}={#[[1]],#[[2]]}&[Dimensions[array]];
arr=PadLeft[PadRight[array,{nx,ny+1}],{nx,ny+2}];
diffL=arr[[All,2;;ny+1]]-arr[[All,1;;ny]];
diffR=arr[[All,3;;ny+2]]-arr[[All,2;;ny+1]];
diff=(diffL+diffR)/2;
diff/div
];


Unevaluated[y]/:grad[Unevaluated[y],array_,div_]:=Block[{arr,nx,ny,diffL,diffR,diff},
{nx,ny}={#[[1]],#[[2]]}&[Dimensions[array]];
arr=PadLeft[PadRight[array,{nx+1,ny}],{nx+2,ny}];
diffL=arr[[2;;nx+1,All]]-arr[[1;;nx,All]];
diffR=arr[[3;;nx+2,All]]-arr[[2;;nx+1,All]];
diff=(diffL+diffR)/2;
diff/div
];


velocityGrad[dgrad_,vectorfield_]:=Module[{pts,vel,xmin,xmax,ymin,ymax,array,grid,
ux,vx,uy,vy,DVxx,DVyx,DVxy,DVyy,DVww,trace},
{pts,vel}=Transpose[vectorfield];
{xmin,xmax}=MinMax@pts[[All,1]];
{ymin,ymax}=MinMax@pts[[All,2]];
array = meshgrid[Range[xmin,xmax,dgrad],Range[ymin,ymax,dgrad]];
grid = Replace[Replace[array,Dispatch[Thread[pts -> vel]],{2}],{__Integer}->{0,0},{2}];
ux = grad[Unevaluated[x],grid[[All,All,1]],dgrad];
vx = grad[Unevaluated[x],grid[[All,All,2]],dgrad];
uy = grad[Unevaluated[y],grid[[All,All,1]],dgrad];
vy = grad[Unevaluated[y],grid[[All,All,2]],dgrad];
DVxx = ux/. 0-> None;
DVyx = DVxy =((uy+vx)/2)/. 0 -> None;
DVyy = vy/. 0 -> None;
DVww = (-uy+vx)/. 0 -> None; (*rotation*)
trace = (DVxx+DVyy)/2;
Print[MatrixPlot[#,ColorFunction->"Rainbow",ColorRules->{0->Black}]&/@({DVxx,DVyx,DVyy,DVww,trace}/. None -> 0)];
{DVxx,DVyy,DVxy,DVyx,DVww,trace,array}
];


(* ::Section::Closed:: *)
(*Strain Rate Measures*)


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


(* ::Code::Initialization:: *)
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


SRMeasures[DVxx_,DVyy_,DVxy_,vectorfield_,array_]:=Module[{graphicsPrimitive={},pos,DV,eVa,eVe,vp1,val1,ind1,
\[Phi],rvect,speed,anglevelmean,rotationTrans,rvectTurned,scalar,minDv,maxDv,Tracee},
pos=SparseArray[(DVxx*DVyy*DVxy)/.None->0]["NonzeroPositions"];

Do[
DV={{(DVxx[[Sequence@@i]]-DVyy[[Sequence@@i]])/2,DVxy[[Sequence@@i]]},
{DVxy[[Sequence@@i]],-(DVxx[[Sequence@@i]]-DVyy[[Sequence@@i]])/2}};

If[FreeQ[DV,None],
{eVa,eVe}={DiagonalMatrix@#[[1]],#[[2]]\[Transpose]}&[Re@getEigensystem[DV]];
vp1 = Diagonal[eVa];
{val1,ind1}={Sort@vp1,Ordering@vp1};
\[Phi] = ArcTan[eVe[[2,ind1[[2]] ]]/eVe[[1,ind1[[2]] ]]];
rvect ={eVe[[1,ind1[[2]] ]],eVe[[2,ind1[[2]] ]]};
speed=vectorfield[[All,2]];
anglevelmean = 45 Degree; (* angle from polarization code maybe used instead*)
rotationTrans = RotationTransform[anglevelmean];
rvectTurned = rotationTrans[rvect];
scalar = Abs[First@rvectTurned];
minDv = eVa[[1,1]]~Min~eVa[[2,2]];
maxDv = eVa[[1,1]]~Max~eVa[[2,2]];
Tracee =(DVxx[[Sequence@@i]]+DVyy[[Sequence@@i]])/2;

Block[{prim},
If[Abs[minDv]*Abs[maxDv]>0,
prim = GeometricTransformation[Circle[array[[Sequence@@i]],{Round[1000 Abs[Tracee]],0}],
RotationTransform[\[Phi],array[[Sequence@@i]]] ];

If[Tracee>0,
AppendTo[graphicsPrimitive, Prepend[{prim},XYZColor[0,0,1,0.33]]],
AppendTo[graphicsPrimitive, Prepend[{prim},XYZColor[1,0,0,0.33]]]
];

AppendTo[graphicsPrimitive, {GrayLevel[0.1],
 GeometricTransformation[ Circle[array[[Sequence@@i]],{Round[360*maxDv],Round[Abs[360*minDv]]}],
RotationTransform[\[Phi],array[[Sequence@@i]]]]}]
  ]
 ]
],{i,pos}];
graphicsPrimitive
];


(* ::Section::Closed:: *)
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


(*plot the vector field along with the strain rates *)
plotVectorFieldSR[graphicsPrimitive_,flowfield_,image_,trim_:True,linethickness_:0.006,
vecscale_:{0.075,0.30},imgres_:500,imgsize_:500,arrowstretch_:200,arrowheadsize_:0.025]:=
Module[{pt,dir,graphicsPrimitiveC,cases,posprims,pickprims,pickprimsOrig,plt},
{pt,dir}=Mean/@{N@flowfield[[All,1]],flowfield[[All,2]]};

plt = Image[ImageRotate[
Rasterize[Show[
ImageRotate[ImageAdjust@image, Pi/2],
Graphics[If[FreeQ[#,0,\[Infinity]],{Thick,#},{Thickness[linethickness],#}]&/@graphicsPrimitive],
With[{reg=ConvexHullMesh@flowfield[[All,1]]},
ListVectorPlot[flowfield,VectorColorFunction->"Rainbow",VectorPoints->Fine,
RegionFunction->Function[{x,y,xu,yv},RegionMember[reg,{x,y}]],VectorScale->vecscale]
],
Graphics[{GrayLevel[0.25],Thick,Arrowheads[{{arrowheadsize,1,{Graphics[Line[{{-1, Rational[1, 2]}, {0, 0},
{-1, Rational[-1, 2]}, {-1, Rational[1, 2]}}], ImageSize -> {27.60000000000103, Automatic}],1}}}],Arrow[{pt,(pt+arrowstretch*dir)}]}],
ImageSize->Large],"Image",ImageResolution->imgres],-Pi/2],ImageSize->imgsize];

If[trim,
graphicsPrimitiveC=Replace[graphicsPrimitive, 
 HoldPattern[{t_,u_[Circle[v_,{w_Integer,0}],x_]}]:> {t,u[Circle[v,{w,1}],x]},{1}];
cases=Cases[graphicsPrimitiveC,Circle[_,{_,_}],{3}];
posprims=Position[Map[Length][Quiet@RandomPoint[#,200]&/@cases],x_/;x==200];
pickprims=Extract[graphicsPrimitiveC,posprims];
pickprimsOrig=Extract[graphicsPrimitive,posprims];
pickprimsOrig=Quiet@With[{reg=ConvexHullMesh@flowfield[[All,1]]},
Pick[pickprimsOrig,(Or@@RegionMember[reg,RandomPoint[#,200]]&/@Cases[pickprims,Circle[_,{_,_}],{3}]),True]
];

plt=Image[ImageRotate[
Rasterize[Show[
ImageRotate[ImageAdjust@image,Pi/2],
Graphics[If[FreeQ[#,0,\[Infinity]],{Thick,#},{Thickness[linethickness],#}]&/@pickprimsOrig],
With[{reg=ConvexHullMesh@flowfield[[All,1]]},
ListVectorPlot[flowfield,VectorColorFunction->"Rainbow",VectorPoints->Fine,
RegionFunction->Function[{x,y,xu,yv},RegionMember[reg,{x,y}]],VectorScale->vecscale]
],
Graphics[{GrayLevel[0.25],Thick,Arrowheads[{{arrowheadsize,1,{Graphics[Line[{{-1, Rational[1, 2]}, {0, 0}, {-1, Rational[-1, 2]},
{-1, Rational[1, 2]}}], ImageSize -> {27.60000000000103, Automatic}],1}}}],Arrow[{pt,(pt+arrowstretch*dir)}]}],
ImageSize->Large],"Image",ImageResolution->imgres],-Pi/2],
ImageSize->imgsize]
];
plt
];
