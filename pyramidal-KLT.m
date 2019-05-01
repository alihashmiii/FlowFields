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


(* ::Section:: *)
(*Main *)


Options[KLTracker]= {"dgrid" -> 8, "threshold" -> 0.99, "windowFilter" -> 30, "pyramidNum" -> 2, "winSize" -> {20,20},
"maxIterations" -> {10,10}, "kthreshold" -> {0.05,0.05}, "blurRadius" -> 1};


KLTracker[images_,imagedata_,images_,masks_,istart_,iend_,dt_,OptionsPattern[]]:= Block[{index,time,nX,nY,img,
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

img=images[[i]];mask=masks[[i]];

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


(* ::Section::Closed:: *)
(*Misc F(x)*)


meanFlowField[flowfield_]:=Normal@GroupBy[Flatten[Thread/@flowfield,1],First->Last,Mean]/.Rule-> List;


movingAvgFlowField[flowfield_,block_:5,trans_:1]:=BlockMap[
 Normal@GroupBy[Flatten[Thread/@#,1],First->Last,Mean]/.Rule ->List &,
flowfield,block,trans];


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
