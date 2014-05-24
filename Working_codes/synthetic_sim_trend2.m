%POSSIBILITY TO INCLUDE KNICKPOINT AND IDW PROFILE USING DIFFERENT P VALUES

%JCC 15072013

%=======================================================================%
%======================== IMPORT DATA & ORGANIZE =======================%
%=======================================================================%

%Here we import the aquifer outline and available geophysical data from
%external text files. NOTE: Cross sections should be specified in a
%south-north direction, in other words starting from the left boundary
%(if we look downstream)

top=importdata('topboundary_synth.txt')';
bottom=importdata('bottomboundary_synth.txt')';
cs1=importdata('cs1_synth_mod.txt')';
% cs2=importdata('cs2_synth.txt')';
% cs3=importdata('cs3_synth.txt')';
cs4=importdata('cs4_synth.txt')';

%================================================================================================%
%======= BEGIN DEFINE GLOBAL MODEL VARIABLES & CHOOSE PAIR OF KNOWN CROSS SECTIONS TO USE =======% 
%================================================================================================%

%This defines which cross sections we want to use for the model, and the
%internal discretization of the grid. We choose the number of interior
%cross sections and number of nodes. 

csLeft=cs1;
csRight=cs4;

numberOfPointsInt=100;
numberOfSections=100;
halfsectionPoints=100;

dir_csLeft=normc([csLeft(1,end)-csLeft(1,1);abs(csLeft(2,1)-csLeft(2,end))]);
dir_csRight=normc([csRight(1,end)-csRight(1,1);abs(csRight(2,1)-csRight(2,end))]);

%VALLEY BOUNDARY MODEL EXTENTS X DIRECTION
[Cmin,Imin]=min(top(1,:));
[Cmax,Imax]=max(top(1,:));

minX_Top=Cmin;
maxX_Top=Cmax;

[Cmin,Imin]=min(bottom(1,:));
[Cmax,Imax]=max(bottom(1,:));

minX_Bottom=Cmin;
maxX_Bottom=Cmax;

%XYZ VECTORS FOR TOP AND BOTTOM BOUNDARIES CONSTRUCTED FROM GIS DATA
x_Top=top(1,:);
y_Top=top(2,:);
z_Top=top(3,:);

x_Bottom=bottom(1,:);
y_Bottom=bottom(2,:);
z_Bottom=bottom(3,:);

n=10000;

%TOP BOUNDARY
Xtop=[(1:1:size(x_Top,2));x_Top];
Ytop=[(1:1:size(y_Top,2));y_Top];
Ztop=[(1:1:size(z_Top,2));z_Top];

XtopInterp=interp1(Xtop(1,:),Xtop(2,:),linspace(1,Xtop(1,end),n));
YtopInterp=interp1(Ytop(1,:),Ytop(2,:),linspace(1,Ytop(1,end),n));
ZtopInterp=interp1(Ztop(1,:),Ztop(2,:),linspace(1,Ztop(1,end),n));

%BOTTOM BOUNDARY
Xbottom=[(1:1:size(x_Bottom,2));x_Bottom];
Ybottom=[(1:1:size(y_Bottom,2));y_Bottom];
Zbottom=[(1:1:size(z_Bottom,2));z_Bottom];

XbottomInterp=interp1(Xbottom(1,:),Xbottom(2,:),linspace(1,Xbottom(1,end),n));
YbottomInterp=interp1(Ybottom(1,:),Ybottom(2,:),linspace(1,Ybottom(1,end),n));
ZbottomInterp=interp1(Zbottom(1,:),Zbottom(2,:),linspace(1,Zbottom(1,end),n));

%REASSIGN INTERPOLATED VALUES
x_Top=XtopInterp;
y_Top=YtopInterp;
z_Top=ZtopInterp;

x_Bottom=XbottomInterp;
y_Bottom=YbottomInterp;
z_Bottom=ZbottomInterp;

%PREPARE VECTOR DEFINING POLYGONAL REGION
xPoly=[x_Top fliplr(x_Bottom) x_Top(1,1)];
yPoly=[y_Top fliplr(y_Bottom) y_Top(1,1)];


%============================================================================%
%=============== CHANNEL CENTERLINE BINARY IMAGE PROCEDURES =================%
%============================================================================%

%This binarizes the domain, so we can apply the skeleotinization algorithm
%to obtain an appropriate centerline. 

%FIND POINTS IN POLYGON
crust=0.3;
dmesh=0.1;
[maxY,I]=max(top(2,:));
[minY,I]=min(bottom(2,:));

[maxX,I]=max(bottom(1,:));
[minX,I]=min(bottom(1,:));
[xx,yy]=meshgrid(minX-crust:dmesh:maxX+crust,minY-crust:dmesh:maxY+crust);

[IN] = inpolygon(xx(:),yy(:),xPoly,yPoly);
IN=reshape(IN,size(xx,1),size(yy,2));


%============================================================================%
%========================== APPLY SKELETON ALGORITHM ========================%
%============================================================================%

%After binarizing the aquifer domain, we extract the centerline. This
%centerline will allow us to guide the construction of an adaptive grid to
%the local geology.

%PREPROCESS
S1=double(IN);

%DO SKELETON
S1=bwmorph(S1,'skel',inf);
S2=bwmorph(S1,'spur',15);

%CROP CRUST
%[xx,yy]=meshgrid(minX:dmesh:maxX,minY:dmesh:maxY);
crustpixels=round(crust/dmesh);
S2=S2(crustpixels:end-crustpixels,crustpixels:end-crustpixels);
xx2=xx(crustpixels:end-crustpixels,crustpixels:end-crustpixels);
yy2=yy(crustpixels:end-crustpixels,crustpixels:end-crustpixels);

%CREATE AN XY VECTOR CONTAINING THE CHANNEL CENTERLINE
polyx=xx2(S2==1);
polyy=yy2(S2==1);

polyxs=polyx(1:10:end);
polyys=polyy(1:10:end);

A=linspace(2,numel(polyxs),1000);
polyxs=interp1(1:numel(polyxs),polyxs,A,'spline');
polyys=interp1(1:numel(polyys),polyys,A,'spline');


%================================================================================================%
%======================== PROCESS KNOWN DATA, ATTACH TO AQUIFER BOUNDARIES ======================%
%================================================================================================%

%This procedure attaches cross sections to boundaries by including closest
%points to the pointsXYZ_ matrix in order to achieve this, we need to 
%extend the known cross sections as close to the boundaries as possible 
%(using inpolygon) and then select the closest value.

% MATRICES FOR KNOWN CROSS SECTIONS 1 & 2
pointsXYZ_csLeft=zeros(3,size(csLeft,2)+2);
pointsXYZ_csRight=zeros(3,size(csRight,2)+2);

pointsXYZ_csLeft(1:3,2:size(csLeft,2)+1)=csLeft(1:3,:);
pointsXYZ_csRight(1:3,2:size(csRight,2)+1)=csRight(1:3,:);

%ATTACH BOTTOM BOUNDARY, CROSS SECTION RIGHT   
step=0.1;
in=1;
startPoint=pointsXYZ_csRight(1:2,2);
i=0;

while in>0

    endPoint=startPoint-i*step*dir_csRight;
    xSection=endPoint(1,1);
    ySection=endPoint(2,1);

    in=inpolygon(xSection,ySection,xPoly,yPoly);
    
    if in==1
        endPoint=endPoint+step*dir_csRight;
    end
    
    i=i+1;

end

bottomPoint=endPoint;

 bottomBoundaryXY=[x_Bottom;y_Bottom];
    distance=zeros(1,size(bottomBoundaryXY,2));
    
    for n=1:size(bottomBoundaryXY,2)
        distance(1,n)=((bottomBoundaryXY(1,n)-bottomPoint(1,1))^2+(bottomBoundaryXY(2,n)-bottomPoint(2,1))^2)^0.5;
    end
    
    [C,I]=min(distance);
    closestXY=[bottomBoundaryXY(1,I);bottomBoundaryXY(2,I)];

pointsXYZ_csRight(1:2,1)=closestXY;
pointsXYZ_csRight(3,1)=z_Bottom(I);

%ATTACH BOTTOM BOUNDARY, CROSS SECTION LEFT
in=1;
startPoint=pointsXYZ_csLeft(1:2,2);
i=0;

while in>0

    endPoint=startPoint-i*step*dir_csLeft;
    xSection=endPoint(1,1);
    ySection=endPoint(2,1);

    in=inpolygon(xSection,ySection,xPoly,yPoly);
    
    if in==1
        endPoint=endPoint+step*dir_csLeft;
    end
    
    i=i+1;

end

bottomPoint=endPoint;

 bottomBoundaryXY=[x_Bottom;y_Bottom];
    distance=zeros(1,size(bottomBoundaryXY,2));
    
    for n=1:size(bottomBoundaryXY,2)
        distance(1,n)=((bottomBoundaryXY(1,n)-bottomPoint(1,1))^2+(bottomBoundaryXY(2,n)-bottomPoint(2,1))^2)^0.5;
    end
    
    [C,I]=min(distance);
    closestXY=[bottomBoundaryXY(1,I);bottomBoundaryXY(2,I)];

pointsXYZ_csLeft(1:2,1)=closestXY;  
pointsXYZ_csLeft(3,1)=z_Bottom(I);

%ATTACH TOP BOUNDARY, CROSS SECTION LEFT
in=1;
startPoint=pointsXYZ_csLeft(1:2,size(csLeft,2)+1);
i=0;

while in>0

    endPoint=startPoint+i*step*dir_csLeft;
    xSection=endPoint(1,1);
    ySection=endPoint(2,1);

    in=inpolygon(xSection,ySection,xPoly,yPoly);
    
    if in==1
        endPoint=endPoint-step*dir_csLeft;
    end
    
    i=i+1;

end

topPoint=endPoint;

 topBoundaryXY=[x_Top;y_Top];
    distance=zeros(1,size(topBoundaryXY,2));
    
    for n=1:size(topBoundaryXY,2)
        distance(1,n)=((topBoundaryXY(1,n)-topPoint(1,1))^2+(topBoundaryXY(2,n)-topPoint(2,1))^2)^0.5;
    end
    
    [C,I]=min(distance);
    closestXY=[topBoundaryXY(1,I);topBoundaryXY(2,I)];

pointsXYZ_csLeft(1:2,size(csLeft,2)+2)=closestXY;
pointsXYZ_csLeft(3,end)=z_Top(I);

%ATTACH TOP BOUNDARY, CROSS SECTION RIGHT
in=1;
startPoint=pointsXYZ_csRight(1:2,size(csRight,2)+1);
i=0;

while in>0

    endPoint=startPoint+i*step*dir_csRight;
    xSection=endPoint(1,1);
    ySection=endPoint(2,1);

    in=inpolygon(xSection,ySection,xPoly,yPoly);
    
    if in==1
        endPoint=endPoint-step*dir_csRight;
    end
    
    i=i+1;

end

topPoint=endPoint;

 topBoundaryXY=[x_Top;y_Top];
    distance=zeros(1,size(topBoundaryXY,2));
    
    for n=1:size(topBoundaryXY,2)
        distance(1,n)=((topBoundaryXY(1,n)-topPoint(1,1))^2+(topBoundaryXY(2,n)-topPoint(2,1))^2)^0.5;
    end
    
    [C,I]=min(distance);
    closestXY=[topBoundaryXY(1,I);topBoundaryXY(2,I)];

pointsXYZ_csRight(1:2,size(csRight,2)+2)=closestXY;
pointsXYZ_csRight(3,end)=z_Top(I);

%Direction vectors need to be redefined for consistency with the clipping
%region. Until here, the direction vectors where the ones extracted from
%the input data. We have tweaked them a little to match nodes of the
%aquifer boundary.

dir_csRight=normc(pointsXYZ_csRight(1:2,end)-pointsXYZ_csRight(1:2,1));
dir_csLeft=normc(pointsXYZ_csLeft(1:2,end)-pointsXYZ_csLeft(1:2,1));

%====================================================================================================%
%============= UNIFORM DISCRETIZATION -> FIRST REMAP KNOWN TRANSECT DATA TO THIS GRID ===============%
%====================================================================================================%

%This uniform grid will be used to position the deepest part of the valley
%within the model. The known transect data is mapped this this grid. The
%transects are named here as the "extended" versions because they are now
%attached to the aquifer boundary.


%RIGHT CROSS SECTION
bottom_csRightext=pointsXYZ_csRight(1:2,1);
top_csRightext=pointsXYZ_csRight(1:2,size(csRight,2)+2);
dir_csRightext=normc(top_csRightext-bottom_csRightext);

width_csRightext=((bottom_csRightext(1,1)-top_csRightext(1,1))^2+(bottom_csRightext(2,1)-top_csRightext(2,1))^2)^0.5;

tExt_csRight=linspace(0,width_csRightext,numberOfPointsInt);

pointsXYZ_csRightext=zeros(3,numberOfPointsInt);

for i=1:1:numberOfPointsInt
pointsXYZ_csRightext(1:2,i)=bottom_csRightext+tExt_csRight(i)*dir_csRightext;
end

tOri_csRight=zeros(1,size(pointsXYZ_csRight,2));

for i=1:size(pointsXYZ_csRight,2)-1
    tOri_csRight(i+1)=((pointsXYZ_csRight(1,1)-pointsXYZ_csRight(1,i+1))^2+(pointsXYZ_csRight(2,1)-pointsXYZ_csRight(2,i+1))^2)^0.5; 
end

pointsXYZ_csRightext(3,:)=interp1(tOri_csRight,pointsXYZ_csRight(3,:),tExt_csRight);

%LEFT CROSS SECTION
bottom_csLeftext=pointsXYZ_csLeft(1:2,1);
top_csLeftext=pointsXYZ_csLeft(1:2,size(csLeft,2)+2);
dir_csLeftext=normc(top_csLeftext-bottom_csLeftext);

width_csLeftext=((bottom_csLeftext(1,1)-top_csLeftext(1,1))^2+(bottom_csLeftext(2,1)-top_csLeftext(2,1))^2)^0.5;

tExt_csLeft=linspace(0,width_csLeftext,numberOfPointsInt);

pointsXYZ_csLeftext=zeros(3,numberOfPointsInt);

for i=1:1:numberOfPointsInt
pointsXYZ_csLeftext(1:2,i)=bottom_csLeftext+tExt_csLeft(i)*dir_csLeftext;
end

tOri_csLeft=zeros(1,size(pointsXYZ_csLeft,2));

for i=1:size(pointsXYZ_csLeft,2)-1
    tOri_csLeft(i+1)=((pointsXYZ_csLeft(1,1)-pointsXYZ_csLeft(1,i+1))^2+(pointsXYZ_csLeft(2,1)-pointsXYZ_csLeft(2,i+1))^2)^0.5; 
end

pointsXYZ_csLeftext(3,:)=interp1(tOri_csLeft,pointsXYZ_csLeft(3,:),tExt_csLeft);


%====================================================================================================%
%==================== UNIFORM DISCRETIZATION -> CREATE INTERIOR CROSS SECTIONS ======================%
%====================================================================================================%

directionControl=1; 
loc=round(linspace(directionControl+1,size(polyxs,2)-directionControl-1,numberOfSections));
pointsXYZ_int=zeros(2*numberOfSections,numberOfPointsInt);

dirs_int=zeros(2,numberOfSections);

display('defining interior cross-sections...')

for i=1:numberOfSections
    
    %Find closest XY point on channel centerline
    oriXY_csi=[polyxs(1,loc(1,i));polyys(1,loc(1,i))];
    
    %Find the desired normal direction using the directionControl parameter
    centerlineDir=normc([polyxs(loc(1,i)+directionControl)-polyxs(loc(1,i)-directionControl);polyys(loc(1,i)+directionControl)-polyys(loc(1,i)-directionControl)]);
    sectionDir=[-centerlineDir(2,1);abs(-centerlineDir(1,1))];
    dirs_int(:,i)=sectionDir;
    
    %Find pivot point in top boundary 
    step=0.1;
    in=1;
    startPoint=oriXY_csi;
    d=0;

    while in>0

        endPoint=startPoint+d*step*sectionDir;
        xSection=endPoint(1,1);
        ySection=endPoint(2,1);

        in=inpolygon(xSection,ySection,xPoly,yPoly);
        
        if in==1
            endPoint=endPoint-step*sectionDir;
        end
    
        d=d+1;
        
    end
    
    topPoint=endPoint;

    topBoundaryXY=[x_Top;y_Top];
    distance=zeros(1,size(topBoundaryXY,2));
    
    for n=1:size(topBoundaryXY,2)
        
        distance(1,n)=((topBoundaryXY(1,n)-topPoint(1,1))^2+(topBoundaryXY(2,n)-topPoint(2,1))^2)^0.5;
    end
    
    [C,I]=min(distance);
    closestXYtop=[topBoundaryXY(1,I);topBoundaryXY(2,I)];
    
    %Find pivot point in bottom boundary 
    
    step=0.1;
    in=1;
    startPoint=oriXY_csi;
    d=0;

    while in>0

        endPoint=startPoint-d*step*sectionDir;
        xSection=endPoint(1,1);
        ySection=endPoint(2,1);

        in=inpolygon(xSection,ySection,xPoly,yPoly);
    
        d=d+1;

    end
    bottomPoint=endPoint;

    bottomBoundaryXY=[x_Bottom;y_Bottom];
    distance=zeros(1,size(bottomBoundaryXY,2));
    
    for n=1:size(bottomBoundaryXY,2)
        distance(1,n)=((bottomBoundaryXY(1,n)-bottomPoint(1,1))^2+(bottomBoundaryXY(2,n)-bottomPoint(2,1))^2)^0.5;  
    end
    
    [C,I]=min(distance);
    closestXYbottom=[bottomBoundaryXY(1,I);bottomBoundaryXY(2,I)];
    
    %Distance between selected points and corresponding linear space
    
    C=((closestXYbottom(1,1)-closestXYtop(1,1))^2+(closestXYbottom(2,1)-closestXYtop(2,1))^2)^0.5;
    t=linspace(0,C,numberOfPointsInt);
    
    %New direction vector
    
    new_sectionDir=normc(closestXYtop-closestXYbottom);
    dirs_int(:,i)=new_sectionDir;
    
    %Create corresponding cross section
    for j=1:numberOfPointsInt
        
        pointsXYZ_int(2*i-1:2*i,j)=closestXYbottom+t(j)*new_sectionDir;     
        
    end
   
    %i
    
end

%================================================================================================%
%================== DEFINE THE POLIGONAL REGION FOR CLIPPING THE MODEL ==========================%
%================================================================================================%

clipcontrol=0; %allows removing some artifacts in the model border

topLeft=[pointsXYZ_csLeftext(1,end);pointsXYZ_csLeftext(2,end)];
bottomLeft=[pointsXYZ_csLeftext(1,1);pointsXYZ_csLeftext(2,1)];
topRight=[pointsXYZ_csRightext(1,end);pointsXYZ_csRightext(2,end)];
bottomRight=[pointsXYZ_csRightext(1,1);pointsXYZ_csRightext(2,1)];

% find closest point on top/bottom boundary needed to close the polygon

findtopLeft=find(topBoundaryXY(1,:)>topLeft(1,1),1,'first');
findbottomLeft=find(bottomBoundaryXY(1,:)>bottomLeft(1,1),1,'first');
findtopRight=find(topBoundaryXY(1,:)<topRight(1,1),1,'last');
findbottomRight=find(bottomBoundaryXY(1,:)<bottomRight(1,1),1,'last');

top=[topBoundaryXY(1,findtopLeft:findtopRight);topBoundaryXY(2,findtopLeft:findtopRight)];
bottom=[bottomBoundaryXY(1,findbottomLeft:findbottomRight);bottomBoundaryXY(2,findbottomLeft:findbottomRight)];

xPolyClip=[bottomLeft(1,1) topLeft(1,1) top(1,:) topRight(1,1) bottomRight(1,1) fliplr(bottom(1,:)) bottomLeft(1,1)];
yPolyClip=[bottomLeft(2,1) topLeft(2,1) top(2,:) topRight(2,1) bottomRight(2,1) fliplr(bottom(2,:)-clipcontrol) bottomLeft(2,1)];


%==================================================================================================%
%================================== VALLEY MINIMUM PROCEDURES =====================================%
%==================================================================================================%

%LOCATION OF VALLEY BOTTOM ON INTERPOLATED CROSS SECTIONS
[C1,I1]=min(pointsXYZ_csLeftext(3,:));
[C2,I2]=min(pointsXYZ_csRightext(3,:));

%SCAN ALL INTERIOR CROSS SECTIONS, FIND CLOSEST POINT TO MINIMUM AND CHECK
%IF INSIDE CLIPPING POLYGON. SWEEP LR AND RL TO FIND THE FIRST INTERIOR
%CROSS SECTION WITHIN THE MODEL. NOTE THAT WE HAVE CREATED CROSS-SECTIONS
%PAST THE AREA OF INTEREST. THESE WILL BE CLIPPED LATER.

%For the left cross section
d=[];

for i=1:(size(pointsXYZ_int,1))/2
    
    d=[];
    minXY=pointsXYZ_csLeftext(1:2,I1);
    
    for j=1:size(pointsXYZ_int,2)
        d=[d ((pointsXYZ_int(2*i-1,j)-minXY(1,1))^2+(pointsXYZ_int(2*i,j)-minXY(2,1))^2)^0.5];
        [A B]=min(d(1,:));
    end
   
    in=inpolygon(pointsXYZ_int(2*i-1,B),pointsXYZ_int(2*i,B),xPolyClip,yPolyClip);
    minLeftSection=i;
    posLeftSection=B;
    
    if in==1 
        break; end
end

%For the right cross section
d=[];

for i=(size(pointsXYZ_int,1))/2:-1:1
    
    d=[];
    minXY=pointsXYZ_csRightext(1:2,I2);
    
    for j=1:size(pointsXYZ_int,2)
        d=[d ((pointsXYZ_int(2*i-1,j)-minXY(1,1))^2+(pointsXYZ_int(2*i,j)-minXY(2,1))^2)^0.5];
        [A B]=min(d(1,:));
    end
   
    in=inpolygon(pointsXYZ_int(2*i-1,B),pointsXYZ_int(2*i,B),xPolyClip,yPolyClip);
    minRightSection=i;
    posRightSection=B;
    
    if in==1 
        break; end
end

%NUMBER OF SECTIONS WITHIN CLIPPING REGION
newNumberOfSections=minRightSection-minLeftSection;
newlocationAxis=round(linspace(posLeftSection,posRightSection,newNumberOfSections+2));

%MATRIX CONTAINING LOCATION OF MINIMUMS
valleyMins=zeros(2,numberOfSections+2);

%for the known cross sections

valleyMins(:,1)=pointsXYZ_csLeftext(1:2,I1);
valleyMins(:,end)=pointsXYZ_csRightext(1:2,I2);

%FILL THE MATRIX OF MINIMUMS - SPECIAL CONSIDERATION OUTSIDE THE CLIPPING
%REGION

%To the left of clipping region
for i=2:minLeftSection
    valleyMins(:,i)=pointsXYZ_int(2*i-3:2*i-2,posLeftSection); 
end

%To the right of clipping region
for i=minRightSection+2:size(valleyMins,2)-1
    valleyMins(:,i)=pointsXYZ_int(2*i-3:2*i-2,posRightSection); 
end

j=2;

%Whithin the clipping region
for i=minLeftSection+1:minRightSection+1
    valleyMins(:,i)=pointsXYZ_int(2*i-3:2*i-2,newlocationAxis(1,j)); 
    j=j+1;
end


%==================================================================================================%
%=========================== NON-UNIFORM (IRREGULAR) DISCRETIZATION ===============================%
%==================================================================================================%

%NOW WE DISCRETIZE = NUMBER OF POINTS ON EACH SIDE OF THE VALLEY MINIMUM.
%THIS WAY WE CAN PAIR 1-1 AN INTERPOLATE EASILY IN THE AREA OF INTEREST

%CREATE TOP & BOTTOM HALF XYZ VECTORS FOR CROSS SECTION 1
widthBottomHalf_csLeft=((pointsXYZ_csLeftext(1,1)-pointsXYZ_csLeftext(1,I1))^2+(pointsXYZ_csLeftext(2,1)-pointsXYZ_csLeftext(2,I1))^2)^0.5;
widthTopHalf_csLeft=((pointsXYZ_csLeftext(1,end)-pointsXYZ_csLeftext(1,I1))^2+(pointsXYZ_csLeftext(2,end)-pointsXYZ_csLeftext(2,I1))^2)^0.5;

tTop_csLeft=linspace(0,widthTopHalf_csLeft,halfsectionPoints);
tBottom_csLeft=linspace(0,widthBottomHalf_csLeft,halfsectionPoints);

vectorTopHalf_csLeft=zeros(3,halfsectionPoints);
vectorBottomHalf_csLeft=zeros(3,halfsectionPoints);

minXY_csLeft=valleyMins(:,1);

    for i=1:halfsectionPoints
        vectorTopHalf_csLeft(1:2,i)=minXY_csLeft+tTop_csLeft(i)*dir_csLeft;
        vectorBottomHalf_csLeft(1:2,i)=minXY_csLeft-tBottom_csLeft(i)*dir_csLeft;
    end
    
    
%CREATE TOP & BOTTOM HALF XYZ VECTORS FOR CROSS SECTION 2
widthBottomHalf_csRight=((pointsXYZ_csRightext(1,1)-pointsXYZ_csRightext(1,I2))^2+(pointsXYZ_csRightext(2,1)-pointsXYZ_csRightext(2,I2))^2)^0.5;
widthTopHalf_csRight=((pointsXYZ_csRightext(1,end)-pointsXYZ_csRightext(1,I2))^2+(pointsXYZ_csRightext(2,end)-pointsXYZ_csRightext(2,I2))^2)^0.5;

tTop_csRight=linspace(0,widthTopHalf_csRight,halfsectionPoints);
tBottom_csRight=linspace(0,widthBottomHalf_csRight,halfsectionPoints);

vectorTopHalf_csRight=zeros(3,halfsectionPoints);
vectorBottomHalf_csRight=zeros(3,halfsectionPoints);

minXY_csRight=valleyMins(:,end);

    for i=1:halfsectionPoints
        vectorTopHalf_csRight(1:2,i)=minXY_csRight+tTop_csRight(i)*dir_csRight;
        vectorBottomHalf_csRight(1:2,i)=minXY_csRight-tBottom_csRight(i)*dir_csRight;
    end
        
%DISCRETIZE INTERMEDIATE CROSS SECTIONS IN A SIMILAR MANNER
vectorTopHalf_int=zeros(3*numberOfSections,halfsectionPoints);
vectorBottomHalf_int=zeros(3*numberOfSections,halfsectionPoints);

a=2;

for i=1:numberOfSections
    
    if i<minLeftSection
        I=posLeftSection;
    elseif i>minRightSection 
        I=posRightSection; 
    else
        I=newlocationAxis(1,a);
        a=a+1;
    end
    
    widthBottomHalf_int=((pointsXYZ_int(2*i-1,1)-pointsXYZ_int(2*i-1,I))^2+(pointsXYZ_int(2*i,1)-pointsXYZ_int(2*i,I))^2)^0.5;
    widthTopHalf_int=((pointsXYZ_int(2*i-1,end)-pointsXYZ_int(2*i-1,I))^2+(pointsXYZ_int(2*i,end)-pointsXYZ_int(2*i,I))^2)^0.5;

    tTop_int=linspace(0,widthTopHalf_int,halfsectionPoints);
    tBottom_int=linspace(0,widthBottomHalf_int,halfsectionPoints);

    minXY_int=valleyMins(:,i+1);

    for j=1:halfsectionPoints
        vectorTopHalf_int(3*i-2:3*i-1,j)=minXY_int+tTop_int(j)*dirs_int(:,i);
        vectorBottomHalf_int(3*i-2:3*i-1,j)=minXY_int-tBottom_int(j)*dirs_int(:,i);
    end
        
end
    

%==================================================================================================%
%============================ CONSOLIDATE ALL DATA IN FINAL MATRICES ==============================%
%==================================================================================================%

% CREATE RECEIVING MATRICES

interpolationX=zeros(2*halfsectionPoints-1,numberOfSections+2);
interpolationY=zeros(2*halfsectionPoints-1,numberOfSections+2);
interpolationZ=zeros(2*halfsectionPoints-1,numberOfSections+2);

midpoint=round((2*halfsectionPoints-1)/2);

% TRANSFER VALUES TO DESTINATION MATRIX X

interpolationX(1:midpoint,1)=fliplr(vectorTopHalf_csLeft(1,:));
interpolationX(midpoint:end,1)=vectorBottomHalf_csLeft(1,:);

interpolationX(1:midpoint,end)=fliplr(vectorTopHalf_csRight(1,:));
interpolationX(midpoint:end,end)=vectorBottomHalf_csRight(1,:);

for i=1:numberOfSections
  interpolationX(1:midpoint,i+1)=fliplr(vectorTopHalf_int(3*i-2,:)); 
  interpolationX(midpoint:end,i+1)=vectorBottomHalf_int(3*i-2,:); 
end

% TRANSFER VALUES TO DESTINATION MATRIX Y

interpolationY(1:midpoint,1)=fliplr(vectorTopHalf_csLeft(2,:));
interpolationY(midpoint:end,1)=vectorBottomHalf_csLeft(2,:);

interpolationY(1:midpoint,end)=fliplr(vectorTopHalf_csRight(2,:));
interpolationY(midpoint:end,end)=vectorBottomHalf_csRight(2,:);

for i=1:numberOfSections
  interpolationY(1:midpoint,i+1)=fliplr(vectorTopHalf_int(3*i-1,:)); 
  interpolationY(midpoint:end,i+1)=vectorBottomHalf_int(3*i-1,:); 
end

% USE CUBIC SPLINES TO LOCATE Z VALUES ON CROSS SECTIONS 1 AND 2
tnew_csLeft=zeros(1,2*halfsectionPoints-1);
tnew_csRight=zeros(1,2*halfsectionPoints-1);

tnew_csLeft(1,1:halfsectionPoints)=tBottom_csLeft;
tnew_csLeft(1,halfsectionPoints+1:end)=tTop_csLeft(1,2:end)+tnew_csLeft(1,halfsectionPoints);

tnew_csRight(1,1:halfsectionPoints)=tBottom_csRight;
tnew_csRight(1,halfsectionPoints+1:end)=tTop_csRight(1,2:end)+tnew_csRight(1,halfsectionPoints);

interpolationZ(:,1)=fliplr(interp1(tOri_csLeft,pointsXYZ_csLeft(3,:),tnew_csLeft,'cubic'));
interpolationZ(:,numberOfSections+2)=fliplr(interp1(tOri_csRight,pointsXYZ_csRight(3,:),tnew_csRight,'cubic'));


%APPLY INPOLYGON TO OBTAIN A MASK OF NODES INSIDE CLIPPING REGION
mask_node=inpolygon(interpolationX,interpolationY,xPolyClip,yPolyClip);

% DISTANCE MATRICES, LEFT AND RIGHT
distanceL=zeros(2*halfsectionPoints-1,numberOfSections);
distanceR=zeros(2*halfsectionPoints-1,numberOfSections);


display('calculating distances...')


%CALCULATE DISTANCES TO THE TRANSECT ON THE LEFT

for k=1:2*halfsectionPoints-1

    for i=1:numberOfSections
    
        for j=1:i-1
        
          if mask_node(k,j+1)==0 && mask_node(k,j+2)==0 
            
              case_point=1;
          
          elseif mask_node(k,j+1)==0 && mask_node(k,j+2)==1 
        
              case_point=2;
              
          elseif mask_node(k,j+1)==1 && mask_node(k,j+2)==1
            
              case_point=3;
          
          elseif mask_node(k,j+1)==1 && mask_node(k,j+2)==0
            
              case_point=4;
            
          end
        
          
        if case_point==1
            distanceL(k,i)=distanceL(k,i)+0;
        end
        
        if case_point==2  
            distanceL(k,i)=distanceL(k,i)+((interpolationY(k,j+2)-interpolationY(k,1)).^2+(interpolationX(k,j+2)-interpolationX(k,1)).^2).^0.5; 
        end
        
        if case_point==3  
            distanceL(k,i)=distanceL(k,i)+((interpolationY(k,j+2)-interpolationY(k,j+1)).^2+(interpolationX(k,j+2)-interpolationX(k,j+1)).^2).^0.5; 
        end
        
        if case_point==4 
            distanceL(k,i)=distanceL(k,i)+((interpolationY(k,j+1)-interpolationY(k,end)).^2+(interpolationX(k,j+1)-interpolationX(k,end)).^2).^0.5; 
            break
        end
        
        end
        
    end
    
end



%CALCULATE DISTANCES TO THE TRANSECT ON THE RIGHT


for k=1:2*halfsectionPoints-1

    for i=numberOfSections:-1:1
        
        for j=numberOfSections:-1:i+1
        
          if mask_node(k,j)==0 && mask_node(k,j+1)==0 
            
              case_point=1;
          
          elseif mask_node(k,j)==1 && mask_node(k,j+1)==0 
        
              case_point=2;
              
          elseif mask_node(k,j)==1 && mask_node(k,j+1)==1
            
              case_point=3;
          
          elseif mask_node(k,j)==0 && mask_node(k,j+1)==1
            
              case_point=4;
            
          end
        
          
        if case_point==1
            distanceR(k,i)=distanceR(k,i)+0;
        end
        
        if case_point==2  
            distanceR(k,i)=distanceR(k,i)+((interpolationY(k,j)-interpolationY(k,end)).^2+(interpolationX(k,j)-interpolationX(k,end)).^2).^0.5; 
        end
        
        if case_point==3  
            distanceR(k,i)=distanceR(k,i)+((interpolationY(k,j)-interpolationY(k,j+1)).^2+(interpolationX(k,j)-interpolationX(k,j+1)).^2).^0.5; 
        end
        
        if case_point==4 
            distanceR(k,i)=distanceR(k,i)+((interpolationY(k,j+1)-interpolationY(k,1)).^2+(interpolationX(k,j+1)-interpolationX(k,1)).^2).^0.5; 
            break
        end
        
        end
        
    end
    
end


%Calculate percentage matrix for interpolation

distanceLR=zeros(2*halfsectionPoints-1,numberOfSections);

distanceLR=distanceL+distanceR;
distanceLPercent=1-distanceL./distanceLR;
distanceRPercent=1-distanceR./distanceLR;


%Mask the interpolation Z matrix

for i=1:size(interpolationZ,1)
    for j=1:size(interpolationZ,2)
        
        if mask_node(i,j)==0
            interpolationZ(i,j)=NaN;
        end
        
    end
end

%==================================================================================================%
%=========================== FINAL INTERPOLATION & ORGANIZE OUTPUT ================================%
%==================================================================================================%

p=5;
xk=0;
yk=0;
W=0.05;

for j=1:numberOfSections
    
    for i=1:size(interpolationZ,1)
    
    interpolationZ(i,j+1)=fn_IDWknick_interp(p,distanceL(i,j),xk,yk,distanceLR(i,j),interpolationZ(i,1),interpolationZ(i,end),W);
  
    end

end

%ELIMINATE FIRST AND LAST COLUMNS FOR CONSISTENCY (REMOVE KNOWN SECTIONS)

% interpolationX=interpolationX(:,2:end-1);
% interpolationY=interpolationY(:,2:end-1);
% interpolationZ=interpolationZ(:,2:end-1);

% CREATE FINAL VECTORS WITH RESULTS 

x=[];
y=[];
z=[];

for i=1:size(interpolationX,2)
x=[x interpolationX(:,i)'];    
end 

for i=1:size(interpolationY,2)
y=[y interpolationY(:,i)'];    
end 

for i=1:size(interpolationZ,2)
z=[z interpolationZ(:,i)'];    
end 

%==================================================================================================%
%============================= STORE INTERPOLATION IN A REGULAR GRID ==============================%
%==================================================================================================%

display('gridding...')

[maxY,I]=max(top(2,:));
[minY,I]=min(bottom(2,:));

[maxXbottom,I]=max(bottom(1,:));
[minXbottom,I]=min(bottom(1,:));

[maxXtop,I]=max(top(1,:));
[minXtop,I]=min(top(1,:));

[maxX,I]=max([maxXtop maxXbottom]);
[minX,I]=min([minXtop minXbottom]);

delta=0.05;
resolution=0.01;
[xx,yy]=meshgrid(minX-delta:resolution:maxX+delta,minY-delta:resolution:maxY+delta);
ZI=griddata(x(:,2:end-1),y(:,2:end-1),z(:,2:end-1),xx,yy);

display('clipping the area of interest...')

synthetic_mask_grid

%==============================================%
%============ LOAD AND APPLY MASK =============%
%==============================================%

%When doing many realization this mask changes so it is better to calculate
%one mask for each case.

% mask=load('mask_delta1.5_res0.01.mat','mask');
% mask=mask.mask;
% ZI=mask.*ZI;

