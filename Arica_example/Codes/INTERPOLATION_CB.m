%% CLEAR MEMORY

clear all

%% IMPORT DATA & ORGANIZE

top=fliplr(importdata('topboundary.txt'))';
bottom=fliplr(importdata('bottomboundary.txt'))';
csA=importdata('cs1.txt')';
csB=importdata('cs3.txt')';
csC=importdata('cs4.txt')';
csD=importdata('cs5.txt')';
csE=importdata('cscity.txt')';
top(1:2,:)=top(2:-1:1,:);
bottom(1:2,:)=bottom(2:-1:1,:);

% NOTE: Cross sections should be specified in a South-North direction

%% BEGIN DEFINE GLOBAL MODEL VARIABLES & CHOOSE PAIR OF KNOWN CROSS SECTIONS TO USE 

csLeft=csC;
csRight=csB;

numberOfPointsInt=100;
numberOfSections=40;
halfsectionPoints=50;

dir_csLeft=normc([csLeft(1,end)-csLeft(1,1);abs(csLeft(2,1)-csLeft(2,end))]);
dir_csRight=normc([csRight(1,end)-csRight(1,1);abs(csRight(2,1)-csRight(2,end))]);

% VALLEY BOUNDARY MODEL EXTENTS X DIRECTION
[Cmin,Imin]=min(top(1,:));
[Cmax,Imax]=max(top(1,:));

minX_Top=Cmin;
maxX_Top=Cmax;

[Cmin,Imin]=min(bottom(1,:));
[Cmax,Imax]=max(bottom(1,:));

minX_Bottom=Cmin;
maxX_Bottom=Cmax;

% XYZ VECTORS FOR TOP AND BOTTOM BOUNDARIES CONSTRUCTED FROM GIS DATA
x_Top=top(1,:);
y_Top=top(2,:);
z_Top=top(3,:);

x_Bottom=bottom(1,:);
y_Bottom=bottom(2,:);
z_Bottom=bottom(3,:);

n=10000;

% TOP BOUNDARY
Xtop=[(1:1:size(x_Top,2));x_Top];
Ytop=[(1:1:size(y_Top,2));y_Top];
Ztop=[(1:1:size(z_Top,2));z_Top];

XtopInterp=interp1(Xtop(1,:),Xtop(2,:),linspace(1,Xtop(1,end),n));
YtopInterp=interp1(Ytop(1,:),Ytop(2,:),linspace(1,Ytop(1,end),n));
ZtopInterp=interp1(Ztop(1,:),Ztop(2,:),linspace(1,Ztop(1,end),n));

% BOTTOM BOUNDARY
Xbottom=[(1:1:size(x_Bottom,2));x_Bottom];
Ybottom=[(1:1:size(y_Bottom,2));y_Bottom];
Zbottom=[(1:1:size(z_Bottom,2));z_Bottom];

XbottomInterp=interp1(Xbottom(1,:),Xbottom(2,:),linspace(1,Xbottom(1,end),n));
YbottomInterp=interp1(Ybottom(1,:),Ybottom(2,:),linspace(1,Ybottom(1,end),n));
ZbottomInterp=interp1(Zbottom(1,:),Zbottom(2,:),linspace(1,Zbottom(1,end),n));

% REASSIGN INTERPOLATED VALUES
x_Top=XtopInterp;
y_Top=YtopInterp;
z_Top=ZtopInterp;

x_Bottom=XbottomInterp;
y_Bottom=YbottomInterp;
z_Bottom=ZbottomInterp;

% PREPARE VECTOR DEFINING POLYGONAL REGION>
xPoly=[x_Top fliplr(x_Bottom) x_Top(1,1)];
yPoly=[y_Top fliplr(y_Bottom) y_Top(1,1)];

%% CHANNEL CENTERLINE BINARY IMAGE

% FIND POINTS IN POLYGON
crust=3000;
dmesh=50;
[maxY,I]=max(csLeft(2,:));
[minY,I]=min(csRight(2,:));

[maxX,I]=max(csRight(1,:));
[minX,I]=min(csLeft(1,:));
[xx,yy]=meshgrid(minX-crust:dmesh:maxX+crust,minY-crust:dmesh:maxY+crust);

[IN] = inpolygon(xx(:),yy(:),xPoly,yPoly);
IN=reshape(IN,size(xx,1),size(yy,2));

% VISUALIZE
% figure;clf
% surf(xx,yy,double(IN))
% view(0,90)
% axis equal tight

%% APPLY SKELETON ALGORITHM

% PREPROCESS
S1=double(IN);

% DO SKELETON
S1=bwmorph(S1,'skel',inf);
S2=bwmorph(S1,'spur',25);

% CROP CRUST
%[xx,yy]=meshgrid(minX:dmesh:maxX,minY:dmesh:maxY);
%crustpixels=crust/dmesh;
crustpixels=10;
S2=S2(crustpixels:end-crustpixels,crustpixels:end-crustpixels);
xx2=xx(crustpixels:end-crustpixels,crustpixels:end-crustpixels);
yy2=yy(crustpixels:end-crustpixels,crustpixels:end-crustpixels);

% VISUALIZE
clf
subplot(2,1,1)
surf(xx,yy,double(S1))
view(0,90)
axis equal tight
subplot(2,1,2)
surf(xx2,yy2,double(S2))
view(0,90)
axis equal tight

%CREATE AN XY VECTOR CONTAINING THE CHANNEL CENTERLINE
polyx=xx2(S2==1);
polyy=yy2(S2==1);

polyxs=polyx(1:10:end);
polyys=polyy(1:10:end);

A=linspace(2,numel(polyxs),10000);
polyxs=interp1(1:numel(polyxs),polyxs,A,'spline');
polyys=interp1(1:numel(polyys),polyys,A,'spline');

%clf
%plot(polyxs,polyys,'xk')

%% PROCESS KNOWN DATA, ATTACH TO AQUIFER BOUNDARIES

% MATRICES FOR KNOWN CROSS SECTIONS 1 & 2
pointsXYZ_csLeft=zeros(3,size(csLeft,2)+2);
pointsXYZ_csRight=zeros(3,size(csRight,2)+2);

pointsXYZ_csLeft(1:3,2:size(csLeft,2)+1)=csLeft(1:3,:);
pointsXYZ_csRight(1:3,2:size(csRight,2)+1)=csRight(1:3,:);

%Attach cross sections to boundaries by including closest points to the 
%pointsXYZ_ matrix in order to achieve this, we need to extend the known 
%cross sections as close to the boundaries as possible (using inpolygon)
%and then select the closest value

%ATTACH BOTTOM BOUNDARY, CROSS SECTION RIGHT   
step=1;
in=1;
startPoint=pointsXYZ_csRight(1:2,2);
i=0;

while in>0

    endPoint=startPoint-i*step*dir_csRight;
    xSection=endPoint(1,1);
    ySection=endPoint(2,1);

    in=inpolygon(xSection,ySection,xPoly,yPoly);
    
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

%% PLOT BOUNDARIES AND CENTERLINE, USER DEFINES CENTERLINE EXTENTS
 
clf
plot(XtopInterp,YtopInterp,'-k')
hold
plot(XbottomInterp,YbottomInterp,'-k')
plot(pointsXYZ_csLeft(1,:),pointsXYZ_csLeft(2,:),'xr')
plot(pointsXYZ_csRight(1,:),pointsXYZ_csRight(2,:),'xr')

reduce=0;

polyxs=polyxs(1,1+reduce:end-reduce);
polyys=polyys(1,1+reduce:end-reduce);

plot(polyxs,polyys,'-b')
axis equal tight

%% UNIFORM DISCRETIZATION - KNOWN CROSS SECTIONS

%THIS UNIFORM DISCRETIZATION ALLOWS TO LOCATE THE VALLEY MINIMUM ALONG THE 
%STUDIED CROSS SECTION. DEFINE "EXTENDED" VERSIONS OF KNOWN CROSS SECTIONS

%CROSS SECTION RIGHT
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

%CROSS SECTION LEFT
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

%% UNIFORM DISCRETIZATION - CREATE INTERIOR CROSS SECTIONS

directionControl=1;
loc=round(linspace(directionControl+1,size(polyxs,2)-directionControl-1,numberOfSections));
pointsXYZ_int=zeros(2*numberOfSections,numberOfPointsInt);

dirs_int=zeros(2,numberOfSections);

h=0;

for i=1:numberOfSections
    
    %Find closest XY point on channel centerline
    oriXY_csi=[polyxs(1,loc(1,i));polyys(1,loc(1,i))];
    
    %Find the desired normal direction using the directionControl parameter
    centerlineDir=normc([polyxs(loc(1,i)+directionControl)-polyxs(loc(1,i)-directionControl);polyys(loc(1,i)+directionControl)-polyys(loc(1,i)-directionControl)]);
    sectionDir=[-centerlineDir(2,1);abs(-centerlineDir(1,1))];
    dirs_int(:,i)=sectionDir;
    
    %Find pivot point in top boundary 
    step=10;
    in=1;
    startPoint=oriXY_csi;
    d=0;

    while in>0

        endPoint=startPoint+d*step*sectionDir;
        xSection=endPoint(1,1);
        ySection=endPoint(2,1);

        in=inpolygon(xSection,ySection,xPoly,yPoly);
    
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
    
    step=10;
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
    
    %Create corresponding cross section
    for j=1:numberOfPointsInt
        
        pointsXYZ_int(2*i-1:2*i,j)=closestXYbottom+t(j)*sectionDir;     
        
    end
    
 h=h+1
   
end


%% CLIPPING REGION

clipcontrol=0;

topLeft=[pointsXYZ_csLeft(1,end);pointsXYZ_csLeft(2,end)];
bottomLeft=[pointsXYZ_csLeft(1,1);pointsXYZ_csLeft(2,1)];
topRight=[pointsXYZ_csRight(1,end);pointsXYZ_csRight(2,end)];
bottomRight=[pointsXYZ_csRight(1,1);pointsXYZ_csRight(2,1)];

% find closest point on top/bottom boundary needed to close the polygon

findtopLeft=find(topBoundaryXY(1,:)>topLeft(1,1),1,'first');
findbottomLeft=find(bottomBoundaryXY(1,:)>bottomLeft(1,1),1,'first');
findtopRight=find(topBoundaryXY(1,:)<topRight(1,1),1,'last');
findbottomRight=find(bottomBoundaryXY(1,:)<bottomRight(1,1),1,'last');

top=[topBoundaryXY(1,findtopLeft:findtopRight);topBoundaryXY(2,findtopLeft:findtopRight)];
bottom=[bottomBoundaryXY(1,findbottomLeft:findbottomRight);bottomBoundaryXY(2,findbottomLeft:findbottomRight)];

xPolyClip=[bottomLeft(1,1) topLeft(1,1) top(1,:) topRight(1,1) bottomRight(1,1) fliplr(bottom(1,:)) bottomLeft(1,1)];
yPolyClip=[bottomLeft(2,1) topLeft(2,1) top(2,:) topRight(2,1) bottomRight(2,1) fliplr(bottom(2,:)-clipcontrol) bottomLeft(2,1)];

% plot(xPolyClip,yPolyClip)

%% VALLEY MINIMUM PROCEDURES

%LOCATION OF VALLEY BOTTOM ON INTERPOLATED CROSS SECTIONS
[C1,I1]=min(pointsXYZ_csLeftext(3,:));
[C2,I2]=min(pointsXYZ_csRightext(3,:));

%SCAN ALL INTERIOR CROSS SECTIONS, FIND CLOSEST POINT TO MINIMUM AND CHECK IF INSIDE CLIPPING POLYGON

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

%% NEW DISCRETIZATION

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
    
%% CONSOLIDATE ALL DATA IN FINAL MATRICES

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

% DISTANCE MATRICES, LEFT AND RIGHT
distanceL=zeros(2*halfsectionPoints-1,numberOfSections);
distanceR=zeros(2*halfsectionPoints-1,numberOfSections);

% DISTANCE L MATRIX
xDiffSqL=zeros(2*halfsectionPoints-1,numberOfSections);   %(x-xi)^2
yDiffSqL=zeros(2*halfsectionPoints-1,numberOfSections);   %(y-yi)^2

% Calculate (x-xi)^2 matrix

for i=1:size(interpolationX,2)-2
    xDiffSqL(:,i)=(interpolationX(:,i+1)-interpolationX(:,1)).^2; 
end

% Calculate (y-yi)^2 matrix

for i=1:size(interpolationY,2)-2
    yDiffSqL(:,i)=(interpolationY(:,i+1)-interpolationY(:,1)).^2; 
end

distanceL=(xDiffSqL+yDiffSqL).^0.5;


% DISTANCE R MATRIX

xDiffSqR=zeros(2*halfsectionPoints-1,numberOfSections);   %(x-xi)^2
yDiffSqR=zeros(2*halfsectionPoints-1,numberOfSections);   %(y-yi)^2

% Calculate (x-xi)^2 matrix

for i=1:size(interpolationX,2)-2
    xDiffSqR(:,i)=(interpolationX(:,i+1)-interpolationX(:,numberOfSections+2)).^2; 
end

% Calculate (y-yi)^2 matrix

for i=1:size(interpolationY,2)-2
    yDiffSqR(:,i)=(interpolationY(:,i+1)-interpolationY(:,numberOfSections+2)).^2; 
end

distanceR=(xDiffSqR+yDiffSqR).^0.5;

% Calculate percentage matrix for interpolation

distanceLR=zeros(2*halfsectionPoints-1,numberOfSections);

distanceLR=distanceL+distanceR;
distanceLPercent=1-distanceL./distanceLR;
distanceRPercent=1-distanceR./distanceLR;

%FINAL INTERPOLATION
for i=1:numberOfSections
    interpolationZ(:,i+1)=interpolationZ(:,1).*distanceLPercent(:,i)+interpolationZ(:,numberOfSections+2).*distanceRPercent(:,i);   
end

%% ELIMINATE FIRST AND LAST COLUMNS FOR CONSISTENCY (REMOVE KNOWN SECTIONS)

interpolationX=interpolationX(:,2:end-1);
interpolationY=interpolationY(:,2:end-1);
interpolationZ=interpolationZ(:,2:end-1);

%% FINAL VECTORS WITH RESULTS 

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

%_________________________________________________________________________________%

%% PLOT AQUIFER BOUNDARIES

figure
hold on

plot(x_Top,y_Top,'-k')
plot(x_Bottom,y_Bottom,'-k')

axis equal tight

%% PLOT CLIPPING POLYGON

plot(xPolyClip,yPolyClip)

%% PLOT KNOWN CROSS SECTIONS - ORIGINAL DATA

plot(pointsXYZ_csLeft(1,:),pointsXYZ_csLeft(2,:),'ok',pointsXYZ_csRight(1,:),pointsXYZ_csRight(2,:),'ok')
axis equal tight

%% PLOT KNOWN CROSS SECTIONS - NEW DISCRETIZATION

plot(pointsXYZ_csLeftext(1,:),pointsXYZ_csLeftext(2,:),'-r')
plot(pointsXYZ_csRightext(1,:),pointsXYZ_csRightext(2,:),'-r')
axis equal tight

%% PLOT INTERIOR CROSS SECTIONS

for i=1:numberOfSections
plot(vectorTopHalf_int(3*i-2,:),vectorTopHalf_int(3*i-1,:),'-b')
plot(vectorBottomHalf_int(3*i-2,:),vectorBottomHalf_int(3*i-1,:),'-b')
end
axis equal tight

%% PLOT LOCATION OF VALLEY MINIMUM

for i=1:numberOfSections
    plot(valleyMins(1,2:end-1),valleyMins(2,2:end-1),'-r','LineWidth',1.2,'MarkerSize',40)
end
axis equal tight

%% PLOT INTERPOLATED SECTIONS WITH 3D LINES

figure
plot3(x,y,z,'xr')
axis tight

%% 2D PLOT - EXTRACT INTERPOLATED CROSS SECTION 

distance=zeros(1,2*halfsectionPoints-1);
extract=10; %number of cross section to extract

for i=1:2*halfsectionPoints-1
    
    d=((interpolationX(1,extract)-(interpolationX(i,extract)))^2+(interpolationY(1,extract)-(interpolationY(i,extract)))^2)^0.5;
    distance(1,i)=d;
    
end

plot(distance,interpolationZ(:,extract))

%% PLOT VALLEY MINIMUM DEPTH

l=[];

for i=1:size(interpolationZ,2)
   [C I]=min(interpolationZ(:,i));
   l=[l C];
end

plot(l)


%% MORPH MOVIE

[a b]=max(l);
[c d]=min(l);

distance=[];
figure

for j=b:d
    
    extract=j; %number of cross section to extract

    for i=1:2*halfsectionPoints-1
   
        d=((interpolationX(1,extract)-(interpolationX(i,extract)))^2+(interpolationY(1,extract)-(interpolationY(i,extract)))^2)^0.5;
        distance(1,i)=d;
    
end

plot(distance,interpolationZ(:,extract));
axis ([0 2 -150 0])
M(j) = getframe;

end

%% PLAY MOVIE

axis ([0 2 -150 0])
movie(M)

%% CREATE AVI MOVIE

movie2avi(M,morph.avi);


%% STORE INTERPOLATION IN A REGULAR GRID (GRIDDATA)

[maxY,I]=max(top(2,:));
[minY,I]=min(bottom(2,:));

[maxX,I]=max(top(1,:));
[minX,I]=min(bottom(1,:));

delta=10;
resolution=5;
[xx,yy]=meshgrid(minX-delta:resolution:maxX+delta,minY-delta:resolution:maxY+delta);
ZI=griddata(x(:,2:end-1),y(:,2:end-1),z(:,2:end-1),xx,yy);


%% CLIP MODEL TO AREA INSIDE CLIPPING POLYGON

i=0;j=0;in=0;h=0;
counter=size(ZI,1)*size(ZI,2);

for i=1:size(xx,1)
     
    for j=1:size(xx,2)
        in=inpolygon(xx(i,j),yy(i,j),xPolyClip,yPolyClip);
        counter-h
        h=h+1;
    
        if in==0
            ZI(i,j)=NaN;
        end 
end

end


%% AQUIFER BOUDARIES

%figure
hold on
plot(x_Top,y_Top,'-k', 'LineWidth',1)
plot(x_Bottom,y_Bottom,'-k','LineWidth',1)
axis equal tight

%% PLOT FINAL CLIPPED SURFACE

surf(xx,yy,ZI/1);
view(0,90)
shading flat