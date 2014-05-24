%% IMPORT WORKSPACES AND PLOT FULL TREND MODEL

xxCB=load('/Users/juancarloscastilla/Dropbox/00 PhD/05 Reconstruction of Paleochannel Geomorphology - Paper/Matlab/Real Data/Workspaces/Stitched/xxCB.mat');
xxDC=load('/Users/juancarloscastilla/Dropbox/00 PhD/05 Reconstruction of Paleochannel Geomorphology - Paper/Matlab/Real Data/Workspaces/Stitched/xxDC.mat');
xxBA=load('/Users/juancarloscastilla/Dropbox/00 PhD/05 Reconstruction of Paleochannel Geomorphology - Paper/Matlab/Real Data/Workspaces/Stitched/xxBA.mat');
xxED=load('/Users/juancarloscastilla/Dropbox/00 PhD/05 Reconstruction of Paleochannel Geomorphology - Paper/Matlab/Real Data/Workspaces/Stitched/xxED.mat');

yyCB=load('/Users/juancarloscastilla/Dropbox/00 PhD/05 Reconstruction of Paleochannel Geomorphology - Paper/Matlab/Real Data/Workspaces/Stitched/yyCB.mat');
yyDC=load('/Users/juancarloscastilla/Dropbox/00 PhD/05 Reconstruction of Paleochannel Geomorphology - Paper/Matlab/Real Data/Workspaces/Stitched/yyDC.mat');
yyBA=load('/Users/juancarloscastilla/Dropbox/00 PhD/05 Reconstruction of Paleochannel Geomorphology - Paper/Matlab/Real Data/Workspaces/Stitched/yyBA.mat');
yyED=load('/Users/juancarloscastilla/Dropbox/00 PhD/05 Reconstruction of Paleochannel Geomorphology - Paper/Matlab/Real Data/Workspaces/Stitched/yyED.mat');


ZICB=load('/Users/juancarloscastilla/Dropbox/00 PhD/05 Reconstruction of Paleochannel Geomorphology - Paper/Matlab/Real Data/Workspaces/Stitched/ZICB.mat');
ZIDC=load('/Users/juancarloscastilla/Dropbox/00 PhD/05 Reconstruction of Paleochannel Geomorphology - Paper/Matlab/Real Data/Workspaces/Stitched/ZIDC.mat');
ZIBA=load('/Users/juancarloscastilla/Dropbox/00 PhD/05 Reconstruction of Paleochannel Geomorphology - Paper/Matlab/Real Data/Workspaces/Stitched/ZIBA.mat');
ZIED=load('/Users/juancarloscastilla/Dropbox/00 PhD/05 Reconstruction of Paleochannel Geomorphology - Paper/Matlab/Real Data/Workspaces/Stitched/ZIED.mat');


xxCB=xxCB.xx;
xxDC=xxDC.xx;
xxBA=xxBA.xx;
xxED=xxED.xx;

yyCB=yyCB.yy;
yyDC=yyDC.yy;
yyBA=yyBA.yy;
yyED=yyED.yy;

ZICB=ZICB.ZI;
ZIDC=ZIDC.ZI;
ZIBA=ZIBA.ZI;
ZIED=ZIED.ZI;

hold on 
surf(xxDC,yyDC,ZIDC);
surf(xxCB,yyCB,ZICB);
surf(xxBA,yyBA,ZIBA);
surf(xxED,yyED,ZIED);

view(0,90)
shading flat


%% RESHAPE xx,yy,ZI

xBA=reshape(xxBA,1,[]);
yBA=reshape(yyBA,1,[]);
zBA=reshape(ZIBA,1,[]);

xCB=reshape(xxCB,1,[]);
yCB=reshape(yyCB,1,[]);
zCB=reshape(ZICB,1,[]);

xDC=reshape(xxDC,1,[]);
yDC=reshape(yyDC,1,[]);
zDC=reshape(ZIDC,1,[]);

xED=reshape(xxED,1,[]);
yED=reshape(yyED,1,[]);
zED=reshape(ZIED,1,[]);

indexBA=isnan(zBA);
indexCB=isnan(zCB);
indexDC=isnan(zDC);
indexED=isnan(zED);

xBAreduced=[];
yBAreduced=[];
zBAreduced=[];

xCBreduced=[];
yCBreduced=[];
zCBreduced=[];

xDCreduced=[];
yDCreduced=[];
zDCreduced=[];

xEDreduced=[];
yEDreduced=[];
zEDreduced=[];
    
for i=1:size(indexBA,2)
    
    if indexBA(1,i)~=1
        xBAreduced=[xBAreduced xBA(1,i)];
        yBAreduced=[yBAreduced yBA(1,i)];
        zBAreduced=[zBAreduced zBA(1,i)];
        size(indexBA,2)-i
    end
    
end

for i=1:size(indexCB,2)
    
    if indexCB(1,i)~=1
        xCBreduced=[xCBreduced xCB(1,i)];
        yCBreduced=[yCBreduced yCB(1,i)];
        zCBreduced=[zCBreduced zCB(1,i)];
        size(indexCB,2)-i
    end
    
end

for i=1:size(indexDC,2)
    
    if indexDC(1,i)~=1
        xDCreduced=[xDCreduced xDC(1,i)];
        yDCreduced=[yDCreduced yDC(1,i)];
        zDCreduced=[zDCreduced zDC(1,i)];
        size(indexDC,2)-i
    end
    
end

for i=1:size(indexED,2)
    
    if indexED(1,i)~=1
        xEDreduced=[xEDreduced xED(1,i)];
        yEDreduced=[yEDreduced yED(1,i)];
        zEDreduced=[zEDreduced zED(1,i)];
        size(indexED,2)-i
    end
    
end
%%
x=[];
y=[];
z=[];

x=[xBAreduced xCBreduced xDCreduced xEDreduced];
y=[yBAreduced yCBreduced yDCreduced yEDreduced];
z=[zBAreduced zCBreduced zDCreduced zEDreduced];

x=[xCBreduced xDCreduced];
y=[yCBreduced yDCreduced];
z=[zCBreduced zDCreduced];

%% STORE INTERPOLATION IN A REGULAR GRID (GRIDDATA)

maxY=7960000;
minY=7942000;

maxX=395000;
minX=358000;

delta=0;
resolution=20;
[xx,yy]=meshgrid(minX-delta:resolution:maxX+delta,minY-delta:resolution:maxY+delta);
ZI=griddata(x,y,z,xx,yy);

%% LOAD DATA NEEDED FOR PLOTS

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

csLeft=csB;
csRight=csA;

numberOfPointsInt=1000;
numberOfSections=40;
halfsectionPoints=500;

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

%% PLOTS

zboundtop=ones(1,size(x_Top,2))*1000;
zboundbottom=ones(1,size(x_Bottom,2))*1000;

figure
hold on
plot3(x_Top,y_Top,zboundtop,'-k', 'LineWidth',1)
plot3(x_Bottom,y_Bottom,zboundbottom,'-k','LineWidth',1)
axis equal tight

surf(xx,yy,ZI);
view(0,90)
shading flat