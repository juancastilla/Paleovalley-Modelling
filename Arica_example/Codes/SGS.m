%% CALL SGSIM, SET MODEL SIZE

S=sgems_get_par('sgsim');

S.dim.nx=size(ZI,1);
S.dim.ny=size(ZI,2);

S.XML.parameters.Nb_Realizations.value=3;


%% SET VARIOGRAM MODE

S.XML.parameters.Variogram=sgems_variogram_xml('68 Gau(50)');


%% SET CONDITIONING DATA FOR PALEOCHANNEL BORDERS

S.d_obs=[];
findnan=isnan(ZI);

for i=1:size(ZI,2)

    idtop=find(ZI(:,i)==0,1,'first');
    idbottom=find(ZI(:,i)==0,1,'last');
    
    if isempty(idtop)==0
    S.d_obs=[S.d_obs;idtop i 0 ZI(idtop,i)];
    end
    
    if isempty(idbottom)==0
    S.d_obs=[S.d_obs;idbottom i 0 ZI(idbottom,i)];
    end
   
end

%% SET CONDITIONING DATA FOR PALEOCHANNEL MINIMUM

conditionBandPixels=0;
check=0;

%Keep paleochannel boundaries fixed

for i=1:size(ZI,2)

   [C,I]=min(ZI(:,i));
   
   check=isnan(C);
   
   if check~=1
   S.d_obs=[S.d_obs;I i 0 0];
   end
        
end

%% SET CONDITIONING DATA FOR WELLS (ONLY FOR SEGMENT DC)

%load wells data 
wells=importdata('Wells Taylor.txt')';
wells=wells.data;

%extract elevations from trend surface

d=zeros(size(xx,1),size(xx,2));
collect=zeros(size(wells,1),4);

for n=1:size(wells,1)   
    for i=1:size(xx,1)
        for j=1:size(xx,2)
            dist=((wells(n,1)-xx(i,j))^2+(wells(n,2)-(yy(i,j)))^2)^0.5;
            d(i,j)=dist;
        end
    end
    
    [minval,ind]=min(d(:)); 
    [Y,Z]=ind2sub([size(d,1) size(d,2)],ind);
    
    collect(n,:)=[Y Z 0 wells(n,3)-ZI(Y,Z)];
    
    d=[];
    
end

S.d_obs=[S.d_obs;collect];


%% RUN SGSIM

S=sgems_grid(S);


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

%% AQUIFER BOUDARIES

zboundtop=ones(1,size(x_Top,2))*1000;
zboundbottom=ones(1,size(x_Bottom,2))*1000;

figure
hold on
plot3(x_Top,y_Top,zboundtop,'-k', 'LineWidth',1)
plot3(x_Bottom,y_Bottom,zboundbottom,'-k','LineWidth',1)
axis equal tight
%% ADD SGSIM NOISE & PLOT

figure
ZIN=ZI+S.D(:,:,1,1);
surf(xx,yy,ZIN);
view(0,90)
shading flat

figure
surf(xx,yy,ZI);
view(0,90)
shading flat

%% PLOT NOISE

for i=1:S.XML.parameters.Nb_Realizations.value;
  subplot(4,3,i);
  imagesc(S.x,S.y,S.D(:,:,1,i));
end

%% PLOT REALIZATIONS

noiseIntensity=1;

for i=1:S.XML.parameters.Nb_Realizations.value;
  subplot(1,3,i);
  ZIN=ZI+S.D(:,:,1,i)*noiseIntensity;
  hold on
  hold on
plot(x_Top,y_Top,'-k', 'LineWidth',1)
plot(x_Bottom,y_Bottom,'-k','LineWidth',1)
axis equal tight
  surf(xx,yy,ZIN);
  view(0,90)
  shading flat
  axis equal tight
end
