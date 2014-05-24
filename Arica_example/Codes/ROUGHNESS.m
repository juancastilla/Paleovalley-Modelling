%EXTRACT SAMPLE VARIOGRAM OF TOPOGRAPHIC FEATURES FROM DEM DATA

%% IMPORT DATA

clear all

R1=importdata('Roughness 1.txt');
R2=importdata('Roughness 2.txt');

%Choose roughness dataset to extract

R=R1;

%% PLOT DEM DATA

dmesh=10;
[maxY,I]=max(R(:,2));
[minY,I]=min(R(:,2));
[maxX,I]=max(R(:,1));
[minX,I]=min(R(:,1));

[xx,yy]=meshgrid(minX:dmesh:maxX,minY:dmesh:maxY);

RI=griddata(R(:,1),R(:,2),R(:,3),xx,yy,'linear');

%% PLOT DEM ELEVATIONS

figure
surf(xx,yy,RI);
view(0,90)
axis equal
shading flat

%% MOVING AVERAGE PROCEDURE

window=30;
SmoothedSurface=MovAvg(RI,window);

%% PLOT MOVING AVERAGE SURFACE

figure
surf(xx,yy,SmoothedSurface);
view(0,90)
axis equal
shading flat

%% RESIDUAL SURFACE

Residual=RI-SmoothedSurface;

figure
surf(xx,yy,Residual);
view(0,90)
axis equal
shading flat

%% PREPARE DATA FOR EXPORT TO SGEMS

%Reshape the grid into a single row

G=reshape(Residual',[],1);   

%Display max/min values to construct SGEMS cartesian grid

maxY
minY
maxX
minX

cellsx=size(RI,2)
cellsy=size(RI,1)

Gnan=isnan(G);

for i=1:size(G,1)
    
    if Gnan(i,1)==1
        G(i,1)=-9999;
    end
    
end

% save -ascii G.txt G


