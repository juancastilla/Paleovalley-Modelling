cerro=importdata('Cerro Data.txt');

dmesh=0.1;
[maxY,I]=max(cerro(:,3));
[minY,I]=min(cerro(:,3));
[maxX,I]=max(cerro(:,2));
[minX,I]=min(cerro(:,2));

[xx,yy]=meshgrid(minX:dmesh:maxX,minY:dmesh:maxY);

ZI=griddata(cerro(:,2),cerro(:,3),cerro(:,1),xx,yy);

surf(xx,yy,ZI);
view(0,90)
shading flat