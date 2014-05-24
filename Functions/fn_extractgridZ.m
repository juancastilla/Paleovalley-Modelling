%Extracts the Z value from a grid generated by meshgrid, given a desired
%(x,y) location

%JCC 11072013

function [z]=fn_extractgridZ(x,y,xx,yy,ZIN)

x_xx_diff=abs(xx(1,:)-x);
y_yy_diff=abs(yy(:,1)-y);

[Cx,Ix]=min(x_xx_diff);
[Cy,Iy]=min(y_yy_diff);

z=ZIN(Iy,Ix);

end