%Generates a cross section between two given points p1(x1,y1) and p2(x2,y2) 
%with a specified resolution. NOTE: p2 has to be the northernmost point.
%Output is a two row vector of the form: 
%[x1 xa xb xc....x2
% y1 ya yb yc....y2]

%JCC 11072013

function [cs]=fn_crosssectiongen(p1,p2,res)

spaceX=linspace(p1(1,1),p2(1,1),res);
spaceY=linspace(p1(1,2),p2(1,2),res);

cs=[spaceX;spaceY];

end