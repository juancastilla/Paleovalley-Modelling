function z=fn_IDWknick_interp(p,x,xk,yk,DT,zL,zR,W)

z=[];

%x corresponds to the value of DL
%xk is entered as a percentage = xk*DT is the location of the knickpoint

if x==0
    z=zL;
end

if x==DT
    z=zR;
end

% First interval
if x<xk*DT && x~=0
z=(W*zL*(x)^-p+(1-W)*zR*(DT-x)^-p)/(W*x^-p+(1-W)*(DT-x)^-p)+x*yk/(xk*DT);
end

if x>=xk*DT && x~=DT 
z=(W*zL*(x)^-p+(1-W)*zR*(DT-x)^-p)/(W*x^-p+(1-W)*(DT-x)^-p)+((DT*yk)/(DT-xk*DT))-((x*yk)/(DT-xk*DT));
end

end