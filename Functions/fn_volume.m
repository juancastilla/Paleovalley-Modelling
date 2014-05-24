function [vol]=fn_volume(ZIN,boundary_left,boundary_right)

    model=ZIN(:,boundary_left:boundary_right);
    model(isnan(model)==1)=0;
    vol=abs(sum(sum(model)));

end
