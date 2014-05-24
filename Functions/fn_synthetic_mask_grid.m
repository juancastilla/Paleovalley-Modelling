%CLIPS MODEL TO AREA INSIDE CLIPPING POLYGON - NO OUTPUT MASK

%JCC 16072013


function [ZI_clipped]=fn_synthetic_mask_grid(ZI,xx,yy,xPolyClip,yPolyClip)

i=0;j=0;in=0;

in=inpolygon(xx,yy,xPolyClip,yPolyClip);

display('finishing clipping procedures...')

for i=1:size(xx,1)
     
    for j=1:size(xx,2)
      
        if in(i,j)==0
            ZI(i,j)=NaN;

        end 

    end

end

ZI_clipped=ZI;

end