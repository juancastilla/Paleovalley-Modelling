%CLIPS MODEL TO AREA INSIDE CLIPPING POLYGON - WITH OUTPUT MASK

%JCC 16072013


function [ZI_clipped,mask]=fn_synthetic_mask_grid2(ZI,xx,yy,xPolyClip,yPolyClip)

i=0;j=0;in=0;h=0;
mask=ones(size(ZI,1),size(ZI,2));

in=inpolygon(xx,yy,xPolyClip,yPolyClip);

display('finishing clipping procedures...')

for i=1:size(xx,1)
     
    for j=1:size(xx,2)
      
        if in(i,j)==0
            ZI(i,j)=NaN;
            mask(i,j)=NaN;   %creates a mask that we can save and output
        end 
       
    end
end

ZI_clipped=ZI;

end