function [ZI_stitched]=fn_stitch_pieces(ZI_piece,xx_piece,yy_piece,mask,mesh_info)
   
    x=[];
    y=[];
    z=[];
    
    vx_reduced=[];
    vy_reduced=[];
    vz_reduced=[];

for i=1:size(ZI_piece,2)
    
    %Remove the NaNs from the piece. This allows us to stitch the pieces
    %correctly
   
    %ZI_piece{i}(isnan(ZI_piece{i})==1)=0;
   
    %Reshape grids into vectors 
    vx=reshape(xx_piece{i}(:,:),1,[]);
    vy=reshape(yy_piece{i}(:,:),1,[]); 
    vz=reshape(ZI_piece{i}(:,:),1,[]); 
    
    %Remove any cell in the data that is a NaN
    index=isnan(vz);
    
    for j=1:size(index,2)
    
        if index(1,j)~=1
            vx_reduced=[vx_reduced vx(1,j)];
            vy_reduced=[vy_reduced vy(1,j)];
            vz_reduced=[vz_reduced vz(1,j)];
        end
    
    end
    
    %Collect the reduced piece vectors in a single array
    x=[x vx_reduced];
    y=[y vy_reduced];
    z=[z vz_reduced];
   
end

%Grid the merged data

minX=mesh_info(1,1);
maxX=mesh_info(1,2);
minY=mesh_info(1,3);
maxY=mesh_info(1,4);
resolution=mesh_info(1,5);

[xx,yy]=meshgrid(minX:resolution:maxX,minY:resolution:maxY);

ZI_stitched=griddata(x,y,z,xx,yy);

ZI_stitched=mask.*ZI_stitched;

end
