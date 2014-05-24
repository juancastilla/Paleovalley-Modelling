
function [cs_all,ZI_realz,volumes_realizations]=fn_realizations_uniform(nrealz,sampling_sites,sections_inside,sampling_res,sill,range,control_length,boundary_left,boundary_right,mask,mesh_info)


%=======================================================================================%
%============================= LOAD DATA TO FUNCTION WORKSPACE =========================%
%=======================================================================================%

load('base_surface_CS1mod_CS4.mat');


%=====================================================================================================%
%==================== SELECT SAMPLING TRANSECTS AND EXTRACT Z DATA FROM REFERENCE ====================%
%=====================================================================================================%

display(['selecting reconstruction data for n=' num2str(sampling_sites) ' transects...'])

%We uniformly select a subset of these cross-sections for sampling

if sampling_sites<size(sections_inside,2)

    vector_sections=round(linspace(sections_inside(1,1),sections_inside(1,end),sampling_sites+2));
    sections_selected=vector_sections(1,2:end-1);
  
else
    
    %If we want to use all the cross sections we don't need to calculate 
    sections_selected=sections_inside;
    
end
    
%Take the first and last points that define the selected sampling
%cross-sections. We put this information in the cs_all matrix. 
%Format here is [x1 y1 x2 y2]

cs_all=[];

for i=1:size(sections_selected,2)
    
    x1=interpolationX(size(interpolationX,1)-control_length,sections_selected(1,i));
    y1=interpolationY(size(interpolationX,1)-control_length,sections_selected(1,i));
    x2=interpolationX(control_length,sections_selected(1,i));
    y2=interpolationY(control_length,sections_selected(1,i));
    
    cs_all=[cs_all;x1 y1 x2 y2];
    
end

%Add leftmost and rightmost sampling sites. These are needed to constrain
%the model appropriately. By inspecting ZIN we choose appropriate
%boundaries for the model

x1_L=xx(1,boundary_left(1,3));
y1_L=yy(boundary_left(1,1),1);
x2_L=xx(1,boundary_left(1,3));
y2_L=yy(boundary_left(1,2),1);

x1_R=xx(1,boundary_right(1,3));
y1_R=yy(boundary_right(1,1),1);
x2_R=xx(1,boundary_right(1,3));
y2_R=yy(boundary_right(1,2),1);

cs_all=[x1_L y1_L x2_L y2_L;cs_all;x1_R y1_R x2_R y2_R];

%Define points that will be used for sampling within the cross section.
%First row are "x" values, second row "y", third row "z" values. We will
%place the z values in the next step extracting them from the base surface

cs_all_full=zeros(3,sampling_res,size(cs_all,1));

for i=1:size(cs_all,1)
    cs_all_full(1:2,:,i)=fn_crosssectiongen(cs_all(i,1:2),cs_all(i,3:4),sampling_res);
end


%This loop samples de reference surface in each point and extracts the depth
%to bedrock

for i=1:size(cs_all,1)
    for j=1:sampling_res
        cs_all_full(3,j,i)=fn_extractgridZ(cs_all_full(1,j,i),cs_all_full(2,j,i),xx,yy,ZIN);
    end
end

display(['reconstruction data for n=' num2str(sampling_sites) ' created!'])

%=====================================================================================================%
%================================ GENERATE TREND SURFACE PIECE BY PIECE ==============================%
%=====================================================================================================%

display(['generating piecewise trend surface for n=' num2str(sampling_sites) ' transects...'])

%matlabpool

parfor i=1:size(cs_all_full,3)-1
    
    %First define the piece we are reconstructing
    csLeft=cs_all_full(:,:,i);
    csRight=cs_all_full(:,:,i+1);
    
    %Run the simulations and obtain a cell structure with the data for the
    %different pieces. We need the ZI matrix, xx and yy. The index "i"
    %denotes the number of the piece
    [ZI_piece{i},xx_piece{i},yy_piece{i}]=fn_synthetic_sim_trend(csLeft,csRight);
    
end

%matlabpool close

display(['piecewise trend surface for n=' num2str(sampling_sites) ' finalized!'])

%==================================================================%
%======================= STITCH PIECES ============================%
%==================================================================%

display(['stitching pieces for n=' num2str(sampling_sites) ' transects...'])

[ZI_stitched]=fn_stitch_pieces(ZI_piece,xx_piece,yy_piece,mask,mesh_info);

display(['stitching for n=' num2str(sampling_sites) ' finalized!'])


%=====================================================================================================%
%=============================== ADD SGS SURFACE AND OBTAIN REALIZATIONS =============================%
%=====================================================================================================%

display(['SGS realizations for n=' num2str(sampling_sites) ' transects...'])

[S]=fn_synthetic_sim_sgs(nrealz,ZI_stitched,sill,range);

for i=1:nrealz
    
    [ZI_realz{i}]=ZI_stitched+S.D(:,:,nrealz)';

end

display(['SGS realizations for n=' num2str(sampling_sites) ' finalized!'])

%=====================================================================================================%
%=================================== CALCULATE VOLUMES FOR OUTPUT ====================================%
%=====================================================================================================%

display(['calculating volumes for n=' num2str(sampling_sites) ' transects...'])

volumes_realizations=[];

for i=1:nrealz

[volumes_realizations]=[volumes_realizations fn_volume(ZI_realz{i},boundary_left(1,3),boundary_right(1,3))];

end

display(['volume calculations for n=' num2str(sampling_sites) ' finalized!'])

end