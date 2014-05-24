%=========================================================================%
%================================ DESCRIPTION ============================%
%=========================================================================%

%GENERATES A MASK USING THE INFORMATION EXTRACTED FROM THE BASE REFERENCE
%SURFACE GRID, AND THE LATERAL TRANSECTS THAT CONSTRAIN THE MODEL. THIS
%FUNCTION IS USEFUL TO CLIP REALIZATIONS FASTER BY MULTIPLYING THIS MASK
%WITH A RECONSTRUCTION GRID.

%JCC 18072013



function [mask]=fn_generate_mask(boundary_left,boundary_right,sampling_res,mesh_info)


%=======================================================================================%
%============================= LOAD DATA TO FUNCTION WORKSPACE =========================%
%=======================================================================================%

load('base_surface_CS1mod_CS4.mat');


%=====================================================================================================%
%================= SELECT LATERAL SAMPLING TRANSECTS AND EXTRACT Z DATA FROM REFERENCE ===============%
%=====================================================================================================%

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

cs_all=[x1_L y1_L x2_L y2_L;x1_R y1_R x2_R y2_R];

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


%=====================================================================================================%
%============= GENERATE TREND SURFACE USING LATERAL SECTIONS TO OBTAIN A CLIPPING MASK ===============%
%=====================================================================================================%

display(['working on the clipping mask for the model...'])

%Choose the lateral cross-sections
csLeft=cs_all_full(:,:,1);
csRight=cs_all_full(:,:,end);
    
[~,~,~,mask]=fn_synthetic_sim_trend2(csLeft,csRight,mesh_info);

display(['clipping mask finalized...'])

%NOTE: We apply this mask within the stitch procedure to clean all the realizations later on