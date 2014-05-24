%=========================================================================%
%================================ DESCRIPTION ============================%
%=========================================================================%

%THIS FUNCTION GENERATES SEVERAL TREND RECONSTRUCTIONS FOR DIFFERENT
%COMBINATIONS OF N CONDITIONING TRANSECTS. HERE WE ONLY FOCUS ON THE
%UNCERTAINTY OF THE TREND.

%JCC 18072013 06082013


function [cs,ZI_realz,volumes_realizations]=fn_realizations_iterative(sampling_sites,sections_inside,sampling_res,control_length,boundary_left,boundary_right,mask,mesh_info,num_repetitions,threshold)


%=======================================================================================%
%============================= LOAD DATA TO FUNCTION WORKSPACE =========================%
%=======================================================================================%

load('base_surface_CS1mod_CS4.mat');


%=====================================================================================================%
%==================== SELECT SAMPLING TRANSECTS AND EXTRACT Z DATA FROM REFERENCE ====================%
%=====================================================================================================%

display(['selecting reconstruction data for N=' num2str(sampling_sites)])

%We randomly select n=sampling_sites for the reconstructions. The output of
%this function is a matrix, and the rows are different combinations of
%sampling sites for a fixed n.

%Make a list of possible sampling combinations, filtered to remove repeated
%ones
[sampling_combinations]=fn_sampling_design(sections_inside,sampling_sites,num_repetitions,threshold);

%NOTE: sampling_combinations is a matrix. The rows indicate a combination of
%cross-sections that will be used to generate a realizations and calculate its volumean.

display(['begin reconstructions for combinations of N=' num2str(sampling_sites)])
    
%NOTE: here, s iterates over the possible sampling combinations:

for s=1:size(sampling_combinations,1)
    
    display(['working on reconstruction ' num2str(s) ' of ' num2str(size(sampling_combinations,1)) '(N=' num2str(sampling_sites) ')'])
    
    sections_selected=sampling_combinations(s,:);
    
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

    display(['reconstruction data for N=' num2str(sampling_sites) ' created!'])

    
    cs{s}=cs_all;   %copy list of cross-sections to variable cs for output
    
    %=====================================================================================================%
    %================================ GENERATE TREND SURFACE PIECE BY PIECE ==============================%
    %=====================================================================================================%

    display(['generating piecewise trend surface for N=' num2str(sampling_sites)])

    parfor i=1:size(cs_all_full,3)-1
    
        %First define the piece we are reconstructing
        csLeft=cs_all_full(:,:,i);
        csRight=cs_all_full(:,:,i+1);
    
        %Run the simulations and obtain a cell structure with the data for the
        %different pieces. We need the ZI matrix, xx and yy. The index "i"
        %denotes the number of the piece
        [ZI_piece{i},xx_piece{i},yy_piece{i}]=fn_synthetic_sim_trend(csLeft,csRight);
    
    end

    display(['piecewise trend surface for N=' num2str(sampling_sites) ' finalized!'])

    
    %==================================================================%
    %======================= STITCH PIECES ============================%
    %==================================================================%

    display(['stitching pieces for N=' num2str(sampling_sites) ' transects...'])
    
    [ZI_stitched{s}]=fn_stitch_pieces(ZI_piece,xx_piece,yy_piece,mask,mesh_info);
    
    display(['stitching for N=' num2str(sampling_sites) ' finalized!'])
    
    ZI_realz{s}=ZI_stitched{s};
    

    %=====================================================================================================%
    %=============================== ADD SGS SURFACE AND OBTAIN REALIZATIONS =============================%
    %=====================================================================================================%
    
    %     display(['SGS realizations for n=' num2str(sampling_sites) ' transects...'])
    %
    %     [S]=fn_synthetic_sim_sgs(nrealz,ZI_stitched,sill,range);
    %
    %     for i=1:nrealz
    %
    %         [ZI_realz{i}]=ZI_stitched+S.D(:,:,nrealz)';
    %
    %     end
    %
    %     display(['SGS realizations for n=' num2str(sampling_sites) ' finalized!'])
    
    %=====================================================================================================%
    %=================================== CALCULATE VOLUMES FOR OUTPUT ====================================%
    %=====================================================================================================%
    
    display(['calculating volumes for N=' num2str(sampling_sites) ' transects...'])
    
    [volumes_realizations{s}]=fn_volume(ZI_realz{s},boundary_left(1,3),boundary_right(1,3));
        
    display(['volume calculations for N=' num2str(sampling_sites) ' finalized!'])
    
end