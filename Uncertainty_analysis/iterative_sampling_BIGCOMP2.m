%=========================================================================%
%================================ DESCRIPTION ============================%
%=========================================================================%

%UNCERTAINTY CHARACTERIZATION - ITERATIVE SAMPLING

%THE LOGIC:

%(1) Load data from a reference surface. This is deemed to be the "reality"
%(2) Extract a series of cross-sections and core data depicting the depth
%    to bedrock
%(3) Build a paleovalley model with different subsets of this input data
%(4) Calculate the volume of each realization. Obtain probability
%    distribution functions
%(5) Compare these distributions with the volume of the real surface

%JCC 18072013 06082013



%======================================================================================%
%=================================== PARAMETERS =======================================%
%======================================================================================%

sampling_res=200;           %density of cross section sampling of the real base surface
control_length=10;         %controls the length of sampling sections
num_repetitions=30;        %total location combinations for each number of transects
threshold=35;              %threshold value for the sampling_design function 

%NOTE: no SGS is done here as we are interested in characterizing the
%uncertainty of the trend reconstruction.

%Housekeeping. We need to remove some artifacts left behind by the
%inpolygon procedure that are located in the left and right boundaries of
%the model. For this, we inspect the matrix ZIN and locate the first and
%last columns where this artifact dissapears. Usually it is around 3 grid
%nodes wide. The sections below are the leftmost and rightmost sampling 
%transects that we will use. Coordinates for these points are obtained from 
%xx(:,section_left) and yy(:,section_right)
%boundary_nnn=[rowtop rowbottom column]

boundary_left=[20 200 9];
boundary_right=[20 200 1203];


%=================================================%
%============ READ BASE SURFACE SURFACE ==========%
%=================================================%

load('base_surface_CS1mod_CS4.mat');

%Extract meshgrid info
%format:[minX maxX minY maxY resolution]
mesh_info=[xx(1,1) xx(1,end) yy(1,1) yy(end,1) (xx(1,2)-xx(1,1))];


%======================================================%
%======== CALCULATE VOLUME OF REFERENCE SURFACE =======%
%======================================================%

[volume_reference]=fn_volume(ZIN,boundary_left(1,3),boundary_right(1,3));


%======================================================%
%====== BUILD THE CLIPPING MASK FOR REALIZATIONS ======%
%======================================================%

[mask]=fn_generate_mask(boundary_left,boundary_right,sampling_res,mesh_info);


%=========================================================================%
%================ IDENTIFY CROSS SECTIONS 100% WITHIN MODEL ==============%
%=========================================================================%

%Find those cross sections that lie completely within the clipping region 
%(i.e the area we want to model). We have this data in the mask_node matrix 
%when calculating the reference trend surface. These cross-sectiosn where
%produced from the skeletonization centerline 

mask_node_double=double(mask_node);
mask_node_double(1,48)=1;                                       %fixes problem in this location
sections_inside=[];

for i=2:size(mask_node_double,2)-1
    
    if sum(mask_node_double(:,i)) == size(mask_node_double,1)   %cross-section is fully inside the clipping region
        case_section=1; 
        sections_inside=[sections_inside i];
    else                                                        %cross-section is fully or partially outside 
        case_section=0; 
    end
   
end

%===========================================================%
%=================== SELECT SIMULATIONS ====================%
%===========================================================%

simulations=1:3:25;

%===========================================================%
%===================== RUN SIMULATIONS =====================%
%===========================================================%

%Run the simulations, the result is a vector of volumes, surfaces and other 
%information relevant for each realization for each ensemble of realizations
%that correspond to the different sampling densities. 

%Open a local parallel computing session
matlabpool 10

%IMPORTANT NOTE: ZI_realz is a cell array containing multiple realizations!!!

parfor i=1:size(simulations,2)
    
    j=simulations(1,i);  %read the simulation parameter for this iteration
    
    display(['Calculating realizations for N=' num2str(j) ' transects...'])
    
    [cs_all{i},ZI_realz{i},volumes_realizations{i}]=fn_realizations_iterative(j,sections_inside,sampling_res,control_length,boundary_left,boundary_right,mask,mesh_info,num_repetitions,threshold);
    
    display(['REALIZATIONS FOR N=' num2str(j) ' TRANSECTS FINALIZED!'])
    
end

matlabpool close

save('iterative_D200_R30_1:3:25_070813_IMAC','-v7.3')











