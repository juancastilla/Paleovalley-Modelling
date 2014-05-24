%ALLOWS CREATING REFERENCE BASE SURFACES WITH A KNICKPOINT AND BEDROCK
%SCOUR HIGHS AND LOWS. BOTH MORPHOLOGICAL FEATURES ARE PARAMETRIZED AND
%SEVERAL BASE SURFACES CAN BE OBTAINED SIMULTANEOUSLY. 

%OUTPUTS THE NECESSARY VARIABLES TO CHARACTERIZE UNCERTAINTY:
%ZIN(interpolated surface including noise),xx,yy,interpolationX,interpolationY,mask_node

%NOTES:
%-OUTPUTS ARE CONTAINED IN CELLS
%-ZIN IS SELECTED FROM ONE OF THE ZI_realz (ONE OF THE REALIZATIONS)
%-ZI_realz IS ITSELF A CELL THAT CONTAINS SEVERAL COLUMNS WITH REALIZATIONS
%-FOR LOOPS SHOULD BE PARALLELIZED FOR MORE COMBINATIONS OF PARAMETERS.

%JCC 15072013 18072013

%=====================================================================================================%
%======================================= SIMULATION PARAMETERS =======================================%
%=====================================================================================================%

% W=linspace(0.001,0.999,12);
% p=linspace(1,8,8);
% sill=linspace(5,40,8);

W=[0.50];
p=5;
sill=10;
range=40;
nrealz=1;
numberOfSections=100;

%=====================================================================================================%
%======================================= GENERATE TREND SURFACE ======================================%
%=====================================================================================================%

%matlabpool

for k=1:size(sill,2)

    for j=1:size(p,2)
    
        for i=1:size(W,2)
                    
            %Create trend surfaces with each W and p
            [ZI{i,j,k},xx{i,j,k},yy{i,j,k},interpolationX{i,j,k},interpolationY{i,j,k},mask_node{i,j,k}]=fn_createbasesurface(p(1,j),W(1,i),numberOfSections);
   
            %Run SGS
            [S{i,j,k}]=fn_synthetic_sim_sgs(nrealz,ZI{i,j,k},sill(1,k),range);
      
            %This is the array that stores the realizations for each
            %reference surface considering combinations of parameters
            ZI_realz{i,j,k}=fn_add_sgs(ZI{i,j,k},S{i,j,k},nrealz);
    
        end

    end

end
    
%matlabpool close


