%=====================================================================%
%============================ DESCRIPTION ============================%
%=====================================================================%

%DESIGNS A SUITE OF SAMPLING COMBINATIONS FOR DIFFERENT NUMBER OF
%CONDITIONING TRANSECT DATA. THIS ENSEMBLE OF COMBINATIONS IS USED TO
%GENERATE MULTIPLE TREND REALIZATIONS AND CHARACTERIZE THE ASSOCIATED
%VOLUME UNCERTAINTY.

%JCC 18072013 08072013


function [sampling_combinations]=fn_sampling_design(sections_inside,num_sections,num_repetitions,threshold)

sampling_combinations=[];

pool=0;

if num_sections==size(sections_inside,2)
    
    sampling_combinations=sections_inside;
    
elseif num_sections<threshold

    while pool<num_repetitions

        combination=randi(size(sections_inside,2),1,num_sections);
        combination_sorted=unique(combination);
    
        if size(combination_sorted,2)==size(combination,2)
        
            sampling_combinations=[sampling_combinations;combination_sorted];
            pool=pool+1;
      
        end
 
    end

else
    
    while pool<num_repetitions

        combination=randi(size(sections_inside,2),1,size(sections_inside,2)-num_sections);
        combination_sorted=unique(combination);
    
        if size(combination_sorted,2)==size(combination,2)
        
            temp=sections_inside;
            temp(combination_sorted)=[];
            sampling_combinations=[sampling_combinations;temp];
            pool=pool+1;
      
        end
 
    end
    
end  
  
sampling_combinations=unique(sampling_combinations,'rows');

%correction for cross-section indexing

if num_sections<threshold
    
    sampling_combinations=sampling_combinations+16;
    
end


end