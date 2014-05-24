
%% CALCULATE AQUIFER VOLUMES FOR ALL REALIZATIONS

for i=1:size(ZI_realz,2)
    
    for j=1:size(ZI_realz{1,i},2)
        
        volumes_realizations_new{1,i}{1,j}=fn_volume(ZI_realz{1,i}{1,j},boundary_left(1,3),boundary_right(1,3));
    
    end
    
end

%% ENSEMBLE STATISTICS (FIGURE 3 LAB NOTEBOOK)

s1=1;
s2=2;
s3=4;
s4=9;

selected_s=[s1 s2 s3 s4];

for i=1:size(selected_s,2)
    
    s=selected_s(1,i);
    ensemble=[];

    for j=1:size(ZI_realz{1,s},2)
        
        ensemble=cat(3,ensemble,ZI_realz{1,s}{1,j});    %creates 3D ensemble array
        ensemble_stats{1,i}{1,1}=mean(ensemble,3);      %calculates ensemble mean  
        ensemble_stats{1,i}{1,2}=std(ensemble,0,3);     %calculates ensemble standard deviation
        %ensemble_stats{1,i}{1,3}=;                     %calculates ensemble range
        
    end
    
end


%%