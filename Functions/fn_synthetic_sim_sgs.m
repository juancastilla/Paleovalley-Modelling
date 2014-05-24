function [S]=fn_synthetic_sim_sgs(nrealz,ZI_stitched,sill,range)


%=======================================%
%===== CALL SGSIM, SET MODEL SIZE ======%
%=======================================%

S=sgems_get_par('sgsim');

S.dim.nx=size(ZI_stitched,2);
S.dim.ny=size(ZI_stitched,1);

S.XML.parameters.Nb_Realizations.value=nrealz;


%=======================================%
%========= SET VARIOGRAM MODE ==========%
%=======================================%

variogram=[num2str(sill) ' Gau(' num2str(range) ')'];
S.XML.parameters.Variogram=sgems_variogram_xml(variogram);


%=====================================================================%
%=========== SET CONDITIONING DATA FOR PALEOCHANNEL BORDERS ==========%
%=====================================================================%

S.d_obs=[];
findnan=round(double(isnan(ZI_stitched)));
findnan(findnan==1)=2;
findnan(findnan==0)=1;
findnan(findnan==2)=0;

for i=1:size(ZI_stitched,2)

    idtop=find(findnan(:,i),1,'first');
    
    if isempty(idtop)==0
    %S.d_obs=[S.d_obs;i idtop 0 ZI_stitched(idtop,i)];
    S.d_obs=[S.d_obs;i idtop 0 0];
    end
    
end

for i=1:size(ZI_stitched,2)

    idbottom=find(findnan(:,i),1,'last');
    
    if isempty(idbottom)==0
    %S.d_obs=[S.d_obs;i idbottom 0 ZI_stitched(idbottom,i)];
    S.d_obs=[S.d_obs;i idbottom 0 0];
    end
    
end


%============================================================================%
%============ SET CONDITIONING DATA FOR PALEOCHANNEL MINIMUM ================%
%============================================================================%

% conditionBandPixels=0;
% check=0;
% 
% %Keep paleochannel boundaries fixed
% 
% for i=1:size(ZI,2)
% 
%    [C,I]=min(ZI(:,i));
%    
%    check=isnan(C);
%    
%    if check~=1
%    S.d_obs=[S.d_obs;I i 0 0];
%    end
%         
% end


%==================================%
%=========== RUN SGSIM ============%
%==================================%

S=sgems_grid(S);

end
