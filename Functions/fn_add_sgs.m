
function ZI_realizations=fn_add_sgs(ZI,S,nrealz)

ZI_realizations=cell(1,nrealz);

    for n=1:nrealz
        ZI_realizations{1,nrealz}=ZI+S.D(:,:,n)';     
    end
    
end