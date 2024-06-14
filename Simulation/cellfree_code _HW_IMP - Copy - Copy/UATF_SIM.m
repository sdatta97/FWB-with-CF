function SINR = UATF_SIM(V_sq,Delta,V,add_INT,n,UE)
SINR=zeros(2,2,UE);
for ue=1:UE
    SINR(:,:,ue) =abs(mean(V_sq(:,:,ue,n,:),5)).^2*pinv(eye(2)+mean(Delta(:,:,ue,n,:),5)+mean(V(:,:,ue,n,:),5)-abs(mean(V_sq(:,:,ue,n,:),5)).^2+mean(add_INT(:,:,ue,n,:),5));
end