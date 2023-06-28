function p = pLoS4(BS_locations, UE_location, theta,omega,psi, idx_max)
M_eff = size(idx_max,1);
M = size(BS_locations,1);
r = zeros(M_eff,1);
for ap = 1:M_eff
    r(ap) = sqrt(sum((UE_location-BS_locations(idx_max(ap),:)).^2)); 
end
[~, bsPriorities] = sort(r,'ascend');
% for ap = 1:M
%     prob(ap) = find(bsPriorities==ap); 
% end
% prob = (1./prob)./sum(1./prob);
% prob = (1./(r.^2))/sum((1./(r.^2)));
prob = (1./r)/sum((1./r));
% prob = (theta.*r)./(theta.*r+omega);
% prob = (1-prob)./sum(1-prob);
r_avg = r'*prob;
% p_lb = (theta*r(1)/(theta*r(1)+psi))^M + (theta*r(1)/(2*omega))*M*((psi*omega)/(theta*r(1)*(theta*r(1)+psi+omega)))*((theta*r(1)/(theta*r(1)+omega))*(1+(omega/(theta*r(1)+psi))))^M+0.5*M*(theta*r(1)/(theta*r(1)+psi))^M*(psi/(theta*r(1)+omega));
% p_ub = (theta*r(1)/(theta*r(1)+psi))^M + (theta*r(1)/(2*omega))*2*((psi*omega)/(theta*r(1)*(theta*r(1)+psi+omega)))^2*((theta*r(1)/(theta*r(1)+omega))*(1+(omega/(theta*r(1)+psi))))^2+0.5*M*(theta*r(1)/(theta*r(1)+psi))^M*(psi/(theta*r(1)+omega));
% p_lb = (theta*r(1)/(theta*r(1)+psi))^M + (theta*r(1)/(2*omega))*M*((psi*omega)/((theta*r(1)+psi)*(theta*r(1)+omega)))*((theta*r(1)/(theta*r(1)+omega))*(1+(omega/(theta*r(1)+psi))))^M+0.5*M*(theta*r(1)/(theta*r(1)+psi))^M*(psi/(theta*r(1)+omega));
% p_ub = (theta*r(1)/(theta*r(1)+psi))^M + (theta*r(1)/(2*omega))*2*((psi*omega)/((theta*r(1)+psi)*(theta*r(1)+omega)))*((theta*r(1)/(theta*r(1)+omega))*(1+(omega/(theta*r(1)+psi))))+0.5*2*(theta*r(1)/(theta*r(1)+psi))^2*(psi/(theta*r(1)+omega));
p_ub = (theta*r_avg/(theta*r_avg+psi))^M + (theta*r_avg/(2*omega))*2*((psi*omega)/((theta*r_avg+psi)*(theta*r_avg+omega)))*((theta*r_avg/(theta*r_avg+omega))*(1+(omega/(theta*r_avg+psi))))+0.5*2*(theta*r_avg/(theta*r_avg+psi))^2*(psi/(theta*r_avg+omega));
% p = 0.5*(p_lb+p_ub);
p = p_ub;
% p_act = (theta*r(1)/(theta*r(1)+psi))^M + (theta*r(1)/omega)*2*((psi*omega)/(theta*r(1)*(theta*r(1)+psi+omega)))^2*((theta*r(1)/(theta*r(1)+omega))*(1+(omega/(theta*r(1)+psi))))^2;
% p_act = (theta*r(1)/(theta*r(1)+psi))^M + (theta*r(1)/omega)*M*((psi*omega)/(theta*r(1)*(theta*r(1)+psi+omega)))*((theta*r(1)/(theta*r(1)+omega))*(1+(omega/(theta*r(1)+psi))))^M;
% for i=2:M
%     p_act = p_act - nchoosek(M,i)*(theta*r(1)/(theta*r(1)+psi))^(M-i)*(psi/(theta*r(1)+psi))^i*(theta*r(1)/(theta*r(1)+omega))^i;
% end
% p = p_act;
end