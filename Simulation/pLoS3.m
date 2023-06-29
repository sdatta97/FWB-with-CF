% function p = pLoS3(theta,omega,R)
function [p,t] = pLoS3(BS_locations, UE_location, theta,omega,psi, idx_max)
M_eff = size(idx_max,1);
M = size(BS_locations,1);
% r = zeros(M_eff,1);
% for ap = 1:M_eff
%     r(ap) = sqrt(sum((UE_location-BS_locations(idx_max(ap),:)).^2)); 
% end
r = zeros(M,1);
% prob = zeros(M,1);
for ap = 1:M
    r(ap) = sqrt(sum((UE_location-BS_locations(ap,:)).^2)); 
end
% [~, bsPriorities] = sort(r,'ascend');
% for ap = 1:M
%     prob(ap) = find(bsPriorities==ap); 
% end
% prob = (1./prob)./sum(1./prob);
% prob = (1./(r.^2))/sum((1./(r.^2)));
prob = (1./r)/sum((1./r));
% prob = (theta.*r)./(theta.*r+omega);
% prob = (1-prob)./sum(1-prob);
r_avg = r'*prob;
% r_avg = mean(r);
% p = (theta*r(1)/(theta*r(1)+omega)) + (theta*r(1)/(theta*r(1)+psi))*(omega/(theta*r(1)+omega));
% p = (theta*r(1)/(theta*r(1)+omega)) + (theta*r(1)/(theta*r(1)+psi))*(theta*r(2)/(theta*r(2)+psi))*(omega/(theta*r(1)+omega));
% p = (theta*r(1)/(theta*r(1)+omega)) + (theta*r(1)/(theta*r(1)+psi))*(theta*r(2)/(theta*r(2)+psi))*(theta*r(3)/(theta*r(3)+psi))*(theta*r(4)/(theta*r(4)+psi))*(omega/(theta*r(1)+omega));
% p = (theta*r(1)/(theta*r(1)+omega)) + (theta*r(1)/(theta*r(1)+psi))^M*(omega/(theta*r(1)+omega));
p = (theta*r_avg/(theta*r_avg+omega)) + (theta*r_avg/(theta*r_avg+psi))^M*(omega/(theta*r_avg+omega));
% p = sum(((theta.*r)./(theta.*r+omega) +  (theta.*r./(theta.*r+psi)).*(omega./(theta.*r+omega))).*prob);
% p = sum(((theta.*r)./(theta.*r+omega)).*prob);
% alpha = theta*(2*R/3);
% p = alpha/(alpha+omega);
p_cov = omega/(theta*r_avg+omega) - (theta*r_avg/(theta*r_avg+psi))^M*(omega/(theta*r_avg+omega));
t = p/(theta*r_avg*p_cov);
end