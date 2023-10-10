function p = pLoS2(locationsBS, UE_locations, theta, omega, psi)
K = size(UE_locations,1);
M = size(locationsBS,1);
p = zeros(M,K);
for ap = 1:M
    for ue = 1:K
        r = sqrt(sum((UE_locations(ue,:)-locationsBS(ap,:)).^2)); 
        % p (ap,ue) = (theta*r/(theta*r+psi))*(1 + (psi/(theta*r+omega)));
        p (ap,ue) = (theta*r/(theta*r+omega)) + (theta*r/(theta*r+psi))*(omega/(theta*r+omega));
    end
end
end