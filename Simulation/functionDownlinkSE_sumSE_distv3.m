function SE = functionDownlinkSE_sumSE_distv3(beta,preLogFactor,L,K,K_mmW,N_AP,N_UE,D,rhomax,tau_p)
%Compute downlink SE according to Corollary 6.3 with sum SE maximizing power allocation
%in Algorithm 7.6
%
%INPUT:
%bk             = Matrix with dimension L x K where (1:La,k) is the non-zero
%                 portion of \tilde{b}_k in (7.25) (La is the number of
%                 APs serving UE k)
%Ck             = Matrix with dimension L x L x K x K where (1:La,1:La,k,i) is the 
%                 non-zero portion of \tilde{C}_{ki} in (7.26) (La is the
%                 number of APs serving UE i)
%preLogFactor   = Pre-log factor in the SE expression of Corollary 6.3
%L              = Number of APs
%K              = Number of UEs in the network
%rhomax         = Maximum allowable AP transmit power
%tau_p          = Length of pilot sequences
%
%OUTPUT:
%SE             = SEs achieved with sum SE maximizing power allocation algorithm
%
%
%This Matlab function was developed to generate simulation results to:
%
%Ozlem Tugfe Demir, Emil Bjornson and Luca Sanguinetti (2021),
%"Foundations of User-Centric Cell-Free Massive MIMO", 
%Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
%pp 162-472. DOI: 10.1561/2000000109
%
%This is version 1.0 (Last edited: 2021-01-31)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.

%Initialize the current objective function
objec_new = inf;

%Prepare array to store the number of APs serving a specficic UE
La = zeros(K,1);
%Prepare cell to store the AP indices serving a specficic UE
Serv = cell(K,1);
%Prepare cell to store the AP indices not serving a specficic UE
NoServ = cell(K,1);
p_fac=1000; %ratio of mmW to sub-6 powers
%Construc the above array and cells
for k = 1:K
    servingAPs = find(D(:,k)==1);
    NoservingAPs = find(D(:,k)==0);
    
    Serv{k} = servingAPs;
    NoServ{k} = NoservingAPs;
    
    La(k) = length(servingAPs);
    beta(:,k) = beta(:,k).*D(:,k);
end
beta_opt = zeros(sum(La),K);
for k=1:K
    if (k==1)
        beta_opt(1:La(k),k) = beta(Serv{k},k);
    else
        beta_opt(1+sum(La(1:k-1)):sum(La(1:k)),k) = beta(Serv{k},k);
    end
end
%Intialize the difference between current and previous objective values 
diff = 100;
%Initialize iterates
eta_eq = zeros(L,K);
if (K_mmW == 0)
    for l = 1:L
        for k = 1:K
            if ismember(l,Serv{k})
                eta_eq(l,k) = 1./(N_AP*N_UE*sum(beta(l,:)));
            end
        end
    end
else
    for l = 1:L
        for k = 1:K
            if ismember(l,Serv{k})
                if (k<=K_mmW)
                    eta_eq(l,k) = p_fac./(N_AP*N_UE*(p_fac*beta(l,1)+sum(beta(l,2:K))));
                else
                    eta_eq(l,k) = 1./(N_AP*N_UE*(p_fac*beta(l,1)+sum(beta(l,2:K))));
                end
            end
        end
    end
end
lambda_eq = zeros(K,1); %sum((sqrt(eta_eq)*D).*beta,1)';
zeta_eq = zeros(K,1);
for k = 1:K
    lambda_eq(k) = (sqrt(eta_eq(:,k)))'*beta(:,k);
    zeta_eq(k) = (lambda_eq(k)^2)/(1/(rhomax*N_AP*N_AP) + (N_UE/N_AP)*beta(:,k)'*sum(beta.*eta_eq,2));
end
lambda_old = lambda_eq;
zeta_old = zeta_eq;
% eta = eta_eq;
SE_eq = preLogFactor*log(1+zeta_eq)/log(2);
SE = zeros(K,1);
%Initizalize the iteration counter to zero
iterr = 0;
n_sca = 2;
%Go through the algorithm steps if the objective function is improved
%more than 0.1 or not improved at all
while (diff>0.1) || (diff<0) || (iterr > n_sca)
    %Increase iteration counter by one
    iterr = iterr+1;
    %Update the previous objective value by the current objective value
    objec_old = objec_new;
    if (K_mmW == 0)
        %Solve the convex problem in (7.33) with CVX
        cvx_begin quiet
        variable t(K) 
        variable zeta(K)
        variable lambda(K)
        variable c2(sum(La),1)
    %     variable c(L,K)
        maximize sum(t)
        subject to
        
        for k=1:K
            if(La(K) > 0)
                t(k) - preLogFactor*log(1+zeta(k))/log(2)<=0;
        %         ((N_UE/N_AP)*beta(:,k)'*sum(beta.*(c.^2),2) + (1/(rhomax*N_AP*N_AP)))*zeta_old(k)^2 <= 2*lambda_old(k)*zeta_old(k)*(lambda(k)-lambda_old(k)) - lambda_old(k)*(zeta(k)-zeta_old(k))
        %         ((N_UE/N_AP)*beta_opt(:,k)'*sum(beta_opt.*(c2.^2),2) + (1/(rhomax*N_AP*N_AP)))*zeta_old(k)^2 <= 2*lambda_old(k)*zeta_old(k)*(lambda(k)-lambda_old(k)) - lambda_old(k)*(zeta(k)-zeta_old(k));       
                sum1 = cvx_zeros([1,1]);
                sum1 = sum1 + (1/(rhomax*N_AP*N_AP))*zeta_old(k)^2;
                for i = 1:K
                    if(La(i) > 0)
                        sum1 = sum1 + zeta_old(k)^2*((N_UE/N_AP)*beta(Serv{i},k)'*(beta(Serv{i},i).*(c2(1+sum(La(1:i-1)):sum(La(1:i))).^2)));
                    end
                end
                sum1 <= 2*lambda_old(k)*zeta_old(k)*(lambda(k)-lambda_old(k)) - lambda_old(k)*(zeta(k)-zeta_old(k));
                lambda(k) <= c2(1+sum(La(1:k-1)):sum(La(1:k)))'*beta_opt(1+sum(La(1:k-1)):sum(La(1:k)),k);
            end
        end
        for l = 1:L
            sum2 = cvx_zeros([1,1]);
    %     for l = 1:sum(La)
            for k = 1:K
    %             beta(l,:)*(c(l,:).^2)'<= 1/(N_AP*N_UE);            
    %         norm(sqrt(beta(l,:))*(c(l,:))')<= 1/sqrt(N_AP*N_UE); 
                [a,b] = ismember(l,Serv{k});
                if a
                    sum2 = sum2 + beta(l,k)*c2(sum(La(1:k-1))+b)^2;
                end
            end
            sum2 <= 1/(N_AP*N_UE);            
        end
        t >= zeros(K,1);
    %     c >= zeros(L,K);
        c2 >= zeros(sum(La),1);
        lambda>=zeros(K,1);
    %     for k=1:K 
    % %         c(Serv{k},k) == c2(1+sum(La(1:k-1)):sum(La(1:k)),k);
    % %         c(NoServ{k},k) == zeros(length(NoServ{k}),1);
    %         c2(1:sum(La(1:k-1)),k) == zeros(sum(La(1:k-1)),1);
    %         c2(1+sum(La(1:k)):sum(La),k) == zeros(sum(La)-sum(La(1:k)),1);
    %     end
        cvx_end
    else
       %Solve the convex problem in (7.33) with CVX
        cvx_begin quiet
        variable t(K) 
        variable zeta(K)
        variable lambda(K)
        variable c2(sum(La),1)
    %     variable c(L,K)
        maximize sum(t)
        subject to
        
        for k=1:K
            if(La(K) > 0)
                t(k) - preLogFactor*log(1+zeta(k))/log(2)<=0;
        %         ((N_UE/N_AP)*beta(:,k)'*sum(beta.*(c.^2),2) + (1/(rhomax*N_AP*N_AP)))*zeta_old(k)^2 <= 2*lambda_old(k)*zeta_old(k)*(lambda(k)-lambda_old(k)) - lambda_old(k)*(zeta(k)-zeta_old(k))
        %         ((N_UE/N_AP)*beta_opt(:,k)'*sum(beta_opt.*(c2.^2),2) + (1/(rhomax*N_AP*N_AP)))*zeta_old(k)^2 <= 2*lambda_old(k)*zeta_old(k)*(lambda(k)-lambda_old(k)) - lambda_old(k)*(zeta(k)-zeta_old(k));       
                sum1 = cvx_zeros([1,1]);
                sum1 = sum1 + (1/(rhomax*N_AP*N_AP))*zeta_old(k)^2;
                for i = 1:K
                    if(La(i) > 0)
                        sum1 = sum1 + zeta_old(k)^2*((N_UE/N_AP)*beta(Serv{i},k)'*(beta(Serv{i},i).*(c2(1+sum(La(1:i-1)):sum(La(1:i))).^2)));
                    end
                end
                sum1 <= 2*lambda_old(k)*zeta_old(k)*(lambda(k)-lambda_old(k)) - lambda_old(k)*(zeta(k)-zeta_old(k));
                lambda(k) <= c2(1+sum(La(1:k-1)):sum(La(1:k)))'*beta_opt(1+sum(La(1:k-1)):sum(La(1:k)),k);
            end
        end
        for l = 1:L
            sum2 = cvx_zeros([1,1]);
    %     for l = 1:sum(La)
            for k = 1:K
    %             beta(l,:)*(c(l,:).^2)'<= 1/(N_AP*N_UE);            
    %         norm(sqrt(beta(l,:))*(c(l,:))')<= 1/sqrt(N_AP*N_UE); 
                [a,b] = ismember(l,Serv{k});
                if a
                    sum2 = sum2 + beta(l,k)*c2(sum(La(1:k-1))+b)^2;
                end
            end
            sum2 <= 1/(N_AP*N_UE);            
        end
        t >= zeros(K,1);
    %     c >= zeros(L,K);
        c2 >= zeros(sum(La),1);
        lambda>=zeros(K,1);
    %     for k=1:K 
    % %         c(Serv{k},k) == c2(1+sum(La(1:k-1)):sum(La(1:k)),k);
    % %         c(NoServ{k},k) == zeros(length(NoServ{k}),1);
    %         c2(1:sum(La(1:k-1)),k) == zeros(sum(La(1:k-1)),1);
    %         c2(1+sum(La(1:k)):sum(La),k) == zeros(sum(La)-sum(La(1:k)),1);
    %     end
        cvx_end
    end
    try
        if (cvx_status == 'Solved')
            %Update the power allocation coefficients 
            %obtained by CVX
            eta = zeros(L,K);
            for k=1:K
                eta(Serv{k},k) = c2(1+sum(La(1:k-1)):sum(La(1:k))).^2;
            end
            %Update the current objective value
            lambda_old = sum(sqrt(eta).*beta,1)';
            zeta_old = zeros(K,1);
            for k = 1:K
                zeta_old(k) = (lambda_old(k)^2)/(1/(rhomax*N_AP*N_AP) + (N_UE/N_AP)*beta(:,k)'*sum(beta.*eta,2));
            end
        %     zeta_old = zeta;
        %     lambda_old = lambda;
            %Compute SEs
            SE = preLogFactor*log(1+zeta_old)/log(2);
            objec_new = sum(SE);
        
            %Obtain the difference between current and previous objective values
            diff = objec_new - objec_old;
            clear c2 t zeta lambda
        else
            break;
        end
    catch 
        if (cvx_status == 'Inaccurate/Solved')
            %Update the power allocation coefficients 
            %obtained by CVX
            eta = zeros(L,K);
            for k=1:K
                eta(Serv{k},k) = c2(1+sum(La(1:k-1)):sum(La(1:k))).^2;
            end
            %Update the current objective value
            lambda_old = sum(sqrt(eta).*beta,1)';
            zeta_old = zeros(K,1);
            for k = 1:K
                zeta_old(k) = (lambda_old(k)^2)/(1/(rhomax*N_AP*N_AP) + (N_UE/N_AP)*beta(:,k)'*sum(beta.*eta,2));
            end
        %     zeta_old = zeta;
        %     lambda_old = lambda;
            %Compute SEs
            SE = preLogFactor*log(1+zeta_old)/log(2);
            objec_new = sum(SE);
        
            %Obtain the difference between current and previous objective values
            diff = objec_new - objec_old;
            clear c2 t zeta lambda
        else
            break;
        end
    end
end
if (sum(SE) < sum(SE_eq))
    SE=SE_eq;
% SE = objec_new;
end
end