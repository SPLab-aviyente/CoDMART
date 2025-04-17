function [Ut_indv,k] = updateUt_indv_SetNOC(Lt,wt,Ut_c,MaxNOC,MinNOC,lambda2,lambda4,GAP)
% Definition:
%     This code is used to determine the number of communities in each individual layer at time point t
%     The number of communities is determined using the standard eigen gap criterion.


%   References:
%       [1] E. Al sharoa, M. Alwardat and S. Aviyente. "Community Detection 
%        in Multi-Aspect Functional Brain Networks: Robust Tensor 
%        Decomposition Approach"

%   Author: Esraa Al-sharoa 
%   Address: Jordan University of Science and Technoloogy, EE
%   email: emalsharoa@just.edu.jo

%   Author: Mohammad Al-wardat
%   Address: Michigan State University, ECE
%   email: alwardat@msu.edu
%   Website: GitHub
%   Date: 12-Feb-2024; Last revision: 03-April-2025
%
%   Copyright (c) 2025, Esraa Al-sharoa and Mohammad Al-wardat

%   All rights reserved.



[n,~,M]=size(Lt);
Ut_indv=cell(1,M);
k=zeros(1,M);
Phi_hat=zeros(n,n,M);
for m=1:M
        Dt=1./(diag(sum(Lt(:,:,m)))^(0.5)); Dt(Dt==Inf)=0;
        Phi_hat(:,:,m)=lambda2*Dt*Lt(:,:,m)*Dt+lambda4*wt(m)*(Ut_c*Ut_c');
        Phi_hat(:,:,m)=(Phi_hat(:,:,m)+Phi_hat(:,:,m)')/2;
        
        [Ut_m,Um_EV] = eigs(Phi_hat(:,:,m),MaxNOC,'la');
         v=diag(Um_EV);
        
        EVec = Ut_m;
        vs=v(1:MaxNOC);
        
        %% Determining the number of clusters
        %% Relative EigenGap        

         if strcmpi(GAP,'REG')
        D=vs;

            for kk=3:MaxNOC-1
                    gap(kk)=(D(kk+1)-mean(D(1:kk)))/(1e-6+mean(D(1:kk)));%   
            end

        [~,bb]=max(abs(gap));
        k(m)=bb;
        
 %% Standard EigenGap  
         elseif strcmpi(GAP,'SEG')

            k_max=MaxNOC;
            eigengaps = zeros(length(vs)-1,1);
                    for i=1:1:length(eigengaps)
                        if ((i==k_max))
                            eigengaps(i)=-1;
                        else
                            eigengaps(i)=vs(i+1)-vs(i);
                        end
                    end

           [~,km] = max(abs(eigengaps(MinNOC:end))); 
           k(m)=km+MinNOC-1; 
         end
         
        Ut_indv{1,m}=normc(EVec(:,1:k(m)));
end

end
