function [Ut_c,k] = updateUt_c_SetNOC(Lt,Ut_indv,Ut1_c,wt,lambda3,lambda4,MaxNOC,MinNOC,GAP)
% Definition:
%     This code is used to determine the number of common communities across the layers at time point t
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
tf = isa(Ut_indv,'cell');
sumGamma=zeros(n,n);

        for m=1:M
            if tf==1
                 Gamma_hat=wt(m)*(Ut_indv{1,m}*Ut_indv{1,m}');    
            else

                 Gamma_hat=wt(m)*(Ut_indv(:,:,m)*Ut_indv(:,:,m)'); 
            end
        
          sumGamma=sumGamma+Gamma_hat;
        end
        
        sumGamma=lambda4*sumGamma+lambda3*(Ut1_c*Ut1_c');
        
        sumGamma(sumGamma<0)=0;
        sumGamma=(sumGamma+sumGamma')/2;
        sumGamma=sumGamma-diag(diag(sumGamma));
        
        [Ut_c,Uc_EV] = eigs(sumGamma,MaxNOC,'la');

        %%
        v=diag(Uc_EV);
	    EVec = Ut_c;
        vs=v(1:MaxNOC);
        
        
      %% Determining the number of clusters
      if strcmpi(GAP,'REG')
%% Relative EigenGap  
            D=vs;

            for kk=3:MaxNOC-1
                    gap(kk)=(D(kk+1)-mean(D(1:kk)))/(1e-6+mean(D(1:kk)));%   
            end

        [~,bb]=max(abs(gap));
        k=bb;
  
      elseif strcmpi(GAP,'SEG')
%% Standard EigenGap
        k_max=MaxNOC;
        eigengaps = zeros(length(vs)-1,1);
                for i=1:1:length(eigengaps)
                    if ((i==k_max))
                        eigengaps(i)=-1;
                    else
                        eigengaps(i)=vs(i+1)-vs(i);
                    end
                end
       [~,k] = max(abs(eigengaps(MinNOC:end))); 
       k=k+MinNOC-1;  
      end

Ut_c=normc(EVec(:,1:k));
end
