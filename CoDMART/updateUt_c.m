function [Ut_c] = updateUt_c(Lt,Ut_indv,Ut1_c,wt,lambda3,lambda4,numClust)
% Definition:
%     This code is used to update the common subspace across all the layers at time point t

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



tf = isa(Ut_indv,'cell');
[n,~,M]=size(Lt);
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
        [Ut_c,~] = eigs(sumGamma,numClust,'la');
	    EVec = Ut_c;
        
        
        
Ut_c=normc(EVec(:,1:numClust));

end
