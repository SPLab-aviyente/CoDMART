function [Ut_indv] = updateUt_indv_1stPoint(Lt,wt,Ut_c,numClust,lambda2,lambda4)
% Definition:
%     This code is used to update the individual subspaces of the different layers at the first time point 

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
Phi_hat=zeros(n,n,M);
Ut_indv=zeros(n,numClust,M);

    for m=1:M
            Dt=1./(diag(sum(Lt(:,:,m)))^(0.5)); Dt(Dt==Inf)=0;
            Phi_hat(:,:,m)=lambda2*Dt*Lt(:,:,m)*Dt+lambda4*wt(m)*(Ut_c*Ut_c');
            Phi_hat(:,:,m)=(Phi_hat(:,:,m)+Phi_hat(:,:,m)')/2;
            [Ut_indv(:,:,m),~] = eigs(Phi_hat(:,:,m),numClust,'la');
            Ut_indv(:,:,m)=normc(Ut_indv(:,:,m));
    end



end
