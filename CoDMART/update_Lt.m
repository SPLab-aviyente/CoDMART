function [Lt] = update_Lt(Lt,At,St,Xt,Yt,Wt,Ut_indv,lambda2,mu)
% Definition:
%     This code is used to update the low-rank component at time point t

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

%% Check the variable type of Ut_indv
tf = isa(Ut_indv,'cell');


[n,~,M]=size(At);
Psi=zeros(n,n,M);


if tf==1
    
    for m=1:M
            Dt=diag(sum(Lt(:,:,m))); Dt=1./sqrt(Dt); Dt(Dt==Inf)=0;
            Psi(:,:,m)=Dt*(Ut_indv{1,m}*Ut_indv{1,m}')*Dt;
            Lt(:,:,m)=(lambda2*Psi(:,:,m)-mu*St(:,:,m)+mu*At(:,:,m)-Yt(:,:,m)-Xt(:,:,m)+mu*Wt(:,:,m))/(2*mu);
    end    


else
    

    for m=1:M
            Dt=diag(sum(Lt(:,:,m))); Dt=1./sqrt(Dt); Dt(Dt==Inf)=0;
            Psi(:,:,m)=Dt*(Ut_indv(:,:,m)*Ut_indv(:,:,m)')*Dt;
            Lt(:,:,m)=(lambda2*Psi(:,:,m)-mu*St(:,:,m)+mu*At(:,:,m)-Yt(:,:,m)-Xt(:,:,m)+mu*Wt(:,:,m))/(2*mu);    
    end

        

end      
%% make Lt symmetric and nonnegative 
        
        for m=1:M
            L=Lt(:,:,m); 
            L(L<0)=0; 
            L=L-diag(diag(L));
            Lt(:,:,m)=0.5*(L+L');
        end
       
end
