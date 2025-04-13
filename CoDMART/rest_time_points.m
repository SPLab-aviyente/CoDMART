function [Ut_c,Ut_indv,wt,VARs]=rest_time_points(At,Ut1_c,lambda1,lambda2,lambda3,lambda4,mu,MaxNOC,MinNOC,Gap)

% Definition:
%     This code returns the variables update for the time points > 1

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



[n,~,M]=size(At);
Lt=zeros(n,n,M); St=Lt; Wt=Lt; 
Yt=Lt; Xt=Lt;
wt=zeros(1,M);
rho=2;
tol=1e-4;

%% Initialize Individual Subspaces Added by me
Ut_indv=zeros(n,MaxNOC,M);
    for jj=1:size(At,3)
        Ut_indv(:,:,jj) = SC(At(:,:,jj),MaxNOC);   
    end

%%



% Ut_indv=zeros(n,numClust,M);

Ut_c=zeros(n,MaxNOC);
max_iter=200;

NOC=MaxNOC;
NOCm=ones(1,M)*MaxNOC;
for i=1:max_iter


        %% Update Wt

        Wt = reshape(wshrinkObj(Lt+Xt/mu,1/mu,[n n M],0,3),[n n M]);

        %% Update St

        [St] = l21_operator(At-Lt+Yt/mu,lambda1,mu);


        %% Update Lt

        Lt = update_Lt(Lt,At,St,Xt,Yt,Wt,Ut_indv,lambda2,mu);

%% Update Ut_c "common subspace"
      
       if i==2
           [Ut_c,NOC] = updateUt_c_SetNOC(Lt,Ut_indv,Ut1_c,wt,lambda3,lambda4,MaxNOC,MinNOC,Gap);
           
       else

           [Ut_c] = updateUt_c(Lt,Ut_indv,Ut1_c,wt,lambda3,lambda4,NOC);
 
       end
%% Update Ut "Indv. subspace"
       if i==2
           
           [Ut_indv,NOCm] = updateUt_indv_SetNOC(Lt,wt,Ut_c,MaxNOC,MinNOC,lambda2,lambda4,Gap);
       else
           Ut_indv = updateUt_indv(Lt,wt,Ut_c,NOCm,lambda2,lambda4);
       end
       
        
%% Update wt

        [wt] = update_wt(Ut_indv,Ut_c,lambda4);

%% Update X & Y 

        dY = At-Lt-St;
        Yt = Yt + mu*dY;

        dX = Lt-Wt;
        Xt = Xt + mu*dX;

%% Check the stopping coditions
        err1=norm(dY(:),'Inf'); err2=norm(dX(:),'Inf');   
        error1(i)=err1; error2(i)=err2; 
        err = max(err1,err2);

        mu = rho*mu;

        if mod(i,10) == 0
        disp(['iter ' num2str(i) ', err=' num2str(err)])
        end
    if err < tol
       break;
    end
end

VARs.Lt=Lt;
VARs.St=St;
VARs.Wt=Wt;
VARs.NOC=NOC;
VARs.NOCm=NOCm;
VARs.wt=wt;
VARs.Error1=error1;
VARs.Error2=error2;

end