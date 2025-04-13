function [Ut_All,Ut_indv_All,wt_All,Distance,VARS] = jer(A,lambda1,lambda2,lambda3,lambda4,mu,numClust,MaxNOC,MinNOC,Gap)
% 
% Definition:
%     This code returns the varibles update for the multi-aspect network
%     over time.
%
% Inputs:
%       A - 4D double. It includes the Multi-aspect network of size n x n
%       x M x T
%       lambda1,lambda2,lambda3- Regularization parameters where
%       lambda4=lambda2
%       numClust- Integer. Initial number of communities at t=1 
%       MaxNOC - Integer. Maximum number of clusters.
%       MinNOC - Integer. Minimun number of clusters.
%       Gap- string. Eigengap used to determine the number of communities
%         (Gap='SEG')
% 
% 
% Outputs:
%          UtALL- cell array of size 1 x T: Common subspace martices
%          over time.
%          Ut_indv_All- cell array of size 1 x T: individual subspaces
%          of the separate layers.
%          wt_All- cell array of size 1 x T: Weights of the separate
%          layers.
%          over time
%          Distance- Vector of the distance between the consecutive
%          common.
%          subspaces.
%       VARS- cell array of size 1 x T, each cell contains:
%          1. Lt- estimated low-rank approximation for the individual 
%          layers at each time point.
%          2. St- estimated sparse component for the individual 
%          layers at each time point.
%          3. Wt- estimated weights for the individual layers at
%          each time point.
%          4. Error1
%          5. Error2

% 
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
%
%
% required external files:
%  [1] wshrinkObj, "On unifying multi-view self-representations for clustering by tensor multi-rank minimization."
%  [2] Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

  
%%%%


%% Initailization 

[n,~,~,tp]=size(A);
Ut_c=zeros(n,numClust(1));
Ut_All = {[0;0]};                          %% to store all Ut_c
Distance=zeros(1,tp);
VARS=cell(1,tp);

%%

Ut_indv_All = {[0;0]};                     %% to store all indiviual subspaces
wt_All = {[0;0]};


for t=1:1:tp  % for-loop for sliding window
    
    
    
    At = A(:,:,:,t);   %% The multiplex tensor at time t

    Ut1_c=Ut_c;        %% store the prevoius common subspace 


    if t==1
        fprintf('-------- Start, time point number %d--------\n', t);
        [Ut_c,Ut_indv,wt,VARs]=first_time_point(At,lambda1,lambda2,lambda4,mu,numClust(t));
        fprintf('-------- End, time point number %d--------\n', t);
    else
        fprintf('-------- Start, time point number %d--------\n', t);tic
        [Ut_c,Ut_indv,wt,VARs]=rest_time_points(At,Ut1_c,lambda1,lambda2,lambda3,lambda4,mu,MaxNOC,MinNOC,Gap);
        fprintf('-------- End, time point number %d--------\n', t);toc
    end

%% Store All the Results

        Ut_All{t} = Ut_c;           % store the common subspace
        Ut_indv_All{t}=Ut_indv;     % store the indv. subspaces
        wt_All{t}=wt;               % store the weights vector
        VARS{t}=VARs;
        
        
            if t>=2
                UcUcT=(Ut_c*Ut_c');
                Distance(t)=trace(UcUcT*(Ut1_c*Ut1_c')); % Calculate the Distance 
            end
        
        
end
        
end


