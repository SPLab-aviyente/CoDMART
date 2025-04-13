function [FinalVariables,VARS] = CoDMART(TN,GT,MaxNOC,MinNOC,lambda1,lambda2,lambda3)
%CoDMART - Main Code for the Community Detection in Multi-Aspect Networks:
%Robust Tensor Approach (CoDMART)
%
%   Inputs:
%       TN - 4D double. It includes the Multi-aspect network of size n x n
%       x T x M
%       GT - 3D double. It contains the ground truth community structure for
%       time points and layers with size n x T x M .
%       MaxNOC - Integer. Maximum number of clusters.
%       MinNOC - Integer. Minimun number of clusters.
%       lambda1,lambda2,lambda3- Regularization parameters.

%
%   Outputs:
%       FinalVariables - a structure that contains:
%          1. UtALL- cell array of size 1 x T: Common subspace martices
%          over time.
%          2. Ut_indv_All- cell array of size 1 x T: individual subspaces
%          of the separate layers.
%          3. wt_All- cell array of size 1 x T: Weights of the separate
%          layers.
%          over time
%          4. Distance- Vector of the distance between the consecutive
%          common.
%          subspaces.
%          5. Plabel- Matrix of the detected common community labels over
%          time, size n x T.
%          6. LayersLabels- cell array of cell arrays that contains the detected individual community 
%          labels across the diffetent layers, size M x T. 
%          7. nmi- Vector array of the normalized mutual information of the
%          detected common structure.
%          8. nmi_Layer- normalized mutual information of the indivdual
%          structure, size M x T.
%          9. MNL- mean nmi across layers.
%          10. D_NOC- Detected common number of communities over time.
%          11. D_NOCm- Detected individual number of communities across layers
%       VARS- cell array of size 1 x T, each cell contains:
%          1. Lt- estimated low-rank approximation for the individual 
%          layers at each time point.
%          2. St- estimated sparse component for the individual 
%          layers at each time point.
%          3. Wt- estimated weights for the individual layers at
%          each time point.
%          4. Error1
%          5. Error2




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


%% Initializations
GT1=GT(:,:,1);
NumOfTime=size(TN,3);
NumOfLayers=size(TN,4);
NumOfNodes=size(TN,1);
D_NOC=zeros(1,NumOfTime);
D_NOCm=zeros(NumOfTime,NumOfLayers);
Plabel=zeros(NumOfNodes,NumOfTime);
nmi=zeros(1,NumOfTime);
LayersLabels=cell(NumOfLayers,NumOfTime);
nmi_Layer=zeros(NumOfLayers,NumOfTime);
A=zeros(NumOfNodes,NumOfNodes,NumOfLayers,NumOfTime);

%%

for pp=1:NumOfLayers
   A(:,:,pp,:)= TN(:,:,:,pp);
    
end
mu=1;
numClust=5; 
%% Choose the eigengap, we are using the SEG
GAP='SEG';%REG for Relative EigenGap . %SEG for Statdard EigenGap


%% Perform the experiment
[Ut_All,Ut_indv_All,wt_All,Distance,VARS] = jer(A,lambda1,lambda2,lambda3,lambda2,mu,numClust,MaxNOC,MinNOC,GAP);

%% Summarize the results

% Detected NOC for the common subspace (D_NOC) and Detected NOC for 
% the indiv. subspaces (D_NOCm)
for kk=2:NumOfTime
    D_NOC(kk)=VARS{1, kk}.NOC;
    D_NOCm(kk,:)=VARS{1, kk}.NOCm;
end


% NMI for the common structure
for v=2:NumOfTime

          % Evaluation metrics 
        Plabel(:,v) = kmeans(Ut_All{v}, D_NOC(v),'replicates',100,'emptyAction','singleton');
        [~, nmi(v), ~] = compute_nmi(GT1(:,v),Plabel(:,v));
end



% NMI for the indiv. structure
for t=2:NumOfTime 
   for l=1:NumOfLayers 
       LayersLabels{l,t}=kmeans(Ut_indv_All{1,t}{1,l},D_NOCm(t,l),'replicates',100,'emptyAction','singleton');
       [~, nmi_Layer(l,t), ~] =compute_nmi(GT1(:,t),LayersLabels{l,t});

   end 
end



mean_nmi_Layer=mean(nmi_Layer);
MNL=mean_nmi_Layer;

FinalVariables.UtALL=Ut_All;
FinalVariables.Ut_indv_All=Ut_indv_All;
FinalVariables.wt_All=wt_All;
FinalVariables.Distance=Distance;
FinalVariables.Plabel=Plabel;
FinalVariables.LayersLabels=LayersLabels;
FinalVariables.nmi=nmi;
FinalVariables.nmi_Layer=nmi_Layer;
FinalVariables.MNL=MNL;
FinalVariables.D_NOC=D_NOC;
FinalVariables.D_NOCm=D_NOCm;
end