%---------------------------------------------------------------------------------|
% This is a demo for the CoDMART algorithm. If you find it useful in your         |
% research, please cite the paper below. Please first run the run_me_first.m      |
% --------------------------------------------------------------------------------|

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


clearvars
clc
close all
warning off; 

%% Initializations
T=60; % Number of time points.
NumberOfExp=1; % Number of times to rerun the experiment.
DataSet=cell(1,NumberOfExp); % cell variable that stores the data set at each run. 
FinalVariables=cell(1,NumberOfExp);% cell variable that stores the final
%                                    variables for each run. More details
%                                    in CoDMART function.
VARS=cell(1,NumberOfExp);% cell variable that stores the final
%                                    variables for each run. More details
%                                    in CoDMART function.
nmi_CoDMART=zeros(T,NumberOfExp);% Matrix array that contains the NMI of 
%                                the detected common community structure at each run.
nmi_Layer_CoDMART=zeros(T,NumberOfExp);% Matrix array that contains the
%                                   mean NMI of the detected individual community structure at each run. 
Distance=zeros(T,NumberOfExp);% Matrix array that contains the distance between the consecutive common subspaces at each run.
%                                   
%%

for k=1:NumberOfExp

NumOfLayers=10; % Number of layers
SN=10; % Percent of sparse noise added
MaxNOC=8; % Maximum number of communities
MinNOC=2; % Minimum number of communities

%% Generate the multi-aspect network (MAN)
[TN,GT]=createMyNewSlowlyTemporalNetwork1_n100_c434(0.6,0.2,0.1,0.2,0.6,0.2,0.1,0.2,0.6,0.2,0.1,0.2,SN,SN,SN,NumOfLayers);
% [TN,GT]=createMyNewSlowlyTemporalNetwork1_n100_c456(0.6,0.2,0.1,0.2,0.6,0.2,0.1,0.2,0.6,0.2,0.1,0.2,SN,SN,SN,NumOfLayers);

DataSet{k}=TN;


%% CoDMART
disp(['Processing Experiment ' int2str(k) '  CoDMART ' ])


lambda1=0.1;
lambda2=2;
lambda3=0.5;

 
[FinalVariables{k},VARS{k}] = CoDMART(TN,GT,MaxNOC,MinNOC,lambda1,lambda2,lambda3);


nmi_CoDMART(:,k)=FinalVariables{k}.nmi;
nmi_Layer_CoDMART(:,k)=FinalVariables{k}.MNL;
Distance(:,k)=FinalVariables{k}.Distance;

end


%% Results Summary

mean_nmi_CoDMART=mean(nmi_CoDMART,2);
mean_nmi_Layer_CoDMART=mean(nmi_Layer_CoDMART,2);
mean_Distance=mean(Distance,2);

figure,subplot(2,1,1),plot(2:T,mean_nmi_CoDMART(2:T),'-b','LineWidth',2),grid,title('Average NMI for the detected common community structure','interpreter','latex','FontSize', 14)
axis([0 60 0 1]);xlabel('Time Point','interpreter','latex','FontSize', 14), ylabel('Average NMI','interpreter','latex','FontSize', 14)

subplot(2,1,2),plot(3:T,mean_Distance(3:T),'-r','LineWidth',2),grid,title('Average Distance cost function','interpreter','latex','FontSize', 14)
xlabel('Time Point','interpreter','latex','FontSize', 14), ylabel('Average Distance','interpreter','latex','FontSize', 14)


figure,plot(2:T,mean_nmi_Layer_CoDMART(2:T),'-b','LineWidth',2),grid,title('Average NMI for the detected individual community structure','interpreter','latex','FontSize', 14)
axis([0 60 0 1]);xlabel('Time Point','interpreter','latex','FontSize', 14), ylabel('Average NMI','interpreter','latex','FontSize', 14)


