function [M11, GT_all] = createMyNewSlowlyTemporalNetwork1_n100_c434(m1,s1,m2,s2,m11,s11,m22,s22,m1i,s1i,m2i,s2i,PofOnesEdges1,PofOnesEdges2,PofOnesEdges3,num)
%%
%This code creats a temporal network that consists of 100 nodes, 60 time
%points and num layers
%Networks 1-20 consist of 4 communities, Clusters:(30,30,20,20)nodes
%Networks 1-20 parameters: m1,s1,m2,s2
%%%
%Networks 21-40 consist of 3 communities, Clusters:(60,20,20)nodes
%nodes 
%Networks 21-40 parameters: m11,s11,m12,s12
%%%
%Networks 41-60 consist of 4 communities, Clusters:(30,30,20,20)nodes
%Networks 41-60 parameters: m1i,s1i,m2i,s2i

%PofOnesEdges1,PofOnesEdges2,PofOnesEdges3: Percent of edges that will be
%replaced by ones in the off-diagonal of the netwok in Interval 1,
%interval 2 and interval 3, respectively.

%The structure of the network starts changeing slowly at time points 21 
%from 4 clusters to 3 clusters. The structure starts to change slowly 
%again at time point 41 from 3 cluster to 4 clusters

%%%Example:
%[TN,GT]=createMyNewSlowlyTemporalNetwork1_n100_c434(0.8,0.2,0.4,0.3,0.7,0.1,0.3,0.1,0.8,0.2,0.4,0.3,5,5,5,1);

% Contact Information:
% Esraa Al-sharoa
% Jordan University of Science and Technology 
% emalsharoa@just.edu.jo; 
% Date: 09-Feb-2024; Last revision: 03-April-2025
% Copyright (c) 2025, Esraa Al-sharoa 
%   All rights reserved.
%%%%

%%
M11=zeros(100,100,60,num);
GT_all=zeros(100,60,num);




for ii=1:num
TN=zeros(100,100,60);
GT=zeros(100,60);
%% Interval 1
for t=1:20
    
            x1=truncnormrnd(30,m1,s1,0,1);
            x2=truncnormrnd(30,m1,s1,0,1);
            x3=truncnormrnd(20,m1,s1,0,1);
            x5=truncnormrnd(20,m1,s1,0,1);
            
            
            x4=truncnormrnd(100,m2,s2,0,1);
%            
            
            Diago{1}=x1;Diago{2}=x2;Diago{3}=x3;Diago{4}=x5;
            X2=blkdiag(Diago{:});
            R = spones(X2);x4(R==1)=0;
            X1=x4;
            R_L1 = spones(X1);
            L_L1=round((PofOnesEdges1/100)*(nnz(R_L1))); % Calculate the number of edges that will be replaced by ones
            A1_L1=ones(1,nnz(R_L1)); 
            p1_L1 = randperm(length(A1_L1));% permute without repetition
            SS1_L1=find(R_L1==1);%indices of the nonzero values.
            SS1_L1=SS1_L1';
                for kk_L1 = 1 : L_L1
                    % Get a random location in X2.
                   index1_L1 = SS1_L1(p1_L1(kk_L1));
                    % Assign the kkth element of A to it.
                    X1(index1_L1) = A1_L1(kk_L1);
                end
            
            
            
%            
  

X1=X1(1:100,1:100);%offDiagonal
X2=X2(1:100,1:100);%Diagonal

for i=1:100
         for j=i:100
         X1(j,i)=X1(i,j);
         X2(j,i)=X2(i,j);
         end
end

      X=X1+X2;
      X=X- diag(diag(X));
TN(:,:,t)=X;
GT(:,t)=[ones(1,30) 2*ones(1,30) 3*ones(1,20) 4*ones(1,20)];
clear Diago
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interval 2
for t=21:40
    
   
    
            x12=truncnormrnd(20,m11,s11,0,1);
            x22=truncnormrnd(20,m11,s11,0,1);
            
                if t>=21 && t<=25 
                    k=(t-20)*0.2;
                    x25a=truncnormrnd(30,m11,s11,0,1);
                    x25b=truncnormrnd(30,m11,s11,0,1);
                    x25=CombineClustersSlowly(k,30,x25a,x25b,m22,s22);
                else
                x25=truncnormrnd(60,m11,s11,0,1);
                end
            
                x42=truncnormrnd(100,m22,s22,0,1);
                
            Diago2{1}=x25;Diago2{2}=x12;Diago2{3}=x22;
            X22=blkdiag(Diago2{:});
            R2 = spones(X22);x42(R2==1)=0;
            X12=x42;
            
            R_L2 = spones(X12);
            L_L2=round((PofOnesEdges2/100)*(nnz(R_L2))); % Calculate the number of edges that will be replaced by ones
            A1_L2=ones(1,nnz(R_L2)); 
            p1_L2 = randperm(length(A1_L2));% permute without repetition
            SS1_L2=find(R_L2==1);%indices of the nonzero values.
            SS1_L2=SS1_L2';
                for kk_L2 = 1 : L_L2
                    % Get a random location in X2.
                   index1_L2 = SS1_L2(p1_L2(kk_L2));
                    % Assign the kkth element of A to it.
                    X12(index1_L2) = A1_L2(kk_L2);
                end
X12=X12(1:100,1:100);%offDiagonal
X22=X22(1:100,1:100);%Diagonal

for i=1:100
         for j=i:100
         X12(j,i)=X12(i,j);
         X22(j,i)=X22(i,j);
         end
end

      Xn=X12+X22;
      Xn=Xn- diag(diag(Xn));

TN(:,:,t)=Xn;
GT(:,t)=[ones(1,60) 2*ones(1,20) 3*ones(1,20)];
clear Diago2
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Interval 3
for t=41:60
    
            x3i=truncnormrnd(20,m1i,s1i,0,1); 
            x5i=truncnormrnd(20,m1i,s1i,0,1); 
            
            x1ia=truncnormrnd(30,m1i,s1i,0,1);
            x2ib=truncnormrnd(30,m1i,s1i,0,1);
            
            if t>=41 && t<=45 
                    k=(t-40)*0.2;
                    [x1i]=SeparateClustersSlowly(k,30,x1ia,x2ib,m2i,s2i,PofOnesEdges3);
                    Diago3{1}=x1i;Diago3{2}=x3i;Diago3{3}=x5i;
            else
                    Diago3{1}=x1ia;Diago3{2}=x2ib;Diago3{3}=x3i;Diago3{4}=x5i;
            end
            
           

            x4i=truncnormrnd(100,m2i,s2i,0,1);
            
            
%             Diago3{1}=x1i;Diago3{2}=x1i;Diago3{3}=x3i;Diago3{4}=x5i;
            X23=blkdiag(Diago3{:});
            R3 = spones(X23);x4i(R3==1)=0;
            X13=x4i;
            
            

            
            Ri = spones(X13);
            L1=round((PofOnesEdges3/100)*(nnz(Ri))); % Calculate the number of edges that will be replaced by ones
            A1=ones(1,nnz(Ri)); 
            p1 = randperm(length(A1));% permute without repetition
            SS1=find(Ri==1);%indices of the nonzero values.
            SS1=SS1';
                for kk = 1 : L1
                    % Get a random location in X2.
                   index1 = SS1(p1(kk));
                    % Assign the kkth element of A to it.
                    X13(index1) = A1(kk);
                end
            
            
            
            X1i=X13;

X1i=X1i(1:100,1:100);%offDiagonal
X23=X23(1:100,1:100);%Diagonal

for i=1:100
         for j=i:100
         X1i(j,i)=X1i(i,j);
         X23(j,i)=X23(i,j);
         end
end

     Xn3=X1i+X23;
      Xn3=Xn3- diag(diag(Xn3));

TN(:,:,t)=Xn3;
GT(:,t)=[ones(1,30) 2*ones(1,30) 3*ones(1,20) 4*ones(1,20)];
clear Diago3
end

M11(:,:,:,ii)=TN;
GT_all(:,:,ii)=GT;

end

end