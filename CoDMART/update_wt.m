function [wt] = update_wt(Ut_indv,Ut_c,lambda4)
% Definition:
%     This code is used to update the weights vector optimized across the different layers at time point t

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


if tf==1
    M=size(Ut_indv,2);
    wtt=zeros(1,M);
    wt=zeros(1,M);

    for m=1:M
        wtt(m) = lambda4*trace((Ut_indv{1,m}*Ut_indv{1,m}')*(Ut_c*Ut_c'));
       
    end
        
    for m=1:M
          wt(m) = wtt(m)/norm(wtt);
    end
    
    
    
else
    [~,~,M]=size(Ut_indv);
    wtt=zeros(1,M);
    wt=zeros(1,M);

      for m=1:M
        wtt(m) = lambda4*trace((Ut_indv(:,:,m)*Ut_indv(:,:,m)')*(Ut_c*Ut_c'));
        
      end
      for m=1:M
          wt(m) = wtt(m)/norm(wtt);
      end 
     
end
end
