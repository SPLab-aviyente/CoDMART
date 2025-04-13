

function[Xn]=CombineClustersSlowly(k,n_c,c1,c2,m22,s22)
            
%             k=(t-20)*0.2;
            
            x4=truncnormrnd(n_c,m22,s22,0,1);
            
            Xn=[c1  (1-k)*x4+k*c1;(1-k)*x4+k*c2 c2];
            
            
end 
            
           