

function[Xn]=SeparateClustersSlowly(k,n_c,c1,c2,m22,s22,PofOnesEdges)
            
%             k=(t-20)*0.2;
            
            x4=truncnormrnd(n_c,m22,s22,0,1);
            
            
            Ri = spones(x4);
            L1=round((PofOnesEdges/100)*(nnz(Ri))); % Calculate the number of edges that will be replaced by ones
            A1=ones(1,nnz(Ri)); 
            p1 = randperm(length(A1));% permute without repetition
            SS1=find(Ri==1);%indices of the nonzero values.
            SS1=SS1';
                for kk = 1 : L1
                    % Get a random location in X2.
                   index1 = SS1(p1(kk));
                    % Assign the kkth element of A to it.
                    x4(index1) = A1(kk);
                end
            
            
            Xn=[c1  (1-k)*c1+k*x4;(1-k)*c2+k*x4 c2];
            
            
end 
            
           