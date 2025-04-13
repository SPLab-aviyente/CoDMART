function [Ut_c] = updateUt_c(Lt,Ut_indv,Ut1_c,wt,lambda3,lambda4,numClust)

tf = isa(Ut_indv,'cell');
[n,~,M]=size(Lt);
sumGamma=zeros(n,n);

        for m=1:M
            
             if tf==1
                 Gamma_hat=wt(m)*(Ut_indv{1,m}*Ut_indv{1,m}');    
             else
              Gamma_hat=wt(m)*(Ut_indv(:,:,m)*Ut_indv(:,:,m)');        

             end
             
          sumGamma=sumGamma+Gamma_hat;
        end
        
        sumGamma=lambda4*sumGamma+lambda3*(Ut1_c*Ut1_c');
        
        sumGamma(sumGamma<0)=0;
        sumGamma=(sumGamma+sumGamma')/2;
        [Ut_c,~] = eigs(sumGamma,numClust,'la');
	    EVec = Ut_c;
        
        
        
Ut_c=normc(EVec(:,1:numClust));

end