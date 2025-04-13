function [Ut_indv] = updateUt_indv(Lt,wt,Ut_c,numClust,lambda2,lambda4)



[n,~,M]=size(Lt);
Ut_indv=cell(1,M);
Phi_hat=zeros(n,n,M);
for m=1:M
        Dt=1./(diag(sum(Lt(:,:,m)))^(0.5)); Dt(Dt==Inf)=0;
        Phi_hat(:,:,m)=lambda2*Dt*Lt(:,:,m)*Dt+lambda4*wt(m)*(Ut_c*Ut_c');
        Phi_hat(:,:,m)=(Phi_hat(:,:,m)+Phi_hat(:,:,m)')/2;

        [Ut_indv11,~] = eigs(Phi_hat(:,:,m),numClust(m),'la');
        EVec = Ut_indv11;
       
        Ut_indv{1,m}=EVec(:,1:numClust(m));
      
       
        Ut_indv{1,m}=normc(Ut_indv{1,m});
end

end