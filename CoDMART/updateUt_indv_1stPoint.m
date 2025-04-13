function [Ut_indv] = updateUt_indv_1stPoint(Lt,wt,Ut_c,numClust,lambda2,lambda4)



[n,~,M]=size(Lt);
Phi_hat=zeros(n,n,M);
Ut_indv=zeros(n,numClust,M);

    for m=1:M
            Dt=1./(diag(sum(Lt(:,:,m)))^(0.5)); Dt(Dt==Inf)=0;
            Phi_hat(:,:,m)=lambda2*Dt*Lt(:,:,m)*Dt+lambda4*wt(m)*(Ut_c*Ut_c');
            Phi_hat(:,:,m)=(Phi_hat(:,:,m)+Phi_hat(:,:,m)')/2;
            [Ut_indv(:,:,m),~] = eigs(Phi_hat(:,:,m),numClust,'la');
            Ut_indv(:,:,m)=normc(Ut_indv(:,:,m));
    end



end