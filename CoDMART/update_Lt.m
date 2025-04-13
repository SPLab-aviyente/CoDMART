function [Lt] = update_Lt(Lt,At,St,Xt,Yt,Wt,Ut_indv,lambda2,mu)

%% Check the variable type of Ut_indv
tf = isa(Ut_indv,'cell');


[n,~,M]=size(At);
Psi=zeros(n,n,M);


if tf==1
    
    for m=1:M
            Dt=diag(sum(Lt(:,:,m))); Dt=1./sqrt(Dt); Dt(Dt==Inf)=0;
            Psi(:,:,m)=Dt*(Ut_indv{1,m}*Ut_indv{1,m}')*Dt;
            Lt(:,:,m)=(lambda2*Psi(:,:,m)-mu*St(:,:,m)+mu*At(:,:,m)-Yt(:,:,m)-Xt(:,:,m)+mu*Wt(:,:,m))/(2*mu);
    end    


else
    

    for m=1:M
            Dt=diag(sum(Lt(:,:,m))); Dt=1./sqrt(Dt); Dt(Dt==Inf)=0;
            Psi(:,:,m)=Dt*(Ut_indv(:,:,m)*Ut_indv(:,:,m)')*Dt;
            Lt(:,:,m)=(lambda2*Psi(:,:,m)-mu*St(:,:,m)+mu*At(:,:,m)-Yt(:,:,m)-Xt(:,:,m)+mu*Wt(:,:,m))/(2*mu);    
    end

        

end      
%% make Lt symmetric and nonnegative 
        
        for m=1:M
            L=Lt(:,:,m); 
            L(L<0)=0; 
            L=L-diag(diag(L));
            Lt(:,:,m)=0.5*(L+L');
        end
       
end