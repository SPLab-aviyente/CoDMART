function [wt] = update_wt(Ut_indv,Ut_c,lambda4)


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