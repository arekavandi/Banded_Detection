function Rhat=explicit(X,b,K,N)

Rhat=zeros(b+1,b+1);
SS=0;
for i=1:K
   SS=SS+X(1:b+1,i)*X(1:b+1,i)';
end
Rhat=SS/K;
for i=b+2:N

        F=ones(K,1);
        for m=i-b:i-1
            temp=0;
            for n=1:i-1
                M=Rhat;
                M(n,:)=[];
                M(:,m)=[];
%                 Rhat
                temp=temp+(-1)^(n+m)*(det(M)/det(Rhat(1:i-2,1:i-2)))*X(n,:).';
            end
            F=[F temp];
        end
%         Xhat
      lambda=(F'*F)^(-1)*F'*X(i,:).';
      r1=zeros(i-b-1,1);
      r2=lambda*(det(Rhat)/det(Rhat(1:i-2,1:i-2)));
      r2(1)=[];
      r=[r1;r2];
%       i
%       r
%       Rhat
      corner=(1/K)*conj(X(i,:))*(eye(K)-F*(F'*F)^(-1)*F')*X(i,:).'+r'*Rhat^(-1)*r;
      Rhat=[Rhat r;r' corner];
    end

end