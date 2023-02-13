% This code was developed by Aref Miri Rekavandi for the paper: Rekavandi, A. M., Seghouane, 
% A. K., & Evans, R. J. (2022). Adaptive Brain Activity Detection in Structured Interference 
% and Partially Homogeneous Locally Correlated Disturbance. IEEE Transactions on Biomedical 
% Engineering, 69(10), 3064-3073. 
% If you use this code in your study, kindly cite the aforementioned paper.
clc
clear all
close all
%% initialization for non-Gaussian case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K=10;
for it=1:10
K=K*2;
index(it)=K;
        it
for run=1:30
bw=2;
indim=5;
N=100;
p=2;
t=2;
ASNR=10;
sigma=3;

cov1=zeros(indim,indim);
for i=1:indim
    for j=1:indim
         cov1(i,j)=rand();
    end
end
cov=cov1*cov1';
for i=1:indim
    for j=1:indim
        if abs(i-j)>bw
            cov(i,j)=0;
        end
    end
end
for i=1:indim
    cov(i,i)=sum(cov(i,:))+3*rand();
end

a=random('normal',0,1,indim,N+K);
b=random('normal',0,1,indim,N+K);
noise=(cov/(2))^(0.5)*(a+sqrt(-1)*b);

%%%%%%%%%%% Estimating covariance matrix from secondary data%%%%%%

Rband=explicit(noise(:,N+1:N+K),bw,K,indim);

SS=0;
for i=1:K
   SS=SS+noise(:,N+i)*noise(:,N+i)';
end
SC=SS/K;

SCmod=SC;
for i=1:indim
    for j=1:indim
        if abs(i-j)>bw
            SCmod(i,j)=0;
        end
    end
end

errorSC(run)=norm(SC-cov,'fro');

errorRband(run)=norm(Rband-cov,'fro');
    end
    er1(it)=mean(errorSC);

    er3(it)=mean(errorRband);
end
figure
plot(log(index),er1,'--^','LineWidth',2)
hold on

plot(log(index),er3,'--*','LineWidth',2)
xlim([3 8.6])
ylim([0 7.1])
grid on
    ylabel('Frobenius Norm of Error', 'Interpreter', 'LaTeX')
    xlabel('log(K)', 'Interpreter', 'LaTeX')
legend({'Sample Covariance','Banded Covariance'}, ...
        'Interpreter', 'LaTeX')

