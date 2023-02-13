% This code was developed by Aref Miri Rekavandi for the paper: Rekavandi, A. M., Seghouane, 
% A. K., & Evans, R. J. (2022). Adaptive Brain Activity Detection in Structured Interference 
% and Partially Homogeneous Locally Correlated Disturbance. IEEE Transactions on Biomedical 
% Engineering, 69(10), 3064-3073. 
% If you use this code in your study, kindly cite the aforementioned paper.
clc
clear all
close all
%% initialization for non-Gaussian case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bw=2;           %Bandwidth
bwhat=bw;
indim=15;       %dimension in time
K=20;           % Number of senondary set
N=10000;        %Number of test samples
p=2;            %Number of columns in H
t=2;            %Number of columns in B
ASNR=10;        %SNR
sigma=3;        %Noise Power
x=1:indim;
H=zeros(indim,p);
freq=[0.1 0.15 0.2];
greq=[-0.1 -0.08 -0.06];
i=sqrt(-1);
for k=1:p
    for j=1:indim
        H(j,k)=1/sqrt(indim)*exp(-1*i*2*pi*freq(k)*(j-1));
    end
end

B=zeros(indim,t);

for k=1:t
    for j=1:indim
        B(j,k)=1/sqrt(indim)*exp(-1*i*2*pi*greq(k)*(j-1));
    end
end

C=[H B];

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

Rband=explicit(noise(:,N+1:N+K),bwhat,K,indim);

SS=0;
for i=1:K
   SS=SS+noise(:,N+i)*noise(:,N+i)';
end
SC=SS/K;
SCmod=SC;
for i=1:indim
    for j=1:indim
        if abs(i-j)>bwhat
            SCmod(i,j)=0;
        end
    end
end

eigcov=eig(cov);
errorSC=norm(SC-cov,'fro');
errorSCmod=norm(SCmod-cov,'fro');
errorRband=norm(Rband-cov,'fro');

pause(5)
[U S V]=svd(SC);
%%%%%%%%%%%%%%%%% Making observations%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu1=[zeros(1,N/2) ones(1,N/2)];

TETA1=random('uniform',0.3,0.3000000001,p,1);
scale=sqrt(((10^(ASNR/10))/((H*TETA1)'*(sigma^2*cov)^(-1)*(H*TETA1))));

for k=1:N   
    X(:,k)=mu1(k)*H*TETA1;
    phi=random('uniform',0.1,0.1000000001,t,1);
    Y(:,k)=scale*X(:,k)+B*phi+sigma*noise(:,k);
end

SNR=10*log10((scale*H*TETA1)'*(sigma^2*cov)^(-1)*(scale*H*TETA1))

%% Diffrent Tests with plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     Htilda=(SC^(-0.5))*H;    
     Btilda=(SC^(-0.5))*B;
     PBO=eye(indim)-((Btilda)*(((Btilda)')*(Btilda))^(-1)*((Btilda)'));
     PHtilda=((Htilda)*(((Htilda)')*(Htilda))^(-1)*((Htilda)'));
     m1=-0.5+sqrt(1+4*(bw+1)*(2*indim-bw))/2;
     m2=-0.5-sqrt(1+4*(bw+1)*(2*indim-bw))/2;
     m=floor(max([0 m1 m2]))+1;
     T=U(:,1:m);
     SC1=T'*SC*T;
     H1=(SC1^(-0.5))*T'*H;    
     B1=(SC1^(-0.5))*T'*B;
     PB1O=eye(m)-((B1)*(((B1)')*(B1))^(-1)*((B1)'));
     G1=PB1O*H1;
     PG1=((G1)*(((G1)')*(G1))^(-1)*((G1)'));

for i=1:N

 temp=(SC^(-0.5))*Y(:,i);
 
 TS1(i)=temp'*PBO*PHtilda*PBO*temp;
 TS2(i)=temp'*PBO*temp;
 ASD(i)=TS1(i)/TS2(i);
 
 temp=Rband^(-0.5)*Y(:,i);
 Bbar=Rband^(-0.5)*B;
 Cbar=Rband^(-0.5)*C;
 PBbarO=eye(indim)-((Bbar)*(((Bbar)')*(Bbar))^(-1)*((Bbar)'));
 PCbarO=eye(indim)-((Cbar)*(((Cbar)')*(Cbar))^(-1)*((Cbar)'));
 TS1(i)=temp'*PBbarO*temp;
 TS2(i)=temp'*PCbarO*temp;
 Proposed1(i)=TS1(i)/TS2(i);
 
 
 temp=Rband^(-0.5)*Y(:,i);
 Bbar=Rband^(-0.5)*B;
 Hbar=Rband^(-0.5)*H;
 PBbarO=eye(indim)-((Bbar)*(((Bbar)')*(Bbar))^(-1)*((Bbar)'));
 PHbar=(Hbar)*(((Hbar)')*(Hbar))^(-1)*((Hbar)');
 TS1(i)=temp'*PBbarO*PHbar*PBbarO*temp;
 TS2(i)=temp'*PBbarO*temp;
 PRao(i)=TS1(i)/TS2(i);
 
 
 temp=Rband^(-0.5)*Y(:,i);
 Bbar=Rband^(-0.5)*B;
 Hbar=Rband^(-0.5)*H;
 Cbar=Rband^(-0.5)*C;
 PBbarO=eye(indim)-((Bbar)*(((Bbar)')*(Bbar))^(-1)*((Bbar)'));
 PCbarO=eye(indim)-((Cbar)*(((Cbar)')*(Cbar))^(-1)*((Cbar)'));
 PHbarBbar=(Hbar)*(((Hbar)')*PBbarO*(Hbar))^(-1)*((Hbar)')*PBbarO;
 TS1(i)=temp'*PHbarBbar'*PHbarBbar*temp;
 TS2(i)=temp'*PCbarO*temp;
 PWald(i)=TS1(i)/TS2(i);

 
 temp=SC^(-0.5)*Y(:,i);
 Bbar=SC^(-0.5)*B;
 Cbar=SC^(-0.5)*C;
 PBbarO=eye(indim)-((Bbar)*(((Bbar)')*(Bbar))^(-1)*((Bbar)'));
 PCbarO=eye(indim)-((Cbar)*(((Cbar)')*(Cbar))^(-1)*((Cbar)'));
 TS1(i)=temp'*PBbarO*temp;
 TS2(i)=temp'*PCbarO*temp;
 GLRT(i)=TS1(i)/TS2(i);
 
 temp=SCmod^(-0.5)*Y(:,i);
 Bbar=SCmod^(-0.5)*B;
 Cbar=SCmod^(-0.5)*C;
 PBbarO=eye(indim)-((Bbar)*(((Bbar)')*(Bbar))^(-1)*((Bbar)'));
 PCbarO=eye(indim)-((Cbar)*(((Cbar)')*(Cbar))^(-1)*((Cbar)'));
 TS1(i)=temp'*PBbarO*temp;
 TS2(i)=temp'*PCbarO*temp;
 AMFSCmod(i)=TS1(i)/TS2(i);
end

figure
subplot(1,5,1)
 scatter(real(ASD(1:N/2)),imag(ASD(1:N/2)))
  title('ASD', 'FontName', 'Times New Roman', ...
        'FontSize',10,'Color','k', 'Interpreter', 'LaTeX')
 hold on
 scatter(real(ASD(N/2+1:N)),imag(ASD(N/2+1:N)))
 subplot(1,5,2)
 scatter(real(GLRT(1:N/2)),imag(GLRT(1:N/2)))
  title('GLRT', 'FontName', 'Times New Roman', ...
        'FontSize',10,'Color','k', 'Interpreter', 'LaTeX')
 hold on
 scatter(real(GLRT(N/2+1:N)),imag(GLRT(N/2+1:N)))
 subplot(1,5,3)
 scatter(real(Proposed1(1:N/2)),imag(Proposed1(1:N/2)))
  title('Banded AMF', 'FontName', 'Times New Roman', ...
        'FontSize',10,'Color','k', 'Interpreter', 'LaTeX')
 hold on
 scatter(real(Proposed1(N/2+1:N)),imag(Proposed1(N/2+1:N)))
subplot(1,5,4)
 scatter(real(PRao(1:N/2)),imag(PRao(1:N/2)))
  title('Banded Rao', 'FontName', 'Times New Roman', ...
        'FontSize',10,'Color','k', 'Interpreter', 'LaTeX')
 hold on
 scatter(real(PRao(N/2+1:N)),imag(PRao(N/2+1:N)))
 subplot(1,5,5)
 scatter(real(PWald(1:N/2)),imag(PWald(1:N/2)))
  title('Banded Wald', 'FontName', 'Times New Roman', ...
        'FontSize',10,'Color','k', 'Interpreter', 'LaTeX')
 hold on
 scatter(real(PWald(N/2+1:N)),imag(PWald(N/2+1:N)))


figure
subplot(1,5,1)
h1 = histogram(real(ASD(1:N/2)),120);
hold on
h2 = histogram(real(ASD(N/2+1:N)),120);
    xlabel(['$\ell(\textbf{y})$'], 'Interpreter', 'LaTeX')
    ylabel('\#Sample', 'Interpreter', 'LaTeX')
    legend({['$\mathcal{H}_0$'],['$\mathcal{H}_1$']}, ...
        'Interpreter', 'LaTeX')
    title('ASD', 'FontName', 'Times New Roman', ...
        'FontSize',10,'Color','k', 'Interpreter', 'LaTeX')  
subplot(1,5,2)
title('')
h1 = histogram(real(GLRT(1:N/2)),120);
hold on
h2 = histogram(real(GLRT(N/2+1:N)),120);
    xlabel(['$\ell(\textbf{y})$'], 'Interpreter', 'LaTeX')
    ylabel('\#Sample', 'Interpreter', 'LaTeX')
    legend({['$\mathcal{H}_0$'],['$\mathcal{H}_1$']}, ...
        'Interpreter', 'LaTeX')
      title('AMF', 'FontName', 'Times New Roman', ...
        'FontSize',10,'Color','k', 'Interpreter', 'LaTeX')
subplot(1,5,3)
title('')
h1 = histogram(real(Proposed1(1:N/2)),120);
hold on
h2 = histogram(real(Proposed1(N/2+1:N)),120);
    xlabel(['$\ell(\textbf{y})$'], 'Interpreter', 'LaTeX')
    ylabel('\#Sample', 'Interpreter', 'LaTeX')
    legend({['$\mathcal{H}_0$'],['$\mathcal{H}_1$']}, ...
        'Interpreter', 'LaTeX')
      title('Banded AMF', 'FontName', 'Times New Roman', ...
        'FontSize',10,'Color','k', 'Interpreter', 'LaTeX')
 subplot(1,5,4)
h1 = histogram(real(PRao(1:N/2)),120);
hold on
h2 = histogram(real(PRao(N/2+1:N)),120);
    xlabel(['$\ell(\textbf{y})$'], 'Interpreter', 'LaTeX')
    ylabel('\#Sample', 'Interpreter', 'LaTeX')
    legend({['$\mathcal{H}_0$'],['$\mathcal{H}_1$']}, ...
        'Interpreter', 'LaTeX')
    title('Banded Rao', 'FontName', 'Times New Roman', ...
        'FontSize',10,'Color','k', 'Interpreter', 'LaTeX') 
    
 subplot(1,5,5)
h1 = histogram(real(PWald(1:N/2)),120);
hold on
h2 = histogram(real(PWald(N/2+1:N)),120);
    xlabel(['$\ell(\textbf{y})$'], 'Interpreter', 'LaTeX')
    ylabel('\#Sample', 'Interpreter', 'LaTeX')
    legend({['$\mathcal{H}_0$'],['$\mathcal{H}_1$']}, ...
        'Interpreter', 'LaTeX')
    title('Banded Wald', 'FontName', 'Times New Roman', ...
        'FontSize',10,'Color','k', 'Interpreter', 'LaTeX')  

ASD1=real(ASD);
GLRT1=real(GLRT);
Proposed11=real(Proposed1);
PRao1=real(PRao);
PWald1=real(PWald);
AMFSCmod1=real(AMFSCmod);
%% ROC calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
last=max([max(ASD1),max(GLRT1),max(Proposed11),max(AMFSCmod1),max(PRao1),max(PWald1)]);
starting=min([min(ASD1),min(GLRT1),min(Proposed11),min(AMFSCmod1),min(PRao1),min(PWald1)]);
th=linspace(starting,last,50000);
for i=1:50000
    pd1(i)=(200*(sum(ASD1(N/2+1:N)>th(i))))/N;
    pf1(i)=(200*(sum(ASD1(1:N/2)>th(i))))/N;
    
    pd2(i)=(200*(sum(GLRT1(N/2+1:N)>th(i))))/N;
    pf2(i)=(200*(sum(GLRT1(1:N/2)>th(i))))/N;  
    
    pd3(i)=(200*(sum(Proposed11(N/2+1:N)>th(i))))/N;
    pf3(i)=(200*(sum(Proposed11(1:N/2)>th(i))))/N; 
    
    pd4(i)=(200*(sum(AMFSCmod1(N/2+1:N)>th(i))))/N;
    pf4(i)=(200*(sum(AMFSCmod1(1:N/2)>th(i))))/N; 
    
    pd5(i)=(200*(sum(PRao1(N/2+1:N)>th(i))))/N;
    pf5(i)=(200*(sum(PRao1(1:N/2)>th(i))))/N; 
    
    pd6(i)=(200*(sum(PWald1(N/2+1:N)>th(i))))/N;
    pf6(i)=(200*(sum(PWald1(1:N/2)>th(i))))/N; 
    i=i+1;
end

%% ROC plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
loglog((pf1),pd1,'b','LineWidth',2)
hold on
loglog((pf2),pd2,'g','LineWidth',2)
hold on
loglog((pf3),pd3,'r','LineWidth',2)
hold on
loglog((pf5),pd5,'k','LineWidth',2)
hold on
loglog(pf6,pd6,'c','LineWidth',2)

 grid on
 xlim([0.02,100])
ylim([4,100]);

    legend({'ASD','GLRT','Banded AMF','Banded Rao','Banded Wald'}, ...
        'Interpreter', 'LaTeX')
    xlabel('Probability of False Alarm (\%)', 'Interpreter', 'LaTeX')
    ylabel('Probability of Detection (\%)', 'Interpreter', 'LaTeX')
    title('', 'FontName', 'Times New Roman', ...
        'FontSize',10,'Color','k', 'Interpreter', 'LaTeX')
