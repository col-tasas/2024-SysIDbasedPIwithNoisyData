%% Parameters 
% define the true system and the weighting matrix 
A=[1.01,0.01,0;0.01,1.01,0.01;0,0.01,1.01]; % system matrix A
B=eye(3); % input matrix B
Q=0.001*eye(3); % weighting matrix Q
R=eye(3); % weighting matrix R
% some parameters about policy iteration 
episode=50; % total number of episodes
noise=L2bounded(0.5,episode);
Num=1;% \tau_IPI
%K_initial=diag([-1.5,-1,-0.5]); % initial gain K_0
[Pstar,Kstar,~,~]=idare(A,B,Q,R); % K^* and P^* from DARE
xstart=[30;20;10]; % initial state x_0
%A_est=0*eye(3); % initial \hat{A}_0
A_est=A+0.5*eye(3);
%B_est=0*eye(3); % initial \hat{B}_0
B_est=B+0.5*eye(3); % initial \hat{B}_0
[~,Ktemp,L,~]=idare(A_est,B_est,Q,R);
K_initial=-Ktemp;
%% Initialization of RLS
nx=size(A,2); % number of states
nu=size(B,2); % number of control inputs
Ahistory2=zeros(nx*nx,episode);
Bhistory2=zeros(nx*nu,episode);
PK2=10*eye(6);  % initial matrix PK=H_t^{-1}
%% iterations of PE and PI(basic parameters)
% Bounded noise with lim approaching 0
xhistory2=zeros(nx,episode);
uhistory2=zeros(nu,episode);
khistory2=zeros(nx*nu,episode);
xhistory2(:,1)=xstart;
khistory2(:,1)=K_initial(:);
K_current2=K_initial;
phistory2=zeros(nx*nx,episode);
factor=10; % additional excitation
%% Start Algorithm
for j=1:episode
    %% Bounded noise with lim approaching 0
    % Policy Evaluation
    Pmb=dlyap((A_est+B_est*K_current2)',Q+K_current2'*R*K_current2); % policy evalutaion based on estiamtes (\hat{A}_{i-1},\hat{B}_{i-1})
    phistory2(:,j)=Pmb(:);
    % RLS
    %uforward=[factor*randn(1,1);factor*randn(1,1);factor*randn(1,1)]; % e_t from Uniform distribution
    uforward=unifrnd(-factor,factor,[3,1]);
    u=uforward+K_current2*xhistory2(:,j); 
    uhistory2(:,j)=u;
    xhistory2(:,j+1)=realsystem(A,B,xhistory2(:,j),u,noise(1,j)); % Collect the data d_t
    [A_est,B_est,PK2]=RLS(A_est,B_est,PK2,xhistory2(:,j),xhistory2(:,j+1),u); % RLS
    Ahistory2(:,j)=A_est(:);
    Bhistory2(:,j)=B_est(:);
    % Policy improvement
    khistory2(:,j)=vec(K_current2);    
    K_current2=-(R+B_est'*Pmb*B_est)\B_est'*Pmb*A_est;  % policy improvement based on estiamtes (\hat{A}_{i},\hat{B}_{i})

    %%
end
save("BoundedNoise");
subplot(2,1,1);
plot(Normalized_data(phistory2,Pstar),'linewidth',1,'color','green');
xlabel('$\mathrm{iteration~index}$','interpreter','latex','FontSize',12);
ylabel('$\frac{\|\mathbf{P}-\mathbf{\hat{P}}\|}{\|\mathbf{\hat{P}}\|}$','interpreter','latex','FontSize',12);
subplot(2,1,2)
plot(Normalized_data(khistory2,-Kstar),'linewidth',1,'color','green');
xlabel('$\mathrm{iteration~index}$','interpreter','latex','FontSize',12);
ylabel('$\frac{\|\mathbf{K}-\mathbf{\hat{K}}\|}{\|\mathbf{\hat{K}}\|}$','interpreter','latex','FontSize',12);
clearvars
%% required function
% real system dynamic
function xnew=realsystem(A,B,state,controlinput,noise)
    noise=noisegeneration(noise);
    xnew=A*state+B*controlinput+noise;
end
function [A,B,P_new]=RLS(A_old,B_old,P_old,x_old,x_new,u)
%% Inputs:
% A_old B_old: previous system dynamic
% P_old: H_{t-1}^{-1}
% x_old x_new u: x_t,x_{t+1},u_t
%% Outputs:
% A B: updated estimated system dynamic
% P_new: H_t^{-1}
    d=[x_old;u];
    P_new=P_old-P_old*(d)*d'*P_old/(1+d'*P_old*d); % rank 1 update: as stated in Remark. 2
    EST=[A_old,B_old];
    EST_new=EST+(x_new-EST*d)*d'*P_new;
    [a,~]=size(A_old);
    A=EST_new(:,1:a);
    B=EST_new(:,a+1:end);
end
function output=vec(input)
% function vec: vectorize a matrix
    output=input(:);
end
function figuredata = Normalized_data(data1,ref)
    data2=kron(ones(1,size(data1,2)),vec(ref));
    figuredata=(sqrt(sum((data1-data2).*(data1-data2),1)))/sqrt(sum(ref(:).*ref(:)));
end

function noise=noisegeneration(energy)
    noise=randn(3,1);
    current_norm=norm(noise,'fro');
    noise=noise*(energy/current_norm);
end

function noiselevel=L1bounded(energy,episode)
     noiselevel=zeros(1,episode);
     for i=1:episode
         noiselevel(i)=energy/i;
     end
end

function noiselevel=L2bounded(energy,episode)
     noiselevel=zeros(1,episode);
     mismatch=0.5;
     for i=1:episode
         noiselevel(i)=energy/i+mismatch;
     end
end

