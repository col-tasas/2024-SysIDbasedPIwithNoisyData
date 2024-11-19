%% Parameters 
A=[-0.53,0.42,-0.44;0.42,-0.56,-0.65;-0.44,-0.65,0.35]; % system matrix A
B=[0.43,-0.82;0.53,-0.78;0.26,0.40]; % input matrix B
Q=[6.12,1.72,0.53;1.72,6.86,1.72;0.53,1.72,5.73]; % weighting matrix Q
R=[1.15,-0.23;-0.23,3.62]; % weighting matrix R
% some parameters about policy iteration 
episode=100; % total number of episodes
noise=L3bounded(0.5,episode);
%K_initial=diag([-1.5,-1,-0.5]); % initial gain K_0
[Pstar,Kstar,~,~]=idare(A,B,Q,R); % K^* and P^* from DARE
xstart=[30;20;10]; % initial state x_0
A_est=A+0.3*A; % initial \hat{A}_0
B_est=B-0.3*B; % initial \hat{B}_0
[~,Ktemp,L,~]=idare(A_est,B_est,Q,R);
K_initial=-Ktemp;
%% Initialization of RLS
nx=size(A,2); % number of states
nu=size(B,2); % number of control inputs
Ahistory3=zeros(nx*nx,episode);
Bhistory3=zeros(nx*nu,episode);
PK3=10*eye(5);  % initial matrix PK=H_t^{-1}
%% iterations of PE and PI(basic parameters)
% Bounded noise with lim approaching 0
xhistory3=zeros(nx,episode);
uhistory3=zeros(nu,episode);
khistory3=zeros(nx*nu,episode);
xhistory3(:,1)=xstart;
khistory3(:,1)=K_initial(:);
K_current3=K_initial;
phistory3=zeros(nx*nx,episode);
factor=10; % additional excitation
%% Start Algorithm
for j=1:episode
    %% Bounded noise with lim approaching 0
    % Policy Evaluation
    Pmb=dlyap((A_est+B_est*K_current3)',Q+K_current3'*R*K_current3); % policy evalutaion based on estiamtes (\hat{A}_{i-1},\hat{B}_{i-1})
    phistory3(:,j)=Pmb(:);
    % RLS
    %uforward=[factor*randn(1,1);factor*randn(1,1);factor*randn(1,1)]; % e_t from Uniform distribution
    uforward=unifrnd(-factor,factor,[2,1]);
    u=uforward+K_current3*xhistory3(:,j); 
    uhistory3(:,j)=u;
    xhistory3(:,j+1)=realsystem(A,B,xhistory3(:,j),u,noise(1,j)); % Collect the data d_t
    [A_est,B_est,PK3]=RLS(A_est,B_est,PK3,xhistory3(:,j),xhistory3(:,j+1),u); % RLS
    Ahistory3(:,j)=A_est(:);
    Bhistory3(:,j)=B_est(:);
    % Policy improvement
    khistory3(:,j)=vec(K_current3);    
    K_current3=-(R+B_est'*Pmb*B_est)\B_est'*Pmb*A_est;  % policy improvement based on estiamtes (\hat{A}_{i},\hat{B}_{i})

    %%
end
save("RLS_PI_EnergyBounded")
subplot(2,1,1);
plot(Normalized_data(phistory3,Pstar),'linewidth',1,'color','green');
xlabel('$\mathrm{iteration~index}$','interpreter','latex','FontSize',12);
ylabel('$\frac{\|\mathbf{P}-\mathbf{\hat{P}}\|}{\|\mathbf{\hat{P}}\|}$','interpreter','latex','FontSize',12);
subplot(2,1,2)
plot(Normalized_data(khistory3,-Kstar),'linewidth',1,'color','green');
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
     mismatch=1;
     for i=1:episode
         noiselevel(i)=energy/i+mismatch;
     end
end

function noiselevel=L3bounded(energy,episode)
     noiselevel=zeros(1,episode);
     for i=1:episode
         noiselevel(i)=energy/(i^2);
     end
end

