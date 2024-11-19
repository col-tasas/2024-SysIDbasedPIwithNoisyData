load("BoundedNoise.mat");
load("BoundedNoiseApproach0.mat");
load("EnergyBoundedNoise.mat");

% A=[-0.53,0.42,-0.44;0.42,-0.56,-0.65;-0.44,-0.65,0.35]; % system matrix A
% B=[0.43,-0.82;0.53,-0.78;0.26,0.40]; % input matrix B
% Q=[6.12,1.72,0.53;1.72,6.86,1.72;0.53,1.72,5.73]; % weighting matrix Q
% R=[1.15,-0.23;-0.23,3.62]; % weighting matrix R

[Pstar,Kstar,L,info]=idare(A,B,Q,R);

Pstarmatrix=kron(ones(1,episode),vec(Pstar));
Kstarmatrix=kron(ones(1,episode),vec(-Kstar));


% PG
P1=timeseries(Normfigure2(phistory1,Pstarmatrix,Pstar),0:1:episode-1);
K1=timeseries(Normfigure2(khistory1,Kstarmatrix,-Kstar),0:1:episode-1);
% IPI 
P2=timeseries(Normfigure2(phistory2,Pstarmatrix,Pstar),0:1:episode-1);
K2=timeseries(Normfigure2(khistory2,Kstarmatrix,-Kstar),0:1:episode-1);
% DPI
P3=timeseries(Normfigure2(phistory3,Pstarmatrix,Pstar),0:1:episode-1);
K3=timeseries(Normfigure2(khistory3,Kstarmatrix,-Kstar),0:1:episode-1);

figure;
subplot(2,1,1);

plot(P2,'linewidth',1.5,'color','blue');
hold on
plot(P1,'linewidth',1.5,'color','red','LineStyle','--');
hold on
plot(P3,'linewidth',1.5,'color','m','LineStyle',':');
hold on


legend("$\mathrm{PB}1~(42a)$","$\mathrm{PB}2~(42b)$","$\mathrm{EB}~~(42c)$",'interpreter','latex');
xlim([0,15]);
ylim([0,12])
title('')
xlabel('$\mathrm{Timesteps}~{t}$','interpreter','latex','FontSize',12);
ylabel('${\frac{|{\hat{P}_t}-{P^*}|}{|{P^*}|}}$','interpreter','latex','FontSize',12);
%title("Convergence of Cost Function Kernel P")
subplot(2,1,2)

plot(K2,'linewidth',1.5,'color','blue');
hold on
plot(K1,'linewidth',1.5,'color','red','LineStyle','--');
hold on
plot(K3,'linewidth',1.5,'color','m','LineStyle',':');
hold on

legend("$\mathrm{PB}1~(42a)$","$\mathrm{PB}2~(42b)$","$\mathrm{EB}~~(42c)$",'interpreter','latex');
xlim([0,15]);
ylim([0,12])
title('')
xlabel('$\mathrm{Timesteps}~{t}$','interpreter','latex','FontSize',12);
%xlabel('Total \# of Timesteps${(i\cdot\tau)}$','interpreter','latex','FontSize',12);
ylabel('${\frac{|{\hat{K}_t}-{K^*}|}{|{K^*}|}}$','interpreter','latex','FontSize',12);
% title("Convergence of Policy Gain K")
%sgtitle('$n_x=3,~n_u=2$','interpreter','latex','FontSize',14)




function output=vec(input)
% function vec: vectorize a matrix
    output=input(:);
end
function figuredata = Normfigure2(data1,data2,ref)
    %figuredata=(sqrt(sum(data1.*data1,1))-sqrt(sum(data2.*data2,1)))/norm(ref(:));
    figuredata=(sqrt(sum((data1-data2).*(data1-data2),1)))/norm(ref(:));
end