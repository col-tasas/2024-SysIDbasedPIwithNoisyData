load("RLS_PG_Boundednonvanish.mat");
load("RLS_PG_Boundedto0.mat");
load("RLS_PG_EnergyBounded.mat");

load("RLS_PI_Boundednonvanish.mat");
load("RLS_PI_Boundedto0.mat");
load("RLS_PI_EnergyBounded.mat");
% A=[-0.53,0.42,-0.44;0.42,-0.56,-0.65;-0.44,-0.65,0.35]; % system matrix A
% B=[0.43,-0.82;0.53,-0.78;0.26,0.40]; % input matrix B
% Q=[6.12,1.72,0.53;1.72,6.86,1.72;0.53,1.72,5.73]; % weighting matrix Q
% R=[1.15,-0.23;-0.23,3.62]; % weighting matrix R

[Pstar,Kstar,L,info]=idare(A,B,Q,R);

Pstarmatrix=kron(ones(1,episode),vec(Pstar));
Kstarmatrix=kron(ones(1,episode),vec(-Kstar));


% bounded to 0
P1=timeseries(Normfigure2(phistory1,Pstarmatrix,Pstar),0:1:episode-1);
K1=timeseries(Normfigure2(khistory1,Kstarmatrix,-Kstar),0:1:episode-1);
% bounded 
P2=timeseries(Normfigure2(phistory2,Pstarmatrix,Pstar),0:1:episode-1);
K2=timeseries(Normfigure2(khistory2,Kstarmatrix,-Kstar),0:1:episode-1);
% energy bounded
P3=timeseries(Normfigure2(phistory3,Pstarmatrix,Pstar),0:1:episode-1);
K3=timeseries(Normfigure2(khistory3,Kstarmatrix,-Kstar),0:1:episode-1);
% bounded to 0
P4=timeseries(Normfigure2(phistory,Pstarmatrix,Pstar),0:1:episode-1);
K4=timeseries(Normfigure2(khistory,Kstarmatrix,-Kstar),0:1:episode-1);
%
P5=timeseries(Normfigure2(phistory5,Pstarmatrix,Pstar),0:1:episode-1);
K5=timeseries(Normfigure2(khistory5,Kstarmatrix,-Kstar),0:1:episode-1);
%
P6=timeseries(Normfigure2(phistory4,Pstarmatrix,Pstar),0:1:episode-1);
K6=timeseries(Normfigure2(khistory4,Kstarmatrix,-Kstar),0:1:episode-1);

figure;
%subplot(2,1,1);

plot(P2,'linewidth',2,'color','blue','LineStyle','-');
hold on
plot(P1,'linewidth',2,'color','red','LineStyle','-.');
hold on
plot(P3,'linewidth',2,'color','m','LineStyle','--');
hold on
plot(P4,'linewidth',1.5,'color','blue','LineStyle','-.');
hold on
plot(P5,'linewidth',1.5,'color','red','LineStyle','--');
hold on
plot(P6,'linewidth',1.5,'color','m','LineStyle','-');
hold on


legend("$\mathrm{PB}1~(42a)~\mathrm{ORLS+PI}$","$\mathrm{PB}2~(42b)~\mathrm{ORLS+PI}$","$\mathrm{EB}~~(42c)~\mathrm{ORLS+PI}$","$\mathrm{PB}1~(42a)~\mathrm{ORLS+PG}$","$\mathrm{PB}2~(42b)~\mathrm{ORLS+PG}$","$\mathrm{EB}~~(42c)~\mathrm{ORLS+PG}$",'interpreter','latex');
xlim([0,70]);
ylim([0,5.5])
title('')
xlabel('$\mathrm{Timesteps}~{t}$','interpreter','latex','FontSize',12);
ylabel('${\frac{|{\hat{P}_t}-{P^*}|}{|{P^*}|}}$','interpreter','latex','FontSize',12);
%title("Convergence of Cost Function Kernel P")
%subplot(2,1,2)

% plot(K2,'linewidth',1.5,'color','blue');
% hold on
% plot(K1,'linewidth',1.5,'color','red','LineStyle','--');
% hold on
% plot(K3,'linewidth',1.5,'color','m','LineStyle',':');
% hold on
% 
% legend("$\mathrm{Instantaneous~Bounded}~1$","$\mathrm{Instantaneous~Bounded}~2$","$\mathrm{Energy~Bounded}$",'interpreter','latex');
% xlim([0,50]);
% ylim([0,5])
% title('')
% xlabel('$\mathrm{Timesteps}~{t}$','interpreter','latex','FontSize',12);
% %xlabel('Total \# of Timesteps${(i\cdot\tau)}$','interpreter','latex','FontSize',12);
% ylabel('${\frac{\|{\hat{K}_t}-{K^*}\|_F}{\|{K^*}\|_F}}$','interpreter','latex','FontSize',12);
% % title("Convergence of Policy Gain K")
% %sgtitle('$n_x=3,~n_u=2$','interpreter','latex','FontSize',14)
% 



function output=vec(input)
% function vec: vectorize a matrix
    output=input(:);
end
function figuredata = Normfigure2(data1,data2,ref)
    %figuredata=(sqrt(sum(data1.*data1,1))-sqrt(sum(data2.*data2,1)))/norm(ref(:));
    figuredata=(sqrt(sum((data1-data2).*(data1-data2),1)))/norm(ref(:));
end