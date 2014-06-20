
td  = [0:0.001:15];
%td = [td;td;td];
load('results_auto_PP','Xf');
load('results_direct_PP','S');
load('results_finite_PP','sensfinal');

figure(6); %clf;
h=plot(td',Xf(:,1),td',S(:,1),td(:,1:end-1)',sensfinal(:,1));
set(h,'Linewidth',2);
set(gca,'Fontsize',20);
title('Pau/Raup');
xlabel('time (s)');
ylabel('Sensitivity');
legend('Auto','Direct','FD')
grid on