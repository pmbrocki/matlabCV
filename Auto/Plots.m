td  = (0:0.001:5);
A1= load('results_auto_PP','Xf','Sols');
A3= load('results_auto_PP3','Xf');
%F= load('results_finite_PP','sensf','y');
TITLE = {'Pau/Raup','Pau/Ral','Pau/Rvl','Pau/Ralp','Pau/Cau','Pau/Cal','Pau/Cvu','Pau/Cvl','Pau/Tsf','Pau/Trf','Pau/Ed','Pau/Es','Pau/Vd','Pau/Rmvop','Pau/Ravop'};


%%
% Sets figures to dock automatically
set(0,'DefaultFigureWindowStyle','docked')

%
% Set large font size for plots
jts_text_size =18;
originalAxFontSize = get(0,'DefaultAxesFontSize');
set(0,'DefaultAxesFontSize',18);

%
% Set linewidth to 2
jts_linewidth = 2;

% set default linewidth
set(0,'defaultlinelinewidth',2);

set(0,'defaultaxeslinewidth',1);

set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on','DefaultAxesZGrid','on');

% Default figure size is offset from the lower left corner and 800x600
originalFigurePosition = get(0,'DefaultFigurePosition');
set(0,'DefaultFigurePosition', [100,100,800,600]);

%% 
%Plots
for i =1:15
figure(i); clf;
plot(td',A1.Xf(:,i),td',A3.Xf(:,i)); %,td(1:end-1)',F.sensf(:,i));
title(TITLE(i));
xlabel('time (s)');
ylabel('Sensitivity');
legend('1','3')
%grid off

end