%Solves model sensitivities using Automatic Diff.

ODE_TOL  = 1e-8;

%n=1001;
[q, x0] = load_global3_SS; 
q = exp(q); % Note load_global returns the log-scaled parameters

Raup = q(1); % Peripheral resistance upper body
Ral  = q(2); % Resistance between upper and lower body arteries
Rvl  = q(3); % Resistance between upper and lower body veins (open valve)
Ralp = q(4);

Cau = q(5); % Upper body arterial compliance
Cal = q(6); % Arterial lower body compliance
Cvu = q(7); % Venous upper body compliance
Cvl = q(8); % Venous lower body compliance

Tsf  = q(9);
Trf  = q(10);
Ed   = q(11); % Minimum contractility left ventricle
Es   = q(12); % Maximum contractility left ventricle
Vd   = q(13); % Dead-space volume left ventricle
Rmvop = q(14);
Ravop = q(15);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Automatic differentiation %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


T = 1;
td  = (0:0.001:5);
%yf = CVmodel_sens(q,td,x0);
k1 = 1;
k2 = 1001;

tsoptions = []; %   Options for tssolve's ode solver. Ignore for now.
Xf = [];        %   Sensitiviy Matrix 
Sols = [];      %   ODE output
m=16;
n=5;
sens0 = (zeros(m*n,1));
for j=1:5 %15 %Loop for each pulse
    tdc = td(k1:k2);
    T   = tdc(end)-tdc(1);
    disp([tdc(1) tdc(end) T]);
    q(16)= tdc(1);
    %%% Perform the AD calculations. See tssolve.m for structure of Dq and Dx0
    % tdc = tdc-tdc(1);
    [Dq,Dx0] = tssolve3(@cardio_rhs_ad,x0,q,sens0,tdc,1,tsoptions);
   
   
    
    %Current sensitivites - move data to 2D matrix
    X = squeeze(Dq.regsens(:,1,:));     %wrt to Parameters
    Z = squeeze(Dx0.regsensx0(:,1,:));  %wrt to IC
    Y = squeeze(Dq.odeout(:,1));        %ODE output
    for i = 2:n
        X = [X, squeeze(Dq.regsens(:,i,:))]; 
        Z = [Z, squeeze(Dx0.regsensx0(:,i,:))]; 
        Y = [Y, squeeze(Dq.odeout(:,i))];
    end
    
   
    
    Xf = [Xf;X(1:end-1,:)]; %Total sensitivites
    Sols = [Sols;Y(1:end-1,:)];
    
    % Give the final values as the new initial condition
    sens0 = Dq.regsens(end,:,:);
    
    
    x0 = [Dq.odeout(end,1) Dq.odeout(end,2) Dq.odeout(end,3) Dq.odeout(end,4) Dq.odeout(end,5)];

    k1 = k2;
    k2 = k2+1000;
end


Xf = [Xf ; X(end,:)];           %Sens
Sols = [Sols ; Y(end,:)];       %ODE output

save results_auto_PP3.mat;


% figure(6)
% h=plot(td(:,1:end-1)',Xf(:,(1:15)));
% set(h,'Linewidth',2);
% set(gca,'Fontsize',20);
% xlabel('time (s)');
% ylabel('Pau');
% legend('Raup','Ral','Rvl','Ralp','Cau','Cal','Cvu','Cvl','Tsf','Trf','Ed','Es','Vd','Rmvop', 'Ravop')
% grid on;
% 
% figure(7)
% h=plot(td(:,1:end-1)',Xf(:,(17:31)));
% set(h,'Linewidth',2);
% set(gca,'Fontsize',20);
% xlabel('time (s)');
% ylabel('Pvu');
% legend('Raup','Ral','Rvl','Ralp','Cau','Cal','Cvu','Cvl','Tsf','Trf','Ed','Es','Vd','Rmvop', 'Ravop')
% grid on;
% 
% figure(8)
% h=plot(td(:,1:end-1)',Xf(:,(33:47)));
% set(h,'Linewidth',2);
% set(gca,'Fontsize',20);
% xlabel('time (s)');
% ylabel('Pvu');
% legend('Raup','Ral','Rvl','Ralp','Cau','Cal','Cvu','Cvl','Tsf','Trf','Ed','Es','Vd','Rmvop', 'Ravop')
% grid on;
% 
% figure(9)
% h=plot(td(:,1:end-1)',Xf(:,(49:63)));
% set(h,'Linewidth',2);
% set(gca,'Fontsize',20);
% xlabel('time (s)');
% ylabel('Pvu');
% legend('Raup','Ral','Rvl','Ralp','Cau','Cal','Cvu','Cvl','Tsf','Trf','Ed','Es','Vd','Rmvop', 'Ravop')
% grid on;
% 
% figure(10)
% h=plot(td(:,1:end-1)',Xf(:,(65:79)));
% set(h,'Linewidth',2);
% set(gca,'Fontsize',20);
% xlabel('time (s)');
% ylabel('Pvu');
% legend('Raup','Ral','Rvl','Ralp','Cau','Cal','Cvu','Cvl','Tsf','Trf','Ed','Es','Vd','Rmvop', 'Ravop')
% grid on;

% figure(11);clf; % Plot pressure upper body arteries model and data
% h=plot(td,pauS);
% set(h,'Linewidth',2);   %2 indicates how tick a line you want
% set(gca,'Fontsize',20); %20 is the fontsize, gca is "get current axis"
% xlabel('time (s)');     %add xlabel
% ylabel('pau (mmHg)');   %add ylabel
% grid on;                %add grid
% 
% figure(12);clf; % Plot pressure upper body arteries model and data
% h=plot(td,VlvS);
% set(h,'Linewidth',2);   %2 indicates how tick a line you want
% set(gca,'Fontsize',20); %20 is the fontsize, gca is "get current axis"
% xlabel('time (s)');     %add xlabel
% ylabel('Vlv (ml)');   %add ylabel
% grid on;                %add grid





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Finite Differences %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% epsi=1e-4; %Step size
% fin_sens = []; 
% for i=1:15
%   a = zeros(size(q));
%   a(i) = 1;
%   newq = q+q.*a.*epsi;
%   [t,Ysens] = ode15s(@cardio_rhs_ad,td,x0,options,newq,tdc(1),T);
%   A = (Ysens-Y')./(q(i)*epsi);
%   %A  = [A(:,1); A(:,2); A(:,3); A(:,4); A(:,5)];
%   fin_sens = [fin_sens A];
% end;
% 
% %td = [td; td; td; td; td; td; td; td; td; td; td; td; td; td; td];
% %hold on 
% for i=1:15;
% figure(1);clf;
% h=plot(td',fin_sens(:,1), td',fin_sens(:,6), td',fin_sens(:,11),td',fin_sens(:,16),td',fin_sens(:,21),td',fin_sens(:,26));
% set(h,'Linewidth',2);
% set(gca,'Fontsize',20);
% xlabel('time (s)');
% ylabel('Pau Sens ');
% %legend('Raup','Ral','Rvl','Ralp','Cau','Cal','Cvu','Cvl','Tsf','Trf','Ed','Es','Vd','Rmvop', 'Ravop')
% grid on;
% end

