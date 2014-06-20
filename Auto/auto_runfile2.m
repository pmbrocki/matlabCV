%Solves model sensitivities using Automatic Diff.

ODE_TOL  = 1e-8;

n=1001;
[q, x0] = load_global2_SS; 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Automatic differentiation %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


T = 1;
td  = [0:0.001:15];

k1 = 1;
k2 = 1001;

tsoptions = []; %  Options for tssolve's ode solver. Ignore for now.
Xf = [];
for j=1:15 %Loop for each pulse
    tdc = td(k1:k2);
    T   = tdc(end)-tdc(1);
    disp([tdc(1) tdc(end) T]);
    q(14)= tdc(1);
    %%% Perform the AD calculations. See tssolve.m for structure of Dq and Dx0
    % tdc = tdc-tdc(1);
    [Dq,Dx0] = tssolve(@cardio_rhs_ad2,x0,q,tdc,1,tsoptions);
    
    rho = 5;
    
    %Current sensitivites - move data to 2D matrix
    X = squeeze(Dq.regsens(:,1,:));
    for i = 2:rho  
        X = [X, squeeze(Dq.regsens(:,i,:))]; 
    end
    
    Xf = [Xf;X(1:end-1,:)]; %Total sensitivites
    
    
    % Solve ODE's
    options=odeset('RelTol',ODE_TOL, 'AbsTol',ODE_TOL); %,'MaxStep',0.001);
    sol = ode15s(@cardio_rhs_ad2,tdc,x0,options,q);
    sols= deval(sol,tdc);
    
    % Solutions
    pau = sols(1,:)'; % Arterial pressure upper body
    pvu = sols(2,:)'; % Venous pressure upper body
    pal = sols(3,:)'; % Arterial pressure lower body
    pvl = sols(4,:)'; % Venous pressure lower body
    Vlv = sols(5,:)'; % Volume left ventricle
    
    % Factors used for scaling compliances are related to stressed vs.
    % unsterssed volume - definition can be found in load_global
    Vau  = Cau*pau; % Volume upper body arteries
    Val  = Cal*pal; % Volume lower body arteries
    Vvl  = Cvl*pvl; % Volume lower body veins
    Vvu  = Cvu*pvu; % Volume upper body veins
    Vtot = Vau + Vvu + Val + Vvl + Vlv; % Total volume in the system
    
    for k = 1:length(tdc)
        Elv(k) = ElastanceBasic(tdc(k)-tdc(1),T,Trf,Ed,Es,Tsf);
    end;
    plv    = Elv'.*(Vlv-Vd); % Pressure left ventricle
   % Rmvop  = Rop; %Rcl - (Rcl-Rop) ./(1+exp(-beta*(pvu-plv)));
   % Ravop  = Rop; %Rcl - (Rcl-Rop) ./(1+exp(-beta*(plv-pau)));
    
    qav  = (plv - pau)./Ravop;
    qmv  = (pvu - plv)./Rmvop;
    qvl  = (pvl - pvu)./Rvl;      % Flow between lower and upper body veins
    qal  = (pau - pal)./Ral;      % Flow upper body to lower body arteries
    qaup = (pau - pvu)./Raup;     % Flow through upper body peripheral circulation
    qalp = (pal - pvl)./Ralp;     % Flow through lower body peripheral circulation

    
    
    
    x0 = [pau(end) pvu(end) pal(end) pvl(end) Vlv(end)];
    
    k1 = k2;
    k2 = k2+1000;
end
Xf = [Xf ; X(end,:)];

save results_auto_PP2.mat;


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

