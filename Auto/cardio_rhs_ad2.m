function dy = cardio_rhs_ad2(t,y,pars)
T=1;
global Rcl Rop beta 
%%% Locations of our states in y. AD works better with direct references to
%%% vector, not with these types of renames
% Define variables
% pau    = y(1);
% pvu    = y(2);
% pal    = y(3); 
% pvl    = y(4); 
% Vlv    = y(5);

Raup = pars(1); % Peripheral resistance upper body
Ral  = pars(2); % Resistance between upper and lower body arteries
Rvl  = pars(3); % Resistance between upper and lower body veins (open valve)
Ralp = pars(4);

Cau = pars(5); % Upper body arterial compliance
Cal = pars(6); % Arterial lower body compliance
Cvu = pars(7); % Venous upper body compliance
Cvl = pars(8); % Venous lower body compliance

Tsf  = pars(9);
Trf  = pars(10);
Ed   = pars(11); % Minimum contractility left ventricle
Es   = pars(12); % Maximum contractility left ventricle
Vd   = pars(13); % Dead-space volume left ventricle
ts = pars(14);

% Calculate length of current cardiac cycle and time elapsed in the cycle
Elv  = ElastanceBasic(t-ts,T,Trf,Ed,Es,Tsf); % Ventriular Elastance function from Elastance1.m
plv  = Elv.*(y(5)-Vd); % Stressed ventricular volume & Vd(consant), ventricular volume at zero dyastolic pressure

% Calculate lower venous, aortic and mitral valve resistances

 Rmv  = Rcl - (Rcl-Rop) ./(1+exp(-beta*(y(2)-plv)));
 Rav  = Rcl - (Rcl-Rop) ./(1+exp(-beta*(plv-y(1))));

%Calculate flows
% qav  = (plv - y(1))./Ravop;
% qmv  = (y(2) - plv)./Rmvop;
% qvl  = (y(4) - y(2))./Rvl;
% qal  = (y(1) - y(3))./Ral;
% qaup = (y(1) - y(2))./Raup;
% qalp = (y(3) - y(4))./Ralp;


 dy =  [(((plv - y(1))./Rav) - ((y(1) - y(3))./Ral) - ((y(1) - y(2))./Raup))./Cau; ...   %pau 
        (((y(4) - y(2))./Rvl) + ((y(1) - y(2))./Raup) - ((y(2) - plv)./Rmv))./Cvu; ...   %pvu 
        (((y(1) - y(3))./Ral) - ((y(3) - y(4))./Ralp))./Cal; ...                           %pal
        (((y(3) - y(4))./Ralp) - ((y(4) - y(2))./Rvl))./Cvl; ...                           %pvl
        (((y(2) - plv)./Rmv) - ((plv - y(1))./Rav));];                                 %Vlv
