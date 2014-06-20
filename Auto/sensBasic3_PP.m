% This file contains the differential equations of the model and also
% calculates the new left ventricular pressure by calling the function
% Elastance.

function xdot = sensBasic3_PP(t,p,pars,ts,T)

% Define variables
pau    = p(1);
pvu    = p(2);
pal    = p(3); 
pvl    = p(4); 
Vlv    = p(5);

% Define resistances
Raup  = pars(1); 
Ral   = pars(2);
Rvl   = pars(3);
Ralp  = pars(4);

% Define compliances
Cau = pars(5); 
Cal	= pars(6);
Cvu	= pars(7);
Cvl = pars(8); 

% Heart parameters
Tsf = pars(9);
Trf = pars(10);
Ed  = pars(11); 
Es  = pars(12);
Vd  = pars(13);
Rmvop = pars(14);
Ravop = pars(15);

% Calculate length of current cardiac cycle and time elapsed in the cycle
Elv  = ElastanceBasic(t-ts,T,Trf,Ed,Es,Tsf); % Ventriular Elastance function from Elastance1.m
plv  = Elv*(Vlv-Vd); % Stressed ventricular volume & Vd(consant), ventricular volume at zero dyastolic pressure

%beta = 10; %Smoothness parameter of valve resistance 
% Calculate lower venous, aortic and mitral valve resistances
%Rmv = Rop ;%Rcl - (Rcl-Rop) ./(1+exp(-beta*(pvu-plv )));
%Rav = Rop ;%Rcl - (Rcl-Rop) ./(1+exp(-beta*(plv-pau)));

%Calculate flows
qav  = (plv - pau)/Ravop;
qmv  = (pvu - plv)/Rmvop;
qvl  = (pvl - pvu)/Rvl;
qal  = (pau - pal)/Ral;
qaup = (pau - pvu)/Raup;
qalp = (pal - pvl)/Ralp;


% Calculation of pressures, volume, average pressures, Ralp and Cvl
xdot = [(qav - qal - qaup)/Cau; ...   %pau 
        (qvl + qaup - qmv)/Cvu; ...   %pvu 
        (qal - qalp)/Cal; ...         %pal
        (qalp - qvl)/Cvl; ...         %pvl
        (qmv - qav);];                %Vlv
    
        

     
  