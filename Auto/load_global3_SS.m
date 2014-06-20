% This function initializes the parameters for the model and sets initial
% values for the variables.
function [x0, Init, low, hi] = load_global3_SS %(ex)

global DIFF_INC ODE_TOL
global Rcl COd Vold

ODE_TOL  = 1e-8;
DIFF_INC = 1e-4; % SHOULD BE LARGER THAN SQRT(ODE_TOL)

H      = 186; % Subject height cm
W      = 83;  % Subject weight kg
Gender = 1;   % Male 1, female 2

% Define subject specific quantities
BSA = sqrt(H*W/3600);  % Body surface area

if Gender == 2
  TotalVol = (3.47*BSA - 1.954)*1000; % Female
elseif Gender == 1
  TotalVol = (3.29*BSA - 1.229)*1000; % Male
end;
TotFlow = TotalVol/60; % Assume blood circulates in one min
COd     = TotFlow*60/1000;

% Flows (related to subject)
qaup = TotFlow*0.90; % Upper body flow, arteries --> veins
qal  = TotFlow*0.10; % Upper body arteries --> lower body arteries
qalp = qal;          % Lower body flow, arteries --> veins
qvl  = qal;          % Lower body veins --> upper body veins

% Pressures (related to subject)     
pau    = 100;           % Upper body arteries
pvu    = 3.5;           % Upper body veins - Note higher venous pressure
pal    = pau*0.98;      % Lower body arteries
pvl    = 3.75;          % Lower body veins
Vlv    = 120;            % Left ventricular volume

% Resistances (Ohm's law)
Raup = (pau-pvu)/qaup; % Upper body peripheral resistance
Ral  = (pau-pal)/qal;  % Upper body arteries --> lower body arteries
Rvl  = (pvl-pvu)/qvl;  % Lower body veins --> upper body veins
Ralp = (pal-pvl)/qalp; % Lower body peripheral resistance

Rop = 0.001; %(max(pd)-max(pd)*0.9)/TotFlow;
Rcl = 40;    % Closed valve
Rmvop = Rop;
Ravop = Rop;
% Volumes (Beneken and deWit)
Vau = TotalVol*0.11;   % Arterial volume 11% of total volume
Vvu = TotalVol*0.66;   % Venous volume 66% of total volume
Val = TotalVol*0.02;   % Lower body arterial volum 2% of total volume
Vvl = TotalVol*0.06;   % Lower body venous volume 6% of total volume
                       % Pulmonary system 16% of total
                       % blood volume (not modeled).
                       % The heart contains about 5.5% of the total volume
Vold = Vau + Vvu + Val + Vvl;

disp('INITIAL:');
disp(strcat('Volume (ml): ',num2str(TotalVol*0.85,4)));
disp(strcat('CO  (l/min): ',num2str(COd,4)));
disp(' ');
disp('COMPUTED:');

% Compliances, stressed volume percentages from Beneken are weighted averages
Cau = (Vau*0.19)/pau;  % Upper body artery compliance (1.58-0.12)*mean(pd)/(585-75) = 0.19
Cvu = (Vvu*0.05)/pvu;  % Upper body venous compliance (219-38)/(2916-295) = 0.07
Cal = (Val*0.05)/pal;  % Lower body artery compliance 0.12*mean(pd)/75 = 0.11 (with mp=68)
Cvl = (Vvl*0.16)/pvl;  % Lower body venous compliance, the value in Benekin does not work

% Ventricle parameters 
%Average of all segments
Tsf = 0.12;     % Time fraction for systolic phase, increasing elasticity
Trf = 0.18;     % Time fraction for Systolic phase, decreasing elasticity
Ed     = 0.03;     % Minimum elasticity of the heart  [ 4 / (125-10) ]
Es     = 2;        % Maximum elasticity of the heart [ 120/(70-10) ]
Vd     = 10;       % Heart volume at end-diastolic pressure = 0


% Initial conditions
Init = [pau pvu pal pvl Vlv]; % SS BEFORE TILT

% Parameter vector
x0 = [Raup Ral Rvl Ralp ...      %1-4
      Cau Cal Cvu Cvl ...        %5-8
      Tsf Trf Ed Es Vd Rmvop Ravop];   %9-15
x0 = log(x0);

% load optim_par.mat
% x0 = x;
% Upper bound and lower bounds for optimization
hi  = x0+log(4); 
low = x0-log(4);