% Close everything open and and clear all variables
close all;
clear all;

% Load the data we will be optimizing against. We will only use the initial values from the data as well as the time for the final step
 load data_hiv;

odename = @cardio_rhs_ad;

% Write data to some variables
 DATA = data_hiv(:,2:end);
 TIMES = data_hiv(:,1);

rho = 6; % Dimension of ODE system
t0 = TIMES(1); % Initial time
tf = TIMES(end); % Final time
y0 = DATA(1,:);

% Create time points we will use the the analysis 
n=1000;
tvec = t0:tf/(n-1): tf;
times = tvec;

% Parameter values
lam1 = 1e+4;        d1 = 0.01;          epsilon = 0;        k1 = 8.0e-7;
lam2 = 31.98;       d2 = 0.01;          f = 0.34;           k2 = 1e-4;
delta = 0.7;        m1 = 1.0e-5;        m2 = 1.0e-5;        NT = 100;
c = 13;             rho1 = 1;           rho2 = 1;           lamE = 1;
bE = 0.3;           Kb = 100;           dE = 0.25;          Kd = 500;
deltaE = 0.1;

%
q   = [lam1,d1,epsilon,k1,lam2,d2,f,k2,delta,m1,m2,NT,c,rho1,rho2,lamE,bE,Kb,dE,Kd,deltaE];
npar = length(q);

options = odeset('RelTol',1e-8);

[t,Y] = ode15s(@cardio_rhs_ad,times,DATA(1,:),options,q);

k = 4;
indices = [2 8 9 17];
epsi=1e-2;
for i=1:4
  a = zeros(size(q));
  a(indices(i)) = 1;
  newpars = q+q.*a.*epsi;
  [t,Ysens] = ode15s(@hiv_rhs_ad,times,DATA(1,:),options,newpars);
  A = (Ysens-Y)./(q(indices(i))*epsi);
  A = [A(:,1); A(:,2); A(:,3); A(:,4); A(:,5); A(:,6)];
  fin_sens(:,i)=A;
end;
fin_sens = [fin_sens(:,1) fin_sens(:,2) fin_sens(:,3) fin_sens(:,4)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Sensitivity Equations %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial values for sensitivities are set to zero. States of original ODE has regular initial values
ysens0 = [zeros(24,1)', DATA(1,:)];
[t,ysens] = ode15s(@hiv_rhs_delta,tvec,ysens0,options,q);

% We are only interested in some of the equations. This is for the first state
Xsenseqdelta = squeeze(ysens(:,1));
Xsenseqd1 = squeeze(ysens(:,7));
Xsenseqk2 = squeeze(ysens(:,13));
XsenseqbE = squeeze(ysens(:,19));

% Here the next states are added to the matrix corresponding to each parameter
for i = 2:rho
    Xsenseqdelta = [Xsenseqdelta; squeeze(ysens(:,i))];
    Xsenseqd1 = [Xsenseqd1; squeeze(ysens(:,1*rho+i))];
    Xsenseqk2 = [Xsenseqk2; squeeze(ysens(:,2*rho+i))];
    XsenseqbE = [XsenseqbE; squeeze(ysens(:,3*rho+i))];

end
%Finally one big matrix is created with all the sensitivities
Xsenseq = [Xsenseqd1 Xsenseqk2 Xsenseqdelta XsenseqbE];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Automatic differentiation %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Same parameters above, but listed for visual selection of parameters for sensitivity analysis
q    = [lam1,d1,epsilon,k1,lam2,d2,f,k2,delta,m1,m2,NT,c,rho1,rho2,lamE,bE,Kb,dE,Kd,deltaE];
qsel = [   0, 1,      0, 0,   0, 0,0, 1,    1, 0, 0, 0,0,   0,   0,   0, 1, 0, 0, 0,     0]; 

% Number of parameters to do sensitivity analysis for
p=sum(qsel);

tsoptions = []; %  Options for tssolve's ode solver. Ignore for now.

%%% Perform the AD calculations. See tssolve.m for structure of Dq and Dx0

[Dq Dx0] = tssolve(odename,y0,q,tvec,1,tsoptions);

%%% Put the sensitivities into the form of the problem
X = squeeze(Dq.regsens(:,1,:));
for i = 2:rho
    X = [X; squeeze(Dq.regsens(:,i,:))];
end

%%%  Cut out the columns of X corresponding to parameters we're not
%%%  interested in.
qsel = repmat(qsel,n*rho,1);
X = reshape(X(qsel ==1),n*rho,p);

for i=1:6
  figure(i)
  plot(times,X((i-1)*n+1:i*n,1),'-k');
  hold on
  plot(times,Xsenseq((i-1)*n+1:i*n,1),'--k');
  plot(times,fin_sens((i-1)*n+1:i*n,1),'-.k');
  legend('AD', 'Sens EQ', 'Fin. Dif');
end;

for i=1:6
  figure(i+6)
  plot(times,X((i-1)*n+1:i*n,2),'-k');
  hold on
  plot(times,Xsenseq((i-1)*n+1:i*n,2),'--k');
  plot(times,fin_sens((i-1)*n+1:i*n,2),'-.k');
  legend('AD', 'Sens EQ', 'Fin. Dif');
end;
for i=1:6
  figure(i+12)
  plot(times,X((i-1)*n+1:i*n,3),'-k');
  hold on
  plot(times,Xsenseq((i-1)*n+1:i*n,3),'--k');
  plot(times,fin_sens((i-1)*n+1:i*n,3),'-.k');
  legend('AD', 'Sens EQ', 'Fin. Dif');
end;
for i=1:6
  figure(i+18)
  plot(times,X((i-1)*n+1:i*n,4),'-k');
  hold on
  plot(times,Xsenseq((i-1)*n+1:i*n,4),'--k');
  plot(times,fin_sens((i-1)*n+1:i*n,4),'-.k');
  legend('AD', 'Sens EQ', 'Fin. Dif');
end;
covavs = inv(X'*X);
covsenseq = inv(Xsenseq'*Xsenseq);
covfin = inv(fin_sens'*fin_sens);

%%%% The code used to produce the graphs showed in the paper handin. 
%i=1;
%h1 = figure(1)
%set(gca,'FontSize',14);
%plot(times,X((i-1)*n+1:i*n,1),'-k','LineWidth',2);
%hold on
%plot(times,Xsenseq((i-1)*n+1:i*n,1),'--k','LineWidth',2);
%plot(times,fin_sens((i-1)*n+1:i*n,1),'-.k','LineWidth',2);
%legend('AD', 'Sens EQ', 'Fin. Dif');
%title('Sensitivity of T1 to d1')
%
%i =2 ;
%h2 = figure(2)
%set(gca,'FontSize',14);
%plot(times,X((i-1)*n+1:i*n,2),'-k','LineWidth',2);
%hold on
%plot(times,Xsenseq((i-1)*n+1:i*n,2),'--k','LineWidth',2);
%plot(times,fin_sens((i-1)*n+1:i*n,2),'-.k','LineWidth',2);
%legend('AD', 'Sens EQ', 'Fin. Dif');
%title('Sensitivity of T2 to k2')
%
%print(h1,'-dpng',['../sensi1.png'])
%print(h2,'-dpng',['../sensi2.png'])
