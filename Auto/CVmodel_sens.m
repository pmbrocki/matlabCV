% -------------------------------------- %
% This function solves the equations specified
% in the cvrhs.m file by using the built-in ode
% solver called ode15s. It also saves the solutions and
% the time-vector in a file called results.mat.
% The last part of the file calculates flows if needed.
% -------------------------------------- %
% Function when not collecting flow-information and optimizing
function  [yf]  = CVmodel_sens(x0,td,Init)

global ODE_TOL
yf =[];
%pars = exp(x0);
k1 = 1;
k2 = 1001;
for i = 1:15
    tdc = td(k1:k2);
    T   = tdc(end)-tdc(1);
    disp([tdc(1) tdc(end) T]);
    
    %Solve ODE's the usual way
    options = odeset('RelTol',ODE_TOL,'AbsTol',ODE_TOL); % Tolerances
    sol = ode15s(@sensBasic3_PP,tdc,Init,options,x0,tdc(1),T);    % Solve ODE's
    t   = sol.x;
    y   = deval(sol,tdc);
    
    yf = [yf y(:,1:end-1)];
    
    
    pau = y(1,:)'; % Arterial pressure upper body
    pvu = y(2,:)'; % Venous pressure upper body
    pal = y(3,:)'; % Arterial pressure lower body
    pvl = y(4,:)'; % Venous pressure lower body
    Vlv = y(5,:)'; % Volume left ventricle
    
%     for j = 1:length(tdc)
%         Elv(j) = ElastanceBasic(tdc(j)-tdc(1),T,Trf,Ed,Es,Tsf);
%     end;
%    plv  = Elv'.*(Vlv-Vd); % Pressure left ventricle
    
    
    Init = [pau(end) pvu(end) pal(end) pvl(end) Vlv(end)];
    
    k1 = k2;
    k2 = k2+1000;
    
end

yf = [yf y(:,end)];