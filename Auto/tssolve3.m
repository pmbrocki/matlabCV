function [RS,RSx0,GS,GSx0]=tssolve3(odename,x0,q,sens0,tint,sensflag,options)

% TSSOLVE AD Sensitivity Solver.
%   This code will compute regular, relative, and if possible general
%   sensitivities with respect to parameters and initial conditions for a
%   model defined by a set of ODEs of the form x'=f(t,x(q),q;x0). This code
%   requires Martin Fink's myAD code in your working path.
%
%   Usage:
%       [RS,RSx0,GS,GSx0]=tssolve(odename,x0,q,tint,sensflag,options)
%
%   Inputs:
%       odename: a string or fhandle of the right-hand-side ode file.
%       x0: initial condition for the states.
%       q: the parameter value at which the sens analysis is performed.
%       x0s: denotes which x0s to estimate
%       qs: denotes parameters to estimate
%       tint: the interval, or time values that the ode will be solved
%       sensflag: chooses how the relative sensitivities are computed. if
%           sensflag==1, then sens=dx/dq .* q/x,
%           sensflag==2, then sens=dx/dq .* q/max(x),
%           sensflag=else, then sens=dx/dq .* q.
%       options: any ode options you would like the integrator to use.
%
%   Output:
%       RS: structure describing the regular and relative sensitivies.
%       RS.t: time vector.
%       RS.odeout: straight output from the ODE integrator.
%       RS.regsens and RS.relsens is ordered in 3D arrays as  
%       (time, state, parameter), so the sensitivity of state one with
%       parameter 5 would be regsens(:,1,5).
%
%       RSx0: structure containing the regular and relative sensitivities 
%       of the states wrt to initial conditions. 
%       RSx0.regsensx0 again is a 3D array ordered by 
%       (time,state,initial condition). 
%       RS.trimmed are the sensitivies, nonscaled and ordered from the
%       integrator.
%
%       GS: structure containing the generalized sensitivities.
%       GS.gensens are the generalized sensitivities in a 3D array ordered
%       as (time, state, parameter). Be careful on the accuracy given the
%       rank of the information matrix.
%       GS.infomat: the information matrix for a given state: (:,:,state).
%       GS.relinfomat: the relative information matrix for a given state
%       (:,:,state).
%
%       GSx0: structure containing the generalized sensitivities with
%       respect to initial conditions. Ordering is identical to the RSx0
%       structure.
% 
% Altered by RTG
%       
% Written and maintained by Adam Attarian (arattari@unity.ncsu.edu).
%
% Inspiration originally from John David's totalsenssolver2 code from years
% past.
%
% This code comes with no guarantee or warranty of any kind.

n=length(x0); % number of states
m=length(q);  % number of parameters

% ensure col vector inputs
x0=x0(:);
q=q(:);
sensX0 = reshape(diag(ones(n,1)),n^2,1);
sens0 = reshape(sens0,m*n,1) ;

if nargin<6; 
    options=[];
end

[t,sens]=ode15s(@odewsensAD,tint,[x0; sens0;sensX0],options,q,odename);
N=length(t);

RelSens=zeros(N,n,m);
RegSens=zeros(N,n,m);

for i=1:n
    for j=1:m
        % this orders as (:,state,param).
        RegSens(:,i,j)=sens(:,i+j*n);
        
        % need to compute the relative sensitivities.
        
        % first, compute the semirelative sensitivies.
        RelSens(:,i,j)=sens(:,i+j*n).*q(j);
        % compute the relative sensitivities if asked for.
        if sensflag==1
            RelSens(:,i,j)=RelSens(:,i,j)./sens(:,i);
        % compute modified rel sens if div by 0 would happen.
        elseif sensflag==2
            RelSens(:,i,j)=sens(:,i+j*n).*(q(j)./max(sens(:,i)));
        end
    end
end

RS.t=t; % save the time vector
RS.odeout=sens; % save the regular sens for output
RS.regsens=RegSens;
RS.relsens=RelSens;

% moving onto computing the sensitivites with respect to initial condition.
% man, the rest really is just commentary.
if nargout>=2
    strim=[sens(:,1:n) sens(:,n+n*m+1:end)]; % trimming to be just the x0 states
    RelSensx0=zeros(N,n,n);
    RegSensx0=zeros(N,n,n);
    
    % same ordering, folks.
    
    for i=1:n
        for j=1:n
            RegSensx0(:,i,j)=strim(:,i+j*n);
            
            % and now relativise.
            RelSensx0(:,i,j)=strim(:,i+j*n).*x0(j);
            if sensflag==1
                RelSensx0(:,i,j)=RelSens(:,i,j)./strim(:,i);
            elseif sensflag==2
                RelSensx0(:,i,j)=strim(:,i+j*n).*(x0(j)./max(strim(:,i)));
            end
%             disp(['j= ' num2str(j)])
        end
%         disp(['i= ' num2str(i)])
    end
    
    RSx0.t=t;
    RSx0.trimmed=strim;
    RSx0.regsensx0=RegSensx0;
    RSx0.relsensx0=RelSensx0;
end

if nargout>=3
    % now we're going to attempt to solve for the generalized sensitivities.
    BigGS=zeros(N,n,m);
    InfoOut=zeros(m,m,n);
    RelInfoMat=zeros(m,m,n);
    
    for i=1:n
        
        A=reshape(RegSens(:,i,:),N,m);
        S0=A'*A; % this is the regular infomation matrix for this state.
        
        B=reshape(RelSens(:,i,:),N,m);
        S0Rel=B'*B; % this is the relative information matrix for this state.
        
        % using cumsum to sum across the second dimension, represents doing
        % incremental updates to the gensens function as time progresses.
        % if there is linear dependance in parameters the info matrix will
        % be rank deficient, and we switch to using the psudeo inverse. i'm
        % assuming a well-conditioned rank test.
        
        if rank(S0)==m
            gs=cumsum(A'.*(S0\A'),2);  
        else
            fprintf('rank deficiency in gensens comp, using pinv(). results may be inaccurate.\n');
            gs=cumsum(A'.*(pinv(S0)*A'),2);
        end
        
        % (time, state, parameter)
        BigGS(:,i,:)=gs';
        
        % (param,param,state);
        InfoOut(:,:,i)=S0;
        RelInfoMat(:,:,i)=S0Rel;
    end
    
    GS.t=t;
    GS.gensens=BigGS;
    GS.infomat=InfoOut;
    GS.relinfomat=RelInfoMat;
end

if nargout==4
    % now we're going to attempt to solve for the generalized sensitivities.
    BigGSx0=zeros(N,n,n);
    InfoOutx0=zeros(n,n,n);
    RelInfoMatx0=zeros(n,n,n);
    
    for i=1:n
        
        A=reshape(RegSensx0(:,i,:),N,n);
        S0=A'*A; % this is the regular infomation matrix for this state.
        
        B=reshape(RelSensx0(:,i,:),N,n);
        S0Rel=B'*B; % this is the relative information matrix for this state.
        
        if rank(S0)==m
            gs=cumsum(A'.*(S0\A'),2);  
        else
            fprintf('rank deficiency in gensens comp, using pinv(). results may be inaccurate.\n');
            gs=cumsum(A'.*(pinv(S0)*A'),2);
        end
        
        % (time, state IC, state IC)
        BigGSx0(:,i,:)=gs';
        
        % (state IC,state IC,state);
        InfoOutx0(:,:,i)=S0;
        RelInfoMatx0(:,:,i)=S0Rel;
    end
    
    GSx0.t=t;
    GSx0.gensens=BigGSx0;
    GSx0.infomat=InfoOutx0;
    GSx0.relinfomat=RelInfoMatx0;
end



%--------------------------------------------------------------------------
% (AD) ODE function that when solved will give us all of the requisite
% sensitivies to paramaters and IC

function dx=odewsensAD(t,x,q,odename)

m=length(q);
n=.5*((-1-m)+sqrt((1+m)^2+4*length(x)));
dx=zeros(n+n*m+n^2,1); %(length(x),1) will also work

% evaluate the ODE, put the state solutions as the first n entries of dx

dx(1:n)=feval(odename,t,x(1:n),q);

% init states
x0=x(1:5);
q0=q;

% construct the AD variables
bigx=myAD(x0);
bigq=myAD(q0);

% evaluate the ODE with the AD variables and extract the derivatives, which
% are in matrix form (jacobians).
result1=feval(odename,t,bigx,q0);
result2=feval(odename,t,x0,bigq);

dfdx=getderivs(result1);
dfdq=getderivs(result2);



% pull the sensitivity equations out, excluding the states, convert to
% matrix. remember that reshape orders by columns, but in the sorting at
% the top i transpose stuff so that the ordering is different... that may
% be a problem.
sens=reshape(x(n+1:n+n*m),n,m);


% pull the IC sensitivities out. they are stuck at the end of the x vector.
sensx0=reshape(x(n+n*m+1:end),n,n);

% compute the sensitivities. this computation is from the chain rule. 
% seriously, it's real calculus. "sens" is the state that is being computed 
% via the ODE solver.
dsens=dfdx*sens+dfdq;
dsensx0=dfdx*sensx0;

% construct the non-states output. ordered as [states; q-sensitivities; IC sens]

dx(n+1:n+n*m)=dsens(:);
dx(n+n*m+1:end)=dsensx0(:);


    





