function dy = hiv_rhs_ad(t,y,parameters)

%%% Locations of our states in y. AD works better with direct references to
%%% vector, not with these types of renames
% T1 = y(1,1);
% T2 = y(2,1);
% T1s = y(3,1);
% T2s = y(4,1);
% V = y(5,1);
% E = y(6,1);

  lam1 = parameters(1);
  d1 = parameters(2);
  epsilon = parameters(3);
  k1 = parameters(4);
  lam2 = parameters(5);
  d2 = parameters(6);
  f = parameters(7);
  k2 = parameters(8);
  delta = parameters(9);
  m1 = parameters(10);
  m2 = parameters(11);
  NT = parameters(12);
  c = parameters(13);
  rho1 = parameters(14);
  rho2 = parameters(15);
  lamE = parameters(16);
  bE = parameters(17);
  Kb = parameters(18);
  dE = parameters(19);
  Kd = parameters(20);
  deltaE = parameters(21);

  dy = [lam1 - d1*y(1) - (1 - epsilon)*k1*y(5)*y(1);
        lam2 - d2*y(2) - (1 - f*epsilon)*k2*y(5)*y(2);
        (1 - epsilon)*k1*y(5)*y(1) - delta*y(3) - m1*y(6)*y(3);
        (1 - f*epsilon)*k2*y(5)*y(2) - delta*y(4) - m2*y(6)*y(4);
        NT*delta*(y(3)+y(4)) - c*y(5) - ((1 - epsilon)*rho1*k1*y(1) + (1 - f*epsilon)*rho2*k2*y(2))*y(5);
        lamE + (bE*(y(3)+y(4))./(y(3)+y(4)+Kb))*y(6) - (dE*(y(3)+y(4))./(y(3)+y(4)+Kd))*y(6) - deltaE*y(6)];

