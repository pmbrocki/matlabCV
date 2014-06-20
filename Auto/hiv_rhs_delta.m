function dy = hiv_rhs_delta(t,Y,optparameters,p)

  
  dT1=Y(1,1);
  dT2=Y(2,1);
  dT1s=Y(3,1);
  dT2s=Y(4,1);
  dV=Y(5,1);
  delE=Y(6,1);
  d1T1=Y(7,1);
  d1T2=Y(8,1);
  d1T1s=Y(9,1);
  d1T2s=Y(10,1);
  d1V=Y(11,1);
  d1E=Y(12,1);
  k2T1=Y(13,1);
  k2T2=Y(14,1);
  k2T1s=Y(15,1);
  k2T2s=Y(16,1);
  k2V=Y(17,1);
  k2E=Y(18,1);
  bET1=Y(19,1);
  bET2=Y(20,1);
  bET1s=Y(21,1);
  bET2s=Y(22,1);
  bEV=Y(23,1);
  bEE=Y(24,1);

  lam1 = optparameters(1);
  d1 = optparameters(2);
  epsilon = optparameters(3);
  k1 = optparameters(4);
  lam2 = optparameters(5);
  d2 = optparameters(6);
  f = optparameters(7);
  k2 = optparameters(8);
  delta = optparameters(9);
  m1 = optparameters(10);
  m2 = optparameters(11);
  NT = optparameters(12);
  c = optparameters(13);
  rho1 = optparameters(14);
  rho2 = optparameters(15);
  lamE = optparameters(16);
  bE = optparameters(17);
  Kb = optparameters(18);
  dE = optparameters(19);
  Kd = optparameters(20);
  deltaE = optparameters(21);

  
T1=Y(25,1);
T2=Y(26,1);
T1s=Y(27,1);
T2s=Y(28,1);
V=Y(29,1);
E=Y(30,1);
  
  
%delta  
dy(1,1)=-d1*dT1-k1*(T1*dV+V*dT1);
dy(2,1)=-d2*dT2-k2*(dV*T2+V*dT2);
dy(3,1)=k1*(dV*T1+V*dT1)-T1s-delta*dT1s-m1*(delE*T1s+E*dT1s);
dy(4,1)=k2*(dV*T2+V*dT2)-T2s-delta*dT2s-m2*(delE*T2s+E*dT2s);
dy(5,1)=NT*((T1s+T2s)+delta*(dT1s+dT2s))-c*dV-(rho1*k1*dT1+rho2*k2*dT2)*V-(rho1*k1*T1+rho2*k2*T2)*dV;
dy(6,1)=bE*(((dT1s+dT2s)*E/(T1s+T2s+Kb))+((T1s+T2s)*delE/(T1s+T2s+Kb))-((T1s+T2s)*E*(dT1s+dT2s)/((T1s+T2s+Kb)^2)))-dE*(((dT1s+dT2s)*E/(T1s+T2s+Kd))+((T1s+T2s)*delE/(T1s+T2s+Kd))-((T1s+T2s)*E*(dT1s+dT2s)/((T1s+T2s+Kd)^2)))-deltaE*delE;

%d1
dy(7,1)=-T1-d1*d1T1-k1*d1V*T1-k1*V*d1T1;
dy(8,1)=-d2*d1T2-k2*d1V*T2-k2*V*d1T2;
dy(9,1)=k1*d1V*T1+k1*V*d1T1-delta*d1T1s-m1*d1E*T1s-m1*E*d1T1s;
dy(10,1)=k2*d1V*T2+k2*V*d1T2-delta*d1T2s-m2*d1E*T2s-m2*E*d1T2s;
dy(11,1)=NT*delta*(d1T1s+d1T2s)-c*d1V-rho1*k1*d1T1*V-rho1*k1*T1*d1V-rho2*k2*d1T2*V-rho2*k2*T2*d1V;
dy(12,1)=bE*((d1T1+d1T2)*E/(T1s+T2s+Kb)+(T1s+T2s)*d1E/(T1s+T2s+Kb)-(T1s+T2s)*E*(d1T1s+d1T2s)/((T1s+T2s+Kb)^2))-dE*((d1T1+d1T2)*E/(T1s+T2s+Kd)+(T1s+T2s)*d1E/(T1s+T2s+Kd)-(T1s+T2s)*E*(d1T1s+d1T2s)/((T1s+T2s+Kd)^2))-deltaE*d1E;

%k2
dy(13,1)=-d1*k2T1-k1*(k2V*T1+V*k2T1);
dy(14,1)=-d2*k2T2-V*T2-k2*(k2V*T2+V*k2T2);
dy(15,1)=k1*(k2V*T1+V*k2T1)-delta*k2T1s-m1*(k2E*T1s+E*k2T1s);
dy(16,1)=V*T2+k2*(k2V*T2+V*k2T2)-delta*k2T2s-m2*(k2E*T2s+E*k2T2s);
dy(17,1)=NT*delta*(k2T1s+k2T2s)-c*k2V-(rho1*k1*k2T1+rho2*T2+rho2*k2*k2T2)*V-(rho1*k1*T1+rho2*k2*T2)*k2V;
dy(18,1)=bE*((k2T1s+k2T2s)*E/(T1s+T2s+Kb)+((T1s+T2s)*k2E/(T1s+T2s+Kb))-((T1s+T2s)*E*(k2T1s+k2T2s)/((T1s+T2s+Kb)^2)))-dE*((k2T1s+k2T2s)*E/(T1s+T2s+Kd)+((T1s+T2s)*k2E/(T1s+T2s+Kd))-((T1s+T2s)*E*(k2T1s+k2T2s)/((T1s+T2s+Kd)^2)))-deltaE*k2E;

%bE
dy(19,1)=-d1*bET1-k1*(bEV*T1+V*bET1);
dy(20,1)=-d2*bET2-k2*(bEV*T2+V*bET2);
dy(21,1)=k1*(bEV*T1+V*bET1)-delta*bET1s-m1*(bEE*T1s+E*bET1s);
dy(22,1)=k2*(bEV*T2+V*bET2)-delta*bET2s-m2*(bEE*T2s+E*bET2s);
dy(23,1)=NT*delta*(bET1s+bET2s)-c*bEV-(rho1*k1*bET1+rho2*k2*bET2)*V-(rho1*k1*T1+rho2*k2*T2)*bEV;
dy(24,1)=((T1s+T2s)*E/(T1s+T2s+Kb))+bE*(((bET1s+bET2s)*E/(T1s+T2s+Kb))+((T1s+T2s)*bEE/(T1s+T2s+Kb))-((T1s+T2s)*E*(bET1s+bET2s)/((T1s+T2s+Kb)^2)))-dE*(((bET1s+bET2s)*E/(T1s+T2s+Kd))+((T1s+T2s)*bEE/(T1s+T2s+Kd))-((T1s+T2s)*E*(bET1s+bET2s)/((T1s+T2s+Kd)^2)))-deltaE*bEE;

%Original system
dy(25,1) = lam1 - d1*T1 - (1 - epsilon)*k1*V*T1;
dy(26,1) = lam2 - d2*T2 - (1 - f*epsilon)*k2*V*T2;
dy(27,1) = (1 - epsilon)*k1*V*T1 - delta*T1s - m1*E*T1s;
dy(28,1) = (1 - f*epsilon)*k2*V*T2 - delta*T2s - m2*E*T2s;
dy(29,1) = NT*delta*(T1s+T2s) - c*V - ((1 - epsilon)*rho1*k1*T1 + (1 - f*epsilon)*rho2*k2*T2)*V;
dy(30,1) = lamE + (bE*(T1s+T2s)./(T1s+T2s+Kb))*E - (dE*(T1s+T2s)./(T1s+T2s+Kd))*E - deltaE*E;
