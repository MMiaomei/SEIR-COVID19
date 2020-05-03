function ydot = himmelode(t,y,k)

Pop = 6.6*10^7;    % Region actual total population
Iq = y(1);  E = y(2); A = y(3);
I = y(4); R1 = y(5); Eq = y(6);R2 = y(7); Aq = y(8);

eta = k(1); theta = k(2); rho = k(3); beta = k(4);
phi = k(5); epsilon = k(6); alpha = k(7);
gammaA = k(8); gamma = k(9); gammaq = k(10); 
Na = k(11); mu = k(12); zeta=k(19);

ydot = [alpha*eta*Eq+theta*I-gammaq*Iq; 
     phi*(1-rho)*(epsilon*E+I+beta*A)*(Na*Pop-Iq-E-A-I-R1-Eq-R2-Aq)-alpha*E;...
     alpha*(1-eta)*E-mu*A-gammaA*A;...
     alpha*eta*E-gamma*I-theta*I;...
     gammaq*Iq+zeta*Aq;...
     rho*phi*(epsilon*E+I+beta*A)*(Na*Pop-Iq-E-A-I-R1-Eq-R2-Aq)-alpha*Eq;...
     gammaA*A+gamma*I;...
     alpha*(1-eta)*Eq+mu*A-zeta*Aq; 
     ];
 
