function dsdt = wongs(s1,s2,C,w,J,Io,gama,a,b,d,taus,G)
x = w*J*s1 + J*G*C*s2 +Io;
R = (a*x-b)./(1-exp(-d*(a*x-b)));
dsdt = -(1/taus)*s1 + gama*(1-s1).*R;

end