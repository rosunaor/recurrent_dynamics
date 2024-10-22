function R = wongdecoFiringR(s,C,w,J,Io,a,b,d,G)
x = w*J*s + J*G*C*s +Io;
R = (a*x-b)./(1-exp(-d*(a*x-b)));
end