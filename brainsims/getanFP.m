function [FP,fval] = getanFP(xo, C, G)

%% set params
w = 0.95; % one fixed point in isolation
J  = 0.2609;
Io = 0.32;
gama = 0.641;
a = 270;
b = 108;
d = 0.154;
taus = 0.1;

%%

Decomfm = @(x) wongdecoMFM(x,C,w,J,Io,gama,a,b,d,taus,G);

h = 0.001; % time step in seconds
tol = h/100; % convergence tolerance for forward simulation

options = optimoptions('fsolve','FiniteDifferenceType','central','OptimalityTolerance',1e-10,'Display','none');




conved = norm(wongdecoMFM(xo,C,w,J,Io,gama,a,b,d,taus,G));
nops = 0;
X = xo; 
while conved > tol
    xo = X;
    x1 =  xo +  h*wongdecoMFM(xo,C,w,J,Io,gama,a,b,d,taus,G);
    X = xo + 0.5*h*(wongdecoMFM(xo,C,w,J,Io,gama,a,b,d,taus,G)+wongdecoMFM(x1,C,w,J,Io,gama,a,b,d,taus,G));

    %                 conved = max(abs(Decomfm(X)));
    conved = norm(Decomfm(X));
    nops = nops+1;
end
%             disp (nops)
%             disp(conved)

[FP,fval] = fsolve(Decomfm,X,options);











end