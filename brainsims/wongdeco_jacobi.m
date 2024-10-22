function Jaco = wongdeco_jacobi(s,C,varargin)
if nargin < 2.2 
    w = 0.95; % one fixed point
    J  = 0.2609;
    Io = 0.32; 
    gama = 0.641;
    a = 270;
    b = 108;
    d = 0.154;
    taus = 0.1;  
    G = 0.55;
else
    w = varargin{1};
    J  = varargin{2};
    Io = varargin{3};
    gama = varargin{4};
    a = varargin{5};
    b = varargin{6};
    d = varargin{7};
    taus = varargin{8};
    G = varargin{9};
end 
n = size(C,1); 
x = w*J*s + J*G*C*s +Io;
R = (a*x-b)./(1-exp(-d*(a*x-b)));

dxds = w*J*eye(n) + J*G*C; 
so = 1-exp(-d*(a*x-b)); 
dRds = (a*dxds.*so - d*a*dxds.*(a*x-b).*exp(-d*(a*x-b)))./(so.^2); 

Jaco = -(1/taus)*eye(n) - gama*eye(n).*R + gama*(1-s).*dRds; 
% % % % % Jaco
% % % % for k = 1:n
% % % %     for m = 1:n 
% % % %         dRds = (a*dxds(k,m)*so(k) - d*a*dxds(k,m)*(a*x(k)-b)*exp(-d*(a*x(k)-b)))/(so(k).^2); 
% % % %         Jaco(k,m) = -(1/taus)*eq(k,m) + (1-eq(k,m))*gama*R(k) +(1-s(k))*gama*dRds; 
% % % %     end
% % % % end 


end

















