function dsdt = wongdecoMFM(s,C,varargin)
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

x = w*J*s + J*G*C*s +Io;
R = (a*x-b)./(1-exp(-d*(a*x-b)));

dsdt = -(1/taus)*s + gama*(1-s).*R;
end