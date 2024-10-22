clear; clc
%% parameters
C = (ones(2)-eye(2)); % symmetric connectivity
n = size(C,1); % number of neural masses

% % % % % % %%  Mean Field Model (MFM)
% % % % % % % changing G gives 1 - 2 stable equilibria
% % % % % % w = 0.9; % one fixed point in isolation
% % % % % % J  = 0.2609;
% % % % % % Io = 0.3; % this changed
% % % % % % gama = 0.641;
% % % % % % a = 270;
% % % % % % b = 108;
% % % % % % d = 0.154;
% % % % % % taus = 0.1;

%% toti's params
% changing G gives 1 - 4 stable equilibria
w = 1.07;
J  = 0.29;
Io = 0.3;
gama = 0.641;
a = 270;
b = 108;
d = 0.154;
taus = 0.1;



figure(1); clf; hold on
s = linspace(0,0.7);
plot(s,wongdecoMFM(s,0,w,J,Io,gama,a,b,d,taus,0))
plot(s,0*s,'k--','LineWidth',1.5)
box on
xlim([0 s(end)])
xlabel('S'), ylabel('dS/dt')

%% change connection strength
% % G = 0.32; %  2 stable fixed points for MFM parameters
G = 0.03; %  4 stable fixed points for eMFM parameters

%% plot nullclines
figure(2); clf; hold on

fimplicit(@(s1,s2)  wongs(s1,s2,C(1,2),w,J,Io,gama,a,b,d,taus,G),[-.1 1],'k')
fimplicit(@(s1,s2)  wongs(s2,s1,C(2,1),w,J,Io,gama,a,b,d,taus,G),[-.1 1],'r')
[s1,s2] = meshgrid(linspace(-.01,0.7,20)); s1 = s1(:)'; s2 = s2(:)';
s = wongdecoMFM([s1;s2],C,w,J,Io,gama,a,b,d,taus,G);
quiver(s1,s2,s(1,:),s(2,:),2,'color',.1*[1 1 1])
axis equal
pause(.1)

%%
for sigs = 0.5:0.1:1.0
    sigi = sigs; %% noise level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T = 120; % sim time in sec
    % T = 600; % sim time in sec
    %%
    tic
    h = 0.0001; % time step in seconds
    % dt = h;
    t = 0:h:T;
    X = zeros(n,length(t));
    
    xo = 0.35*(1+0.25*randn(n,1)); disp(xo)
    X(:,1) = xo;
    
    
    for k = 2:length(t)
        xo = X(:,k-1);
        z = sqrt(h)*randn(n,1);
        x1 =  xo + sigi*z + h*wongdecoMFM(xo,C,w,J,Io,gama,a,b,d,taus,G);
        X(:,k) = xo + sigi*z + 0.5*h*(wongdecoMFM(xo,C,w,J,Io,gama,a,b,d,taus,G)+wongdecoMFM(x1,C,w,J,Io,gama,a,b,d,taus,G));
    end
    toc
    %%
    figure(2);
    plot(X(1,:),X(2,:),'LineWidth',1,'color',[0 .3 .5 0.1])
    xlabel('$S_1$'); ylabel('$S_2$')
    axis equal
    
    figure(3); clf; hold on
    plot(t,X(1,:),t,X(2,:))
    
    %%
    figure(4); clf; hold on
    hi = histogram2(X(1,:),X(2,:),'DisplayStyle','tile','ShowEmptyBins','on');
    set(hi, 'EdgeColor', 'none');
    fimplicit(@(s1,s2)  wongs(s1,s2,C(1,2),w,J,Io,gama,a,b,d,taus,G),[-.1 1],'w--')
    fimplicit(@(s1,s2)  wongs(s2,s1,C(2,1),w,J,Io,gama,a,b,d,taus,G),[-.1 1],'w--')
    title(['G = ' num2str(G)])
    %% downsample
    t = t(1:20:end); dt = t(2)-t(1)
    X = X(:,1:20:end);
    figure(3);
    plot(t,X(1,:),t,X(2,:))
    
    %% save
    fname = 'deco_totiMFM_2pops_4eq_v';
    k = 1;
    filename = [fname num2str(k)];
    while isfile([filename '.mat'])
         k = k+1; 
         filename = [fname num2str(k)];
    end
    
    disp(filename)
    save(filename,'X','C','w','J','Io','gama','a','b','d','taus','G','h','dt','sigi','T')
    pause(.1)
end

%%
function dsdt = wongs(s1,s2,C,w,J,Io,gama,a,b,d,taus,G)
x = w*J*s1 + J*G*C*s2 +Io;
R = (a*x-b)./(1-exp(-d*(a*x-b)));
dsdt = -(1/taus)*s1 + gama*(1-s1).*R;

end
