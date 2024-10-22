clear; clc
%% parameters
% Connectivity
C = readmatrix('..\structural_connectivity_matrices\S025.csv');
C = C'; % need transpose --> The elements in the first row of the connectivity matrix are proportional to the number of streamlines originating in A and entering the corresponding ROIs.
n = size(C,1); % number of neural masses
%% set params
w = 0.95; % one fixed point in isolation
J  = 0.2609;
Io = 0.32;
gama = 0.641;
a = 270;
b = 108;
d = 0.154;
taus = 0.1;
% % G = 0.08;   
G = 0.55; 
%%

tic
h = 0.0001; % time step in seconds
dt = h;
t = 0:h:601;
X = zeros(n,length(t));
X(:,1) = rand(n,1);
sigi = 0.25; %% noise level %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 2:length(t)
    xo = X(:,k-1);
    z = sqrt(h)*randn(n,1);
    x1 =  xo + sigi*z + h*wongdecoMFM(xo,C,w,J,Io,gama,a,b,d,taus,G);
    X(:,k) = xo + sigi*z + 0.5*h*(wongdecoMFM(xo,C,w,J,Io,gama,a,b,d,taus,G)+wongdecoMFM(x1,C,w,J,Io,gama,a,b,d,taus,G));
end

S = X; clear X
toc
%%
% x = w*J*S + J*G*C*S +Io;
% R = (a*x-b)./(1-exp(-d*(a*x-b)));
% x = x-Io;
%% get rid of transient start 
ind = t>=10; 
S = S(:,ind); 
t = t(ind); 

%% find fixed point at each time
fname = ['..\S025_FPs_G_' num2str(G*1000) '_w_' num2str(w*100)];  
load(fname)
FPs = FPs';
figure(7); clf; hold on 
for k = 1:size(FPs,2)
    plot(FPs(:,k)+k-1)
end


%%
fpd   = zeros(size(S,2),1);
fpidx = zeros(size(S,2),1);
tic
for k = 1:size(S,2)
    [dist,idx] = min(mean((FPs-S(:,k)).^2));
    
    fpd(k) = dist;
    fpidx(k) = idx;
    
end
toc
fpd = sqrt(fpd); 
%% 
figure(11); clf; hold on 
[~,edges] = histcounts(fpd,0.2*linspace(0,1,30)); 
for k = 1:size(FPs,2) 
    Nd = histcounts(fpd(fpidx==k),edges,'normalization','pdf'); 
    plot(0.5*(edges(1:end-1)+edges(2:end)),Nd,'LineWidth',2)
end 
set(gca,'YScale','lin')
%%
figure(12); clf 
plot(t,fpidx,'LineWidth',2)
figure(13); clf 
histogram(fpidx)
%% 
figure(14); clf; hold on 
[~,ind]  = min(abs(t-20.2001)); 
plot(S(:,ind),'k','LineWidth',1.5)
didi = sqrt(mean((S(:,ind)-FPs).^2)); 
[didi, didik] = mink(didi,2); 
plot(FPs(:,didik(1)),'LineWidth',1.5)
plot(FPs(:,didik(2)),'--','LineWidth',1.5)
%%
%% downsample 
S = S(:,1:50:end); 
t = t(1:50:end); 
fpidx = fpidx(1:50:end); 

dt = dt*50; disp(1/dt)

ind = t>1; 
figure(1); clf; hold on
subplot(2,1,1)
imagesc(t(ind),1:n,S(:,ind)); caxis([prctile(S(:),0.01) prctile(S(:),99)])
xlabel('t (s)'); ylabel('neural pop')
h = colorbar;
h.Label.String = 'S';
subplot(2,1,2)
plot(S(:,end))
%% save data set 


figure(12); hold on 
plot(t,fpidx,'--','LineWidth',2)
ylim([0 size(FPs,2)+1])

%% 
fname = ['sim_x_fpidx_S025_G_' num2str(1000*G) '_w_'  num2str(100*w) '_s_'  num2str(100*sigi)]; 
save([fname '_v02'],'S','fpidx','dt','sigi','G') 

%%
% ind = t>1;
% figure(1); clf; hold on
% subplot(2,1,1)
% imagesc(t(ind),1:n,R(:,ind)); caxis([prctile(R(:),0.01) prctile(R(:),99)])
% xlabel('t (s)'); ylabel('neural pop')
% h = colorbar;
% h.Label.String = 'R';
% subplot(2,1,2)
% plot(R(:,end))


% figure(2); clf; hold on
% subplot(2,1,1)
% imagesc(t(ind),1:n,S(:,ind)); caxis([prctile(S(:),0.01) prctile(S(:),99)])
% xlabel('t (s)'); ylabel('neural pop')
% h = colorbar;
% h.Label.String = 'S';
% subplot(2,1,2)
% plot(S(:,end))

% figure(3); clf; hold on
% subplot(2,1,1)
% imagesc(t(ind),1:n,x(:,ind)); caxis([prctile(x(:),0.01) prctile(x(:),99)])
% xlabel('t (s)'); ylabel('neural pop')
% h = colorbar;
% h.Label.String = 'syn input';
% subplot(2,1,2)
% plot(x(:,end))



%%
% figure(4); clf;
% 
% y = S(22,:);
% subplot(2,1,1)
% plot(t,y);
% subplot(2,1,2)
% [om,ps] = ezfft(t,y);
% loglog(om/(2*pi),ps)


%% find average switching time 

% d = diff(fpidx(1e5:end)); d = diff(find(d>0)); 
% 
% [~,mps,bin] = histcounts(d,30); 
% mps = 0.1*0.5*(mps(1:end-1)+mps(2:end)); % midpoints of bins in ms 
% totime = 0*mps; 
% for k = 1:length(mps) 
%     totime(k) = sum(d(bin==k)); 
% end 
% totime = totime/sum(d); 
% figure(14); clf
% bar(mps,totime)
% figure(15); clf
% plot(mps,cumsum(totime))


