clear; clc
%% parameters
% Connectivity
C = readmatrix('../structural_connectivity_matrices/S025.csv');
C = C'; % need transpose --> The elements in the first row of the connectivity matrix are proportional to the number of streamlines originating in A and entering the corresponding ROIs.

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
% % % % load('sim_x_fpidx_S025_G_550_w_95_s_25_v02.mat')
load('sim_x_fpidx_S025_G_80_w_95_s_25_v02.mat')
stepp = 2; 
tic
FPs = zeros(size(S,1), floor(size(S,2)/stepp));

nf = 0; 
for n = 1:stepp:size(S,2)
    nf = nf+1; 
    xo = S(:,n);
    [FP,fval] = getanFP(xo, C, G);
    FPs(:,nf) = FP;
    if mod(nf,400) == 0
        disp(toc)
    end
end

toc

% %%
% figure(1); clf
% imagesc(FPs)

%% find unique fixed points
tic
nfxps = 1;
[~,~,sumd,~] = kmeans(FPs',nfxps,'distance','cityblock');
while max(sumd)>1e-3 && nfxps<1000
    nfxps = nfxps+1;
    [~,~,sumd,~] = kmeans(FPs',nfxps,'distance','cityblock');
end

[idx,FPsuni,sumd,D] = kmeans(FPs',nfxps,'distance','cityblock');
toc 
FPsu = FPsuni'; 
figure(2); clf 
plot(FPsu,'.-')

save('sim_x_fpidx_S025_G_80_w_95_s_25_v02_uniq_FPs', 'FPsu')