clear; clc
%% TVART folders
p = genpath('C:\Users\rosun\Documents\MATLAB\BRAINS\Tensor_Decomps\tvart-master\src');
addpath(p)
p = genpath('C:\Users\rosun\Documents\MATLAB\BRAINS\Tensor_Decomps\tvart-master\unlocbox');
addpath(p)
%% load simulation data
tic
fname = 'sim_x_fpidx_S025_G_80_w_95_s_25_v02';
% % % fname = 'sim_x_fpidx_S025_G_550_w_95_s_25_v02';
load(fname)
toc
X = S; clear S

%%
% N = size(x,1);
t = ((1:(length(X)))-1)*dt;
% T = size(X,2);
% % % % % % X = X-mean(X,2);
N = size(X,1);
disp([num2str(N) ' populations'])
T = size(X,2);
disp([num2str(T) ' frames'])

%% TVART parameters
regularization = 'spline'; % 'TV', 'spline', 'TVL0'
max_iter = 40;

eta  = 10; % penalty 1/eta
beta = 0.1; % temporal smoothing strength
delays = 1:10:61;
%% run TVART

% Rs = 9:9:90; % up to N

M = 100; % window size
R = 30;
% % % % Rs = 9:27:90;
disp(['windows are between ' num2str(M(1)*dt*1e3) ' and '  num2str(M(end)*dt*1e3) ' ms long'])


rmses =  zeros(length(delays),1);
r2s = rmses;
Uall = struct();


for kd = 1:length(delays)
    delay = delays(kd);
    tic
    [~, U1, U2, U3, cost, ~,~, rmse] = ...
        TVART_alt_min_delay(X, M, R, delay, ...
        'center', 1, ... % for affine = 1
        'eta', eta, ...
        'beta', beta, ...
        'regularization', regularization, ...
        'verbosity', -1,...
        'max_iter', max_iter);
    
    Uall(kd).U1 = U1;
    Uall(kd).U2 = U2;
    Uall(kd).U3 = U3;

    %%
    T = floor((size(X,2)-delay) / M);
    N = size(U3, 1);
    m1 = size(U1,1);
    m2 = size(U2, 1);

    A = zeros(m1,m2,N);
    Yhat = zeros(size(X,1), M*T);

    for k = 1:N
        A(:,:,k) = U1*diag(U3(k,:))*(U2');
        % predictions
        x = X(:,(1+(k-1)*M):k*M);
        Yhat(:,(1+(k-1)*M):k*M) = A(:,1:end-1,k)*x+A(:,end,k);
    end

    Y = X(:,(1+delay):M*T+delay);
    rmses(kd) = sqrt(norm(Y-Yhat,'fro')^2/(m1*M*T)); % matches!
    if abs(rmses(kd)-rmse(end))>0.01*rmse(end) % tolerate up to 1% discrepancy
        error('you done messed up')
    else
        disp('rmses match!')
    end

    r2s(kd) = 1-(norm(Y-Yhat,'fro')^2)/(norm(Y-mean(Y,2),'fro')^2);

    disp([rmses(kd) r2s(kd)])
    figure(2); clf
    plot(rmse); drawnow; pause(.01)
    toc
end

%%
t0 = t(1); 
tend = t(end); 
clear A eigs fpidx lambda U1 U2 U3 x X y Yhat t  
save([fname '_R_' num2str(R) '_M_' num2str(M) '_TVART_delays'])

















