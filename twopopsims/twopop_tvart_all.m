clear; clc;
%% TVART folders
p = genpath('C:\Users\rosun\Documents\MATLAB\BRAINS\Tensor_Decomps\tvart-master\src');
addpath(p)
p = genpath('C:\Users\rosun\Documents\MATLAB\BRAINS\Tensor_Decomps\tvart-master\unlocbox');
addpath(p)

%% load data
% % fname = 'deco_MFM_2pops_2eq_v4.mat';
% % fname = 'deco_totiMFM_2pops_4eq_v4.mat';
% % load(fname)


files = dir('deco_MFM_2pops_2eq_v*');
% delays = 1:2:31;
delays = [1,16,31,60,100,150,250,500,750,1000]; % explore up to delays in order of seconds, which is delay of 500 

sigis = zeros(numel(files),length(delays));
rmses = sigis;
r2= sigis;
loglds= sigis;
dmahas= sigis;
posts= sigis;

U1all = struct; 
U2all = struct; 
U3all = struct; 


for kf = 1:numel(files)
    load(files(kf).name)
    %% plot time series
    t = ((1:length(X))-1)*dt;
    figure(1); clf; hold on
    plot(t,X)
    
    
    %% tvart
    %%
    N = size(X,1); disp(N)
    %% TVART parameters
    regularization = 'spline'; % regularization type: 'TV', 'spline,' or 'TVL0'
    max_iter = 200;
    R    = N;
    eta  = 1/2; % Tikhonov regularization ~1/eta
    beta = .5; % temporal smoothing strength
    M = 100; % window size
    disp(['windows are ' num2str(M*dt/taus) ' tau_s long'])
    disp(['windows are ' num2str(M*dt*1000) ' ms long'])
    
    
    
    
    
    
    for kd = 1:length(delays)
        tic
        sigis(kf,kd) = sigi;
        delay = delays(kd);
        
        
        [lambda, U1, U2, U3, cost, ~,~, rmse] = ...
            TVART_alt_min_delay(X, M, R, delay, ...
            'center', 1, ...% for affine
            'eta', eta, ...
            'beta', beta, ...
            'regularization', regularization, ...
            'verbosity', -1,...
            'iter_disp',50,...
            'max_iter', max_iter);
        
        
        %%
        T = floor((size(X,2)-delay) / M);
        N = size(U3, 1);
        m1 = size(U1,1);
        m2 = size(U2, 1);
        
        A = zeros(m1,m2,N);
        eigs = zeros(m1,N);
        normdiffs = zeros(N-1,1);
        Yhat = zeros(size(X,1), M*T);
        for k = 1:N
            A(:,:,k) = U1*diag(U3(k,:))*(U2');
            if k>1
                normdiffs(k-1) = norm(A(:,:,k)-A(:,:,k-1),'fro');
            end
            eigs(:,k) = eig(A(:,1:end-1,k));
            
            % predictions
            x = X(:,(1+(k-1)*M):k*M);
            Yhat(:,(1+(k-1)*M):k*M) = A(:,1:end-1,k)*x+A(:,end,k);
        end
        Y = X(:,(1+delay):M*T+delay);
        rmses(kf,kd) = sqrt(norm(Y-Yhat,'fro')^2/(m1*M*T)); % matches!
        if abs(rmses(kf,kd)-rmse(end))>0.01*rmse(end) % tolerate up to 1% discrepancy
            error('you done messed up')
        else
            disp('rmses match!')
        end
        
        r2(kf,kd) = 1-(norm(Y-Yhat,'fro')^2)/(norm(Y-mean(Y,2),'fro')^2);
        
        if sum(imag(eigs(:))~=0)<.1
            disp('all real eigenvalues')
        else
            disp('some complex eigenvalues')
        end
        
        %%
        figure(7); clf; hold on
        GMModel = fitgmdist(U3,1);
        logld = GMModel.NegativeLogLikelihood;
        GMModel = fitgmdist(U3,2);
        logld = -GMModel.NegativeLogLikelihood+logld;
        loglds(kf,kd) = logld;
        
        idx = cluster(GMModel,U3);
        % posterior probabilities
        pos = posterior(GMModel,U3); pos = max(pos,[],2);
        % mahalanobis distance
        dmaha = mahal(GMModel,U3); dmaha = abs(dmaha(:,1)-dmaha(:,2));
        
        dmahas(kf,kd) = mean(dmaha);
        posts(kf,kd) = mean(pos);
        
        scatter(U3(idx==1,1),U3(idx==1,2),50*pos(idx==1),'filled')
        scatter(U3(idx==2,1),U3(idx==2,2),50*pos(idx==2),'filled')
        
       U3all(kf,kd).U1 = U1;
       U3all(kf,kd).U2 = U2;
       U3all(kf,kd).U3 = U3;
        
        haxis = gca;
        xlim = haxis.XLim;
        ylim = haxis.YLim;
        d = (max([xlim ylim])-min([xlim ylim]))/1000;
        [X1Grid,X2Grid] = meshgrid(xlim(1):d:xlim(2),ylim(1):d:ylim(2));
        
        contour(X1Grid,X2Grid,reshape(pdf(GMModel,[X1Grid(:) X2Grid(:)]),...
            size(X1Grid,1),size(X1Grid,2)),20)
        pause(.1)
        toc
    end
    
end

toc
delays = repmat(delays,size(sigis,1),1); 
fname = 'deco_MFM_2pops_2eq_TVARTclustermetrics';
save(fname,'files','sigis','rmses','r2', 'loglds','dmahas','posts','delays','U3all')

%% 
figure(7); clf; hold on
% subplot(2,2,1)
% surf(sigis,delays,posts)
% xlabel('\sigma'); ylabel('delays'); zlabel('RMSE') 

for k = 1:size(posts,1)
    plot(delays(k,:),posts(k,:))
    pause
end

% %%
% figure(2); clf
% subplot(2,2,1); hold on
% yyaxis left; hold on
% plot(delays,rmses); ylabel('RMSE')
% % plot(delays,0.045+delays*dt*sigi,'--')
% yyaxis right
% plot(delays,r2); ylabel('R^2')
%
% xlabel('delay')
%
% subplot(2,2,2);
% plot(delays,loglds); ylabel('loglikelihood diff')
% xlabel('delay')
%
% subplot(2,2,3);
% plot(delays,dmahas); ylabel('mean mahalanobis diff')
% xlabel('delay')
%
% subplot(2,2,4);
% plot(delays,posts); ylabel('mean posterior probability')
% xlabel('delay')