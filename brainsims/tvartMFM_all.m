clear; clc
%% TVART folders
p = genpath('/Users/rosun/Documents/MATLAB/BRAINS/tvart-master/src');
addpath(p)
p = genpath('/Users/rosun/Documents/MATLAB/BRAINS/tvart-master/unlocbox');
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



% delays = 10:10:30;
delays = [15, 21, 26, 31];
%% run TVART

Rs = 9:9:90; % up to N

Ms = [100, 200, 1000, 1e4]; % window size
% Rs = [1 2 4 8 16];
% % % % Rs = 9:27:90;
disp(['windows are between ' num2str(Ms(1)*dt*1e3) ' and '  num2str(Ms(end)*dt*1e3) ' ms long'])


rmses =  zeros(length(Ms),length(Rs));
r2s = rmses;


for km = 1:length(Ms)
    M = Ms(km);
    for kr = 1:length(Rs)
        for kd = 1:length(delays)
            delay = delays(kd);
            R = Rs(kr);
            tic
            [lambda, U1, U2, U3, cost, ~,~, rmse] = ...
                TVART_alt_min_delay(X, M, R, delay, ...
                'center', 1, ... % for affine = 1
                'eta', eta, ...
                'beta', beta, ...
                'regularization', regularization, ...
                'verbosity', -1,...
                'max_iter', max_iter);
            toc


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
            rmses(km,kr,kd) = sqrt(norm(Y-Yhat,'fro')^2/(m1*M*T)); % matches!
            if abs(rmses(km,kr,kd)-rmse(end))>0.01*rmse(end) % tolerate up to 1% discrepancy
                error('you done messed up')
            else
                disp('rmses match!')
            end

            r2s(km,kr,kd) = 1-(norm(Y-Yhat,'fro')^2)/(norm(Y-mean(Y,2),'fro')^2);

            disp([rmses(km,kr,kd) r2s(km,kr,kd)])
            figure(2); clf
            plot(rmse); drawnow; pause(.01)
        end
    end
end
save([fname '_errors_vs_M_R_v02'])
%%
figure(1); clf; hold on
colores = parula(length(Rs));
for kd = 1:length(delays)
    subplot(2,4,kd); hold on 
    leyenda = [' '' '];
    for kr = 1:length(Rs)
        plot(Ms,r2s(:,kr,kd),'.-','MarkerSize',22,'color',colores(kr,:))
        leyenda = [leyenda 'R = ' num2str(Rs(kr)) ''','''];
    end
    set(gca,'XScale','log')
end
xlabel('window size (frames)'); ylabel('RMSE')
eval(['legend('  leyenda(1:end-2) ')'])


%% 
figure(2); clf; hold on
leyenda = [' '' '];
for kr = 1:length(Rs)
    plot(Ms,r2s(:,kr),'.-','MarkerSize',22)
    leyenda = [leyenda 'R = ' num2str(Rs(kr)) ''','''];
end
xlabel('window size (frames)'); ylabel('$R^2$')
eval(['legend('  leyenda(1:end-2) ')'])
set(gca,'XScale','log')

C = readmatrix('../structural_connectivity_matrices/S025.csv');
C = C'; % need transpose --> The elements in the first row of the connectivity matrix are proportional to the number of streamlines originating in A and entering the corresponding ROIs.
figure(3);
subplot(1,2,1); imagesc(C)
subplot(1,2,2); imagesc(std(A,[],3))

% tm = (0.5:1:size(U3,1))*M*dt;
% [~,s,v] = svd((U3-mean(U3))','econ');
% figure(5); clf; hold on
% yyaxis left
% plot(t,fpidx)
% yyaxis right
% plot(tm,v(:,1))
% for k = 1:R
%     plot(U3(:,k))
%     pause
% end
%
% [u,s,v] = svd(U3','econ');
% figure(6); clf
% [~, lamb] = myspec_cluster(U3,20);
% % ids = kmedoids(v(:,1:3),3);
% plot(abs(lamb))
%
% % ids = spectralcluster(U3,4);
% nclusts  = sum(abs(lamb)<1e-8) ;
% [idu, ~] = myspec_cluster(U3,nclusts);
% % ids = kmedoids(v(:,1:3),3);
% % s = diag(s);
% % v = u(:,1:npcs)*s(1:npcs,1:npcs)*(v(:,1:npcs)');
% figure(7); clf
% scatter3(v(:,1),v(:,2),v(:,3),100,idu,'filled');
% %
% figure(8); clf; hold on
% yyaxis left
% plot(tm,idu,'-ok','LineWidth',2)
% yyaxis right
% plot(dt*(1:T),fpidx)
%
% %%
% % figure(9);
% % for k = 1:N
% %     clf; hold on
% %     yyaxis left
% %     plot(dt*(1:T),fpidx)
% %     yyaxis right
% %     plot(tm,U3(:,k));
% %     pause
% % end
% %%
%
%
