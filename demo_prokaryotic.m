% EMSC-NLTHG: single-dataset example for Pan, Li and Che (2024).
% Run demo_prokaryotic from this folder. Requires Statistics and Machine
% Learning Toolbox. Original ER_pro.m update order and parameters retained.
clear
clc
root_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(root_dir,'core'), ...
    fullfile(root_dir,'third_party','clustering'), ...
    fullfile(root_dir,'third_party','gspbox'));
assert(exist('kmeans','file') == 2 && exist('knnsearch','file') == 2, ...
    'Statistics and Machine Learning Toolbox is required.');
rng(0,'twister');
load(fullfile(root_dir,'data','prokaryotic.mat'));
fprintf('Dataset: prokaryotic, %d samples, %d classes, 3 views.\n', ...
    numel(truth), numel(unique(truth)));
gt  = truth;
cls_num = length(unique(gt));
X{1}=text';
X{2}=proteome_comp';
X{3}=gene_repert';
 for v=1:3
    [X{v}]=NormalizeData(X{v});
end
K = length(X); N = size(X{1},2); %sample number

for k=1:K
    Z{k} = zeros(N,N); %Z{2} = zeros(N,N);
    W{k} = zeros(N,N);
    G{k} = zeros(N,N);
    E{k} = zeros(size(X{k},1),N); %E{2} = zeros(size(X{k},1),N);
    FI{k} = zeros(size(X{k},1),N); %Y{2} = zeros(size(X{k},1),N);
end

w = zeros(N*N*K,1);
g = zeros(N*N*K,1);

sX = [N, N, K];

Isconverg = 0;epson = 1e-7;
lambda1 = 0.1;
lambda2 = 2;
alfa=100;
theta=10;
iter = 0;
mu = 10e-5; max_mu = 10e10; pho_mu = 2;
tic;

start =1;

while(Isconverg == 0)
    fprintf('----processing iter %d--------\n', iter+1);
    for k=1:K
        if start==1
          Weight{k} = constructW_PKN((abs(Z{k})+abs(Z{k}'))./2, 3);
          Diag_tmp = diag(sum(Weight{k}));
          L{k} = Diag_tmp - Weight{k};
        else
          P =  (abs(Z{k})+abs(Z{k}'))./2;
          param.k = 3;
          HG = gsp_nn_hypergraph(P', param);
          L{k} = HG.L;
        end
    end
    start = 0;

    for k=1:K

        tmp = X{k}'*(X{k}-E{k}+FI{k}/mu)+G{k}-W{k}/mu;
        Z{k}=inv(lambda2/mu*(L{k}+L{k}')+X{k}'*X{k}+eye(N,N))*tmp;

        F = [X{1}-X{1}*Z{1}+FI{1}/mu;X{2}-X{2}*Z{2}+FI{2}/mu;X{3}-X{3}*Z{3}+FI{3}/mu];
        [Econcat] = solve_log1(F,lambda1/mu);
        E{1} = Econcat(1:size(X{1},1),:);
        E{2} = Econcat(size(X{1},1)+1:size(X{1},1)+size(X{2},1),:);
        E{3} = Econcat(size(X{1},1)+size(X{2},1)+1:end,:);
        FI{k} = FI{k} + mu*(X{k}-X{k}*Z{k}-E{k});
    end

    Z_tensor = cat(3, Z{:,:});
    W_tensor = cat(3, W{:,:});
    z = Z_tensor(:);
    w = W_tensor(:);

      [g, objV] = wshrinkObj_Geman(z + 1/mu*w,1/mu,sX,0,3,alfa,theta)   ;
    G_tensor = reshape(g, sX);

    w = w + mu*(z - g);

    history.objval(iter+1)   =  objV;

    max_norm_Z_G=0;
    max_norm_Z=0;
    Isconverg = 1;
    for k=1:K
        if (norm(X{k}-X{k}*Z{k}-E{k},inf)>epson)
            history.norm_Z = norm(X{k}-X{k}*Z{k}-E{k},inf);
            fprintf('    norm_Z %7.10f    ', history.norm_Z);
            max_norm_Z=max(max_norm_Z,history.norm_Z);
            Isconverg = 0;
        end

        G{k} = G_tensor(:,:,k);
        W_tensor = reshape(w, sX);
        W{k} = W_tensor(:,:,k);
        if (norm(Z{k}-G{k},inf)>epson)
            history.norm_Z_G = norm(Z{k}-G{k},inf);
            fprintf('norm_Z_G %7.10f    \n', history.norm_Z_G);
            max_norm_Z_G=max(max_norm_Z_G,history.norm_Z_G);
            Isconverg = 0;
        end
    end

    % Record actual residuals, including values below the stopping tolerance.
    Recons(iter+1)=max(cellfun(@(x,z,e) norm(x-x*z-e,inf),X,Z,E));
    Match(iter+1)=max(cellfun(@(z,g) norm(z-g,inf),Z,G));

    if (iter>200)
        Isconverg  = 1;
    end
    iter = iter + 1;
    mu = min(mu*pho_mu, max_mu);
end
S = 0;
for k=1:K
    S = S + abs(Z{k})+abs(Z{k}');
end
time_cost=toc;
clear F P

    for i=1:10
        C{i} = SpectralClustering(S,cls_num);
        [Fi(i),Pi(i),Ri(i)] = compute_f(truth,C{i});
        ACCi(i) = Accuracy(C{i},double(truth));
        [A nmii(i) avgenti(i)] = compute_nmi(truth,C{i});
        if (min(truth)==0)
            [ARi(i),RIi(i),MIi(i),HIi(i)]=RandIndex(truth+1,C{i});
        else
            [ARi(i),RIi(i),MIi(i),HIi(i)]=RandIndex(truth,C{i});
        end
    end
ACC(1)=mean(ACCi);ACC(2) = std(ACCi);
nmi(1) = mean(nmii); nmi(2) = std(nmii);
AR(1) = mean(ARi);AR(2) = std(ARi);
F(1) = mean(Fi); F(2) = std(Fi);
P(1) = mean(Pi); P(2) = std(Pi);
R(1) = mean(Ri); R(2) = std(Ri);
avgent(1) = mean(avgenti); avgent(2) = std(avgenti);
result.ACC(1)=mean(ACCi); result.ACC(2) = std(ACCi);
result.nmi(1) = mean(nmii); result.nmi(2) = std(nmii);
result.AR(1) = mean(ARi); result.AR(2) = std(ARi);
result.F(1) = mean(Fi); result.F(2) = std(Fi);
result.P(1) = mean(Pi); result.P(2) = std(Pi);
result.R(1) = mean(Ri); result.R(2) = std(Ri);
result

Result=[ACCi;nmii;ARi;Fi;Pi;Ri;];

assert(all(isfinite(S(:))) && isreal(S), 'Invalid affinity matrix.');
assert(all(isfinite(Result(:))), 'Invalid clustering metrics.');
output_dir = fullfile(root_dir,'results');
if ~exist(output_dir,'dir'), mkdir(output_dir); end
save(fullfile(output_dir,'prokaryotic_result.mat'), ...
    'result','Result','C','S','Recons','Match','iter','time_cost', ...
    'lambda1','lambda2','alfa','theta');
fprintf('ACC %.4f +/- %.4f; NMI %.4f +/- %.4f; AR %.4f +/- %.4f\n', ...
    result.ACC(1),result.ACC(2),result.nmi(1),result.nmi(2), ...
    result.AR(1),result.AR(2));
fprintf('Iterations: %d; final reconstruction residual: %.3e; tensor residual: %.3e\n', ...
    iter,Recons(end),Match(end));
fprintf('Saved results to %s\n',fullfile(output_dir,'prokaryotic_result.mat'));
