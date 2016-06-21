% min_{S=S', S>=0, S1=1, Tr(S)=k, rank(L_S)=n-k} ||S-A0||_F^2 + r*||S||_F^2
% Xiaoqian Wang, Feiping Nie, Heng Huang. 
% Structured Doubly Stochastic Matrix for Graph Based Clustering. 
% The 22nd ACM SIGKDD Conference on Knowledge Discovery and Data Mining (KDD 2016).
% Author: Xiaoqian Wang
function [Ind, S, F, obj, conv] = bistochasticKcut(A0, k, gamma, r, mu, rho, NITER)
% Input:
%   A0: the input affinity matrix
%   k: number of clusters
%   r: parameter for ||S||_F^2
%   gama: parameter for ||L||*. If gamma=0, non-convex; if gamma>0, convex
%   mu: penalty parameter for ALM
%   rho: increment step for ALM
%   NITER: maximum number of iterations.
% Output:
%   Ind: clustering indicator
%   S: the learnt block diagonal similarity matrix
%   F: eigenvectors of L_S
%   obj: objective function value
%   conv: a parameter to check if converges, small value means convergence.

%% Initialization
%
n = size(A0, 1); % number of samples
Lambda = zeros(n); % Lagrange multiplier
In = eye(n);
obj = zeros(NITER,1);

%
S = A0;
L = In - S; % graph Laplacian


%% Iteration
for iter = 1:NITER
    
    inmu = 1 / mu;
    
    % update L
    if gamma == 0 % non-convex case
        [UU, SS, VV] = svd(In-S+inmu*Lambda);
        [SS_value, SS_index] = sort(diag(SS), 'descend');

        UU = UU(:, SS_index(1: n-k));
        VV = VV(:, SS_index(1: n-k));
        SS = SS(SS_index(1: n-k), SS_index(1: n-k));
    
        L = UU * SS * VV';
    else % convex case
        L = Proximal_tracep(In-S+inmu*Lambda, gamma*inmu, 1);
    end;
    
    % update S
    V = (A0 + 0.5*mu*(In-L) + 0.5*Lambda) / (r+0.5*mu+1);
    S = solveS(V, k);
    
    % update Lambda
    Lambda = Lambda + mu*(In-S-L);
    
    % update mu for ALM
    mu = min(10^10,rho*mu);
    
    % compute the object function value
    if gamma ~= 0
        [UU, SS, VV] = svd(In-S);
    end;
    obj(iter) = norm(S-A0,'fro')^2 + r*norm(S,'fro')^2 + 0.5*mu*norm(In-S-L+inmu*Lambda,'fro')^2 + sum(gamma*diag(SS));
    
    if mod(iter,10) == 0
        fprintf('%d-th iteration, obj=%f\n', iter, obj(iter));
    end
end;

% Update cluster indicators with graph laplacian L_S
[F, temp, ev]=eig1(L, k, 0);
[labv, tem, Ind] = unique(round(0.1*round(1000*F)),'rows');
S = In - L;

% Check if this function converges
da = mu*(In-S-L);
conv = max(abs(da(:)));

end