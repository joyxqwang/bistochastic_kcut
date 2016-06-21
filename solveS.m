% min_{S>=0, S=S', S*1=1, Tr(S)=k}  ||S-V||_F^2
function S = solveS(V, k)

%% Initialization
n = size(V, 1);
S = V;

%% Iteration
for i = 1 : 10
    
    % Lemma 3
    S = max(S, 0);
    [diagS ft] = EProjSimplex_new(diag(S), k);
    for jj = 1 : n
        S(jj, jj) = diagS(jj);
    end;
    
    % Lemma 2
    M = (S + S') / 2;
    %S = M + (n+one'*M*one)*one*one'/(n^2) - 1/n*M*one*one' - 1/n*one*one'*M;
    S = M + (n+sum(sum(M)))/(n^2)*ones(n) - 1/n*repmat(sum(M),n,1) - 1/n*repmat(sum(M,2),1,n);
end;

        