Source code for the following paper:

    Xiaoqian Wang, Feiping Nie, Heng Huang.
    Structured Doubly Stochastic Matrix for Graph Based Clustering. 
    The 22nd ACM SIGKDD Conference on Knowledge Discovery and Data Mining (KDD 2016).

This function optimizes the following problem:

    min_{S=S', S>=0, S1=1, Tr(S)=k, rank(L_S)=n-k} ||S-A0||_F^2 + r*||S||_F^2

Format of input:

    A0: the input affinity matrix
    k: number of clusters
    r: parameter for ||S||_F^2
    gama: parameter for ||L||*. If gamma=0, non-convex; if gamma>0, convex
    mu: penalty parameter for ALM
    rho: increment step for ALM
    NITER: maximum number of iterations.

Simply run the code in matlab as below:

    [Ind, S, F, obj, conv] = bistochasticKcut(A0, k, gamma, r, mu, rho, NITER);
