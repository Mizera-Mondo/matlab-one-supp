function G = laplacianProjection(M)
%LAPLACIANPROJECTION Project the input matrix so as to get a nearest valid Laplacian.
%   Input: M
if ~issymmetric(M)
    M = (M + M')/2;
end

% Initialization
G = M;
K = M;
[N, ~] = size(M);
Phi = zeros(N);
rho = 0.3;
m = mat2vec(M);

% Quadratic constriants
A = [];
b = [];
Aeq = [];
beq = [];

% Gij = Gji <= 0
for i = 1:N
    for j = 1:i - 1

        T = zeros(N);
        T(i, j) = 1;
        A = [A, mat2vec(T)];
        b = [b, 0];

        T(j, i) = -1;
        Aeq = [Aeq, mat2vec(T)];
        beq = [beq, 0];

    end
end

% G1 = 0
for i = 1:N
    T = zeros(N);
    T(i, :) = ones(1, N);
    Aeq = [Aeq, mat2vec(T)];
    beq = [beq, 0];
end

% trace(G) = N
Aeq = [Aeq, mat2vec(eye(N))];
beq = [beq, N];

A = A';
b = b';
Aeq = Aeq';
beq = beq';


% Technical varbs
iter = 1;
maxIter = 100;
isAdmmConverge = false;
isMaxIter = false;

while ~isAdmmConverge && ~isMaxIter
    G_old = G;
    % Non-SPD part w.r.t G, using quadratic programming solver
    g = mat2vec(G);
    k = mat2vec(K);
    u = mat2vec(Phi);
    H = (rho + 1)*eye(length(g));
    f = (u - m - rho*k);
    g = quadprog(H, f, A, b, Aeq, beq, [], [], g, optimoptions(@quadprog, 'Display', 'off'));
    G = vec2squ(g);
    % SPD part w.r.t K, using 1-step SPD projection
    U = G + 1/rho*Phi;
    K = nearestSPD(U);
    
    % Update Phi
    Phi = Phi + rho*(G - K);
    % Convergence Check
    isAdmmConverge = norm(G - G_old, 'fro')/norm(G_old) < 1e-5;
    isMaxIter = iter >= maxIter;
    iter = iter + 1;
end
end
