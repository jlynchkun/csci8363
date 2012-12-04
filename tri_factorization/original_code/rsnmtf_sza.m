% robust sparse non-negative matrix tri-factorization
% 2012. 1. Seung-Jun Kim
% Added a trick to avoid "inadmissible structural zeros"

clear
rng('default')

m = 10;
n = 15;
k1 = 3;
k2 = 4;
F = zeros(m,k1);
S = zeros(k1,k2);
G = zeros(n,k2);
Z = zeros(m,n);     % outlier matrix O (used Z to avoid confusion with zero)
lambF = 1e-6;  % e.g. 1e-6
lambS = 1e-5;  % e.g. 1e-5
lambG = 500;  % e.g. 1e-2;
lambO = 0.6;
kappa = 1e-6;   % Inadmissible structural zero avoidance adjustment (e.g. 1e-6)
kappa_tol = 1e-10;  % Tolerance for identifying a potential structural nonzero (e.g., 1e-10)
epsilon = 1e-10;    % minimum divisor to prevent divide-by-zero (e.g. 1e-10)

% generate test data
F(1:3, 1) = 1;
F(4:7, 2) = 1; F(4:5, 1) = 1;
F(7:10,3) = 1;
G(1:4, 1) = 1;
G(5:9, 2) = 1;
G(9:12,3) = 1;
G(13:15,4) = 1;
S(1,[1 3]) = 1;
S(2,[2 ]) = 1;              %S(2,[2 3]) = 1;
S(3,[3 4]) = 1;             %S(3,4) = 1;

X = F*S*G' > 0.5;
%X = F*S*G';
Z(1,1) = 1; % outliers
Z(4,10) = 1;

Xtrue = X;
Ftrue = F; 
Strue = S;
Gtrue = G;
Ztrue = Z;

X = X - Z;

% apply algorithm
obj = 1e9;
obj_old = 1e10;

% initialize factors
Z = zeros(m,n);
% randomly
F = rand(m,k1);
S = rand(k1,k2);
G = rand(n,k2);

% use nmf
%[F,dummy] = nmf(X,k1,lambF);
%[G,dummy] = nmf(X',k2,lambG);
%S = rand(k1,k2);

obj_traj = [];  % collect sequence of objective values
while (obj_old - obj > 1e-6)
    upd_factor = (F'*(X-Z)*G)./(F'*F*S*G'*G + lambS*sum(sum(S)) + epsilon);
    Q = kappa*(upd_factor > 1 & S < kappa_tol);
    if sum(sum(Q)) > 0 disp(sprintf('inadmissible zeros in S (iter = %d)',length(obj_traj))), end
    S = (S + Q).*upd_factor;
    
    upd_factor = ((X-Z)*G*S')./(F*S*G'*G*S' + lambF*sum(sum(F)) + epsilon);
    Q = kappa*(upd_factor > 1 & F < kappa_tol);
    if sum(sum(Q)) > 0 disp(sprintf('inadmissible zeros in F (iter = %d)',length(obj_traj))), end
    F = (F + Q).*upd_factor;
    
    upd_factor = ((X-Z)'*F*S)./(G*S'*F'*F*S + lambG*sum(sum(G)) + epsilon);
    Q = kappa*(upd_factor > 1 & G < kappa_tol);
    if sum(sum(Q)) > 0 disp(sprintf('inadmissible zeros in G (iter = %d)',length(obj_traj))), end
    G = (G + Q).*upd_factor;    
    
    T = X - F*S*G';
    Z = min(X, sign(T).*max(abs(T) - lambO,0));
    
    obj_old = obj;
    obj = 0.5*(norm(X - F*S*G' - Z,'fro')^2 + lambF*sum(sum(F))^2 + ...
        lambS*sum(sum(S))^2 + lambG*sum(sum(G))^2) + lambO*sum(sum(abs(Z)));
    obj_traj = [obj_traj obj];
end
disp(sprintf('Total # of iteration = %d',length(obj_traj)));

% normalize
% [F,d1] = norm_cols(F,1);
% [G,d2] = norm_cols(G,1);
% S = diag(d1)*S*diag(d2);
