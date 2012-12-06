%function [return1 return2] = main(Gene_expression_data,Image_data,survival,Clinical_numeric)
%% Use tri-matrix factorization to uncover latent relationships between
% genes and image data.

%load Breast_cancer_image.mat

% robust sparse non-negative matrix tri-factorization
% 2012. 1. Seung-Jun Kim
% Added a trick to avoid "inadmissible structural zeros"
% Performs column normalization of F & G per iteration




%%set parameters
% outlier matrix O (used Z to avoid confusion with zero)
lambF = 1e-6;  % e.g. 1e-6
lambS = 0;  % e.g. 1e-5
lambG = 0.005;  % e.g. 1e-2;
lambO = 0.6;
kappa = 1e-6;   % Inadmissible structural zero avoidance adjustment (e.g. 1e-6)
kappa_tol = 1e-10;  % Tolerance for identifying a potential structural nonzero (e.g., 1e-10)
epsilon = 1e-10;    % minimum divisor to prevent divide-by-zero (e.g. 1e-10)



%% Prepare data

%remove some genes???


% avoid this kind of colomn or row: sum == 0
index = find(sum(Gene_expression_data,1) == 0);
Gene_expression_data(:,index) = Gene_expression_data(:,index) + eps;
index = find(sum(Image_data,1) == 0);
Image_data(:,index) = Image_data(:,index) + eps;


%MAKE DATA SAME SCALE - normalize by column
Gene_expression_data = zscore(Gene_expression_data);
Image_data = zscore(Image_data);


%make positive
part1 = [max(Gene_expression_data,0) max(-Gene_expression_data,0)];
part2 = [max(Image_data,0) max(-Image_data,0)];


%SORT BY survival/metastasis status
survival = survival/12; %survival from month to year scale
% correct labels 
labels = Clinical_numeric(:,6); %metastasis status
labels(labels==1)=-1; %metastasis
labels(labels==0)=1; %non-metastasis

%put everything in the right order
[Y I] = sort(labels,'descend'); %sort labels by class - I is indices
labels_sorted = labels(I);
W0 = [labels_sorted > 0 labels_sorted < 0];
survival = survival(I,:);
%part1 = part1(I,:); CHANGE BACK
%part2 = part2(I,:);


%combine into data matrix
X = [part1'; part2'];
G = W0;
% 
% 
% % 
%  X(X<3)=0;
%  return1 = X;
%   return2 = W0; 
%   return
% 
%  
 
 
%  
% 
% %make all rows of X nonzero
% index = find(sum(X,2) == 0);
% X(index,:) = X(index,:) + eps;
% 
% %make all columns nonzero
% G = Clinical_numeric(:,[5 6 7 8]); %this is clinical data for survival
% index = find(sum(G,1) == 0);
% G(:,index) = G(:,index) + eps;




[m n] = size(X);
k1 = 10; %columns of F
k2 = 2; %rows of G - want to be the same as number of labels



%% apply algorithm - need to pass in X
%TODO - Include G-G0 to updates (figure out how)

obj = 1e9;
obj_old = 1e10;
F = rand(m,k1);
S = rand(k1,k2);
%G = rand(n,k2);
Z = zeros(m,n);
obj_traj = [];  % collect sequence of objective values
while (obj_old - obj > 1e-6)
    disp('hi')
    upd_factor = ((X-Z)*G*S')./(F*S*G'*G*S' + lambF*sum(sum(F)) + epsilon);
    Q = kappa*(upd_factor > 1 & F < kappa_tol);
    if sum(sum(Q)) > 0 disp(sprintf('inadmissible zeros in F (iter = %d)',length(obj_traj))), end
    F = (F + Q).*upd_factor;
    [F,d1] = norm_cols(F,2); S = diag(d1)*S;
    
    upd_factor = (F'*(X-Z)*G)./(F'*F*S*G'*G + lambS*sum(sum(S)) + epsilon);
    Q = kappa*(upd_factor > 1 & S < kappa_tol);
    if sum(sum(Q)) > 0 disp(sprintf('inadmissible zeros in S (iter = %d)',length(obj_traj))), end
    S = (S + Q).*upd_factor;
    
%    upd_factor = ((X-Z)'*F*S)./(G*S'*F'*F*S + lambG*sum(sum(G)) + epsilon);
%    Q = kappa*(upd_factor > 1 & G < kappa_tol);
%    if sum(sum(Q)) > 0 disp(sprintf('inadmissible zeros in G (iter = %d)',length(obj_traj))), end
%    G = (G + Q).*upd_factor;
%    [G,d2] = norm_cols(G,2); S = S*diag(d2);
    
    T = X - F*S*G';
    Z = min(X, sign(T).*max(abs(T) - lambO,0));
    
    obj_old = obj;
    obj = 0.5*(norm(X - F*S*G' - Z,'fro')^2 + lambF*sum(sum(F))^2 + ...
        lambS*sum(sum(S))^2 + lambG*sum(sum(G))^2) + lambO*sum(sum(abs(Z)));
    obj_traj = [obj_traj obj];
    plot(obj_traj)
    pause(.1)
end
disp(sprintf('Total # of iteration = %d',length(obj_traj)));

% normalize
%[F,d1] = norm_cols(F,1);
%[G,d2] = norm_cols(G,1);
%S = diag(d1)*S*diag(d2);