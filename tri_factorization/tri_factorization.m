% X1 and X2 should be normalized
function [F S G] = tri_factorization(X1, X2, labels, parameters)
%function [return1 return2] = main(Gene_expression_data,Image_data,survival,Clinical_numeric)
%% Use tri-matrix factorization to uncover latent relationships between
% genes and image data.

% robust sparse non-negative matrix tri-factorization
% 2012. 1. Seung-Jun Kim
% Added a trick to avoid "inadmissible structural zeros"
% Performs column normalization of F & G per iteration

%% STRATEGY
%1. F is labels only (normalized columns? Try different labels?)
%2. Use k-means clustering to initialize G, then remember to transpose...
%3. Calculate S from initial G and F
%4. Update G, S but NOT F

%Minimize ||X - F*S*G^T|| + lambda_G * ||G||_1 + lambda_S * ||S||_1
%where F is binary matrix of labels and is fixed
%G is initialized to k-means clusters of features, but then updated
%S is randomly initialized then updated







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
%{
%remove some genes???
if false
load ../data/Breast_cancer_image.mat
% avoid this kind of colomn or row: sum == 0
index = find(sum(Gene_expression_data,1) == 0);
Gene_expression_data(:,index) = Gene_expression_data(:,index) + eps;
index = find(sum(Image_data,1) == 0);
Image_data(:,index) = Image_data(:,index) + eps;
%MAKE DATA SAME SCALE - normalize by column
Gene_expression_data = zscore(Gene_expression_data);
Image_data = zscore(Image_data);
%}
% Image_data -> X2
% Gene_expression_data -> X1
%make positive
nn_X1 = [max(X1.data,0) max(-X1.data,0)];
nn_X2 = [max(X2.data,0) max(-X2.data,0)];


%SORT BY survival/metastasis status
%survival = survival/12; %survival from month to year scale
% correct labels 
%labels = Clinical_numeric(:,6); %metastasis status
%labels(labels==1)=-1; %metastasis
%labels(labels==0)=1; %non-metastasis

%put everything in the right order
[Y I] = sort(labels.numeric,'descend'); %sort labels by class - I is indices
labels_sorted = labels.numeric(I);
W0 = [labels_sorted > 0 labels_sorted < 0];
%survival = survival(I,:);
nn_X1 = nn_X1(I,:);
if isempty(nn_X2)
  % do not sort
else
  nn_X2 = nn_X2(I,:);
end

figure
subplot(2,2,1)
imagesc(nn_X1)
title({'tri factorization', X1.name})
subplot(2,2,2)
hist(nn_X1(:), 100)
title({'tri factorization', X1.name})

subplot(2,2,3)
imagesc(nn_X2)
title({'tri factorization', X2.name})
subplot(2,2,4)
hist(nn_X2(:), 100)
title({'tri factorization', X2.name})

%combine into data matrix
X = [nn_X1  nn_X2]; %now everything is [G -G I -I
figure
subplot(1,2,1)
imagesc(X)
title({'tri factorization', [X1.name ' and ' X2.name]})
subplot(1,2,2)
hist(X(:), 100)
title({'tri factorization', [X1.name ' and ' X2.name]})

F = W0;

figure
k2 = svd_singular_value_plot(X, ['tri factorization ' X1.name ' and ' X2.name]);

[m n] = size(X);
k1 = 2; %columns of F
%k2 = 50; %rows of G - want to be the same as number of labels




%initialize G
%how many clusters? let's say 50
G = zeros(k2,n);
tic
[idx c] = kmeans(X',k2,'Replicates',3);
toc
[temp indices] = sortrows(idx);


for i=1:n
    G(idx(i),i) = 1;
end
G_original = G;

G = G';

%% apply algorithm - need to pass in X
%TODO - Include G-G0 to updates (figure out how)

obj = 1e9;
obj_old = 1e10;
%F = rand(m,k1);
F=F+0;
[~, d1]= norm_cols(F+0,2); %NORMALIZE F - makes columns sum to 1
S = rand(k1,k2);
S = diag(d1)*S;
%G = rand(n,k2);
Z = zeros(m,n);
obj_traj = [];  % collect sequence of objective values
figure
while (obj_old - obj > 1e-6) && length(obj_traj) < parameters.max_iteration_count
    length(obj_traj)
%    upd_factor = ((X-Z)*G*S')./(F*S*G'*G*S' + lambF*sum(sum(F)) + epsilon);
%    Q = kappa*(upd_factor > 1 & F < kappa_tol);
%    if sum(sum(Q)) > 0 disp(sprintf('inadmissible zeros in F (iter = %d)',length(obj_traj))), end
%    F = (F + Q).*upd_factor;
%    [F,d1] = norm_cols(F,2); S = diag(d1)*S;    
    upd_factor = (F'*(X-Z)*G)./(F'*F*S*G'*G + lambS*sum(sum(S)) + epsilon);
    Q = kappa*(upd_factor > 1 & S < kappa_tol);
    if sum(sum(Q)) > 0 disp(sprintf('inadmissible zeros in S (iter = %d)',length(obj_traj))), end
    S = (S + Q).*upd_factor;
    
   upd_factor = ((X-Z)'*F*S)./(G*S'*F'*F*S + lambG*sum(sum(G)) + epsilon);
   Q = kappa*(upd_factor > 1 & G < kappa_tol);
   if sum(sum(Q)) > 0 disp(sprintf('inadmissible zeros in G (iter = %d)',length(obj_traj))), end
   G = (G + Q).*upd_factor;
   [G,d2] = norm_cols(G,2); S = S*diag(d2); %normalize rows of G, columns of G'
    
    T = X - F*S*G';
    Z = min(X, sign(T).*max(abs(T) - lambO,0));
    
    obj_old = obj
    obj = 0.5*(norm(X - F*S*G' - Z,'fro')^2 + lambF*sum(sum(F))^2 + ...
        lambS*sum(sum(S))^2 + lambG*sum(sum(G))^2) + lambO*sum(sum(abs(Z)));
    obj_traj = [obj_traj obj];
    plot(obj_traj)
    drawnow %pause(.1)
end
title({'tri factorization', [X1.name ' ' X2.name], sprintf('Total # of iteration = %d',length(obj_traj))})

disp(sprintf('Total # of iteration = %d',length(obj_traj)));

% normalize
[F,d1] = norm_cols(F,1);
[G,d2] = norm_cols(G,1); %normalizing by columns of G (rows of G') at the end gives group membership
%nGt = norm_cols(G', 1);
%G = nGt';
S = diag(d1)*S*diag(d2);

figure
subplot(1,3,1)
imagesc(F)
title({'tri factorization F', 'columns normalized', [X1.name ' and ' X2.name]})
colorbar
subplot(1,3,2)
imagesc(S)
title({'tri factorization S', [X1.name ' and ' X2.name]})
colorbar
subplot(1,3,3)
imagesc(G)
title({'tri factorization G', [X1.name ' and ' X2.name]})
colorbar

figure
SG = S*G';
imagesc(SG)
title({'tri factorization SG''', [X1.name ' and ' X2.name]})
colorbar
save('SG.mat', 'SG')

end