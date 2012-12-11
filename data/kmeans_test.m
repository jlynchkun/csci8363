% n=100; %number of data points
% m = 5; %num of columns
% k=2; %num clusters
% 
% 
% x = rand(n,m); %2D data matrix
% F = zeros(n,k);
% G = zeros(k,m);
% 
% %cluster by rows
% [idx c] = kmeans(x,3,'Replicates',10);
% [temp indices] = sortrows(idx);
% for i=1:n
%     F(i,idx(i)) = 1;
% end
% 
% 
% 
% %sort columns
% [idx c] = kmeans(x',3,'Replicates',10);
% [temp indices] = sortrows(idx);
% for i=1:m
%     G(idx(i),i) = 1;
% end
% 
% 
% S = F'*x*G;
% 
% % x_sorted=x(indices,:);
% % figure;imagesc(x_sorted);
% % figure;
% % scatter(x_sorted(:,1),x_sorted(:,2),30,temp);
% 






%create X.
image_transformed = zscore(Image_data);
gene_transformed = zscore(Gene_expression_data);

X_pos = [max(gene_transformed,0) max(-gene_transformed,0) (max(image_transformed,0)) (max(-image_transformed,0))];

%F is based on our labels
%orthogonal - either one label or the other
F = [max(labels,0) max(-labels,0)];

%sort X,F based on labels. Normalize X.
[F_sorted indices] = sortrows(F,1);
X_sorted = X_pos(indices,:);

%initialize G
G = zeros(2,size(X_sorted,2));
[idx c] = kmeans(X_sorted',2,'Replicates',10);
[temp indices] = sortrows(idx);


for i=1:size(X_sorted,2)
    G(idx(i),i) = 1;
end
G_original = G;

%update G
s = size(F_sorted'*F_sorted);

for i=1:1000
G = G.*(F_sorted'*X_sorted)./((F_sorted'*F_sorted+.01*eps)*G);
end
