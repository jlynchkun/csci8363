n=100; %number of data points
m = 5; %num of columns
k=2; %num clusters


x = rand(n,m); %2D data matrix
F = zeros(n,k);
G = zeros(k,m);

%cluster by rows
[idx c] = kmeans(x,3,'Replicates',10);
[temp indices] = sortrows(idx);
for i=1:n
    F(i,idx(i)) = 1;
end



%sort columns
[idx c] = kmeans(x',3,'Replicates',10);
[temp indices] = sortrows(idx);
for i=1:m
    G(idx(i),i) = 1;
end


S = F'*X*G;

% x_sorted=x(indices,:);
% figure;imagesc(x_sorted);
% figure;
% scatter(x_sorted(:,1),x_sorted(:,2),30,temp);
