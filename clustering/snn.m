function [A,ranks,distances] = snn(X,k)
%returns an adjacency matrix giving the number of shared nearest neighbors
%less than or equal to k

%A will contain number of shared nearest neighbors <= k
A = zeros(size(X,1),size(X,1));


distances = zeros(size(X,1),size(X,1));

ranks = zeros(size(X,1),size(X,1));

for i=1:size(X,1)-1
    for j=i+1:size(X,1)
        s = (X(i,:) - X(j,:)).^2;
        s = sum(s);
        distances(i,j) = s;
        distances(j,i) = s;
    end
end

%sort each row by nearness
for i=1:size(X,1)
    [Y I] = sort(distances(i,:));
    ranks(i,:) = I;
end


%find nearest neighbors
for i=1:size(X,1)-1
    for j=i+1:size(X,1)
        temp = length(intersect(ranks(i,2:k+1),ranks(j,2:k+1)));
        A(i,j) = temp;
        A(j,i) = temp;
    end
end