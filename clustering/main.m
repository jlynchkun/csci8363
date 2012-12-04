


%%make matrix X
disp('creating data matrix')
% avoid this kind of colomn or row: sum == 0
index = find(sum(Gene_expression_data,1) == 0);
Gene_expression_data(:,index) = Gene_expression_data(:,index) + eps;
index = find(sum(Image_data,1) == 0);
Image_data(:,index) = Image_data(:,index) + eps;



%clean up gene data
%% Clean up gene data
[Common_gene ai bi] = intersect(GeneName(:,2), gene_name_ge);
Gene_expression_data = Gene_expression_data(:,ai);
preprocess = 0; %remove data with small absolute values and small variance

%preprocess genetic data - OPTIONAL - THIS WAS DONE IN ZHANG 2011
if preprocess
    mask1 = genelowvalfilter(Gene_expression_data','Percentile',60);
    mask2 = genevarfilter(Gene_expression_data','Percentile',30);    
    Gene_expression_data = Gene_expression_data(:,mask1 & mask2); %I THINK this applies both masks
    gene_network = gene_network(mask1&mask2,mask1&mask2);
else
    X2 = Gene_expression_data;
end






%MAKE DATA SAME SCALE - normalize by column
Gene_expression_data = zscore(Gene_expression_data);
Image_data = zscore(Image_data);


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
part1 = part1(I,:);
part2 = part2(I,:);


%combine into data matrix
X = [part1 part2];


%%snn
disp('creating adjacency matrix')
[A ranks distances] =  snn(X,30);

%%cluster
disp('clustering')
%calculate laplacian
L = diag(sum(A))-A;
[eigenvectors eigenvalues] = eig(A);
%number of clusters
k = 5;
eigenvectors_part = eigenvectors(:,[1:k]);


clusters = kmeans(eigenvectors_part,k);







function mask = genelowvalfilter(X,cutoff)



end

function mask = genevarfilter(X,cutoff)


end

