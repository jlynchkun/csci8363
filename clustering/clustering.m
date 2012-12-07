
function [clusters, cluster_entropy, subject_count_per_cluster] = main(Gene_expression_data, Image_data, gene_name_ge, GeneName, gene_network, Clinical_numeric, survival, target_cluster_count)

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
part1 = Gene_expression_data(I,:);
part2 = Image_data(I,:);


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
k = target_cluster_count;
eigenvectors_part = eigenvectors(:,[1:k]);

clusters = kmeans_pp(eigenvectors_part',k)';
cluster_labels = unique(clusters)
cluster_label_count = numel(cluster_labels);
clinical_labels = unique(labels)
clinical_label_count = numel(clinical_labels);
classification_results = zeros(cluster_label_count, clinical_label_count);
cluster_entropy = zeros(cluster_label_count, 1);
subject_count_per_cluster = zeros(cluster_label_count, 1);
for ci = 1:numel(cluster_labels)
   cluster_label = cluster_labels(ci);
   subject_in_cluster_ci = clusters == cluster_label;
   subject_count_per_cluster(ci) = nnz(subject_in_cluster_ci);
   for cj = 1:numel(clinical_labels)
       clinical_label = clinical_labels(cj);
       subject_has_clinical_label_cj = labels == clinical_label;
       % how many subjects have cluster label ci and clinical label cj?
       classification_results(ci, cj) = nnz(subject_in_cluster_ci & subject_has_clinical_label_cj);
       display([num2str(classification_results(ci, cj)) ' subjects have cluster label ' num2str(cluster_label) ' and clinical label ' num2str(clinical_label)])
   end
   cluster_ci_p = classification_results(ci, :) / nnz(subject_in_cluster_ci);
   cluster_entropy(ci) = -sum(cluster_ci_p .* log2(cluster_ci_p));
end
display('cluster entropy: ')
cluster_entropy
end



