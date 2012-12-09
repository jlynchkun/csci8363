% data should be zscored
%function [clusters, cluster_entropy, subject_count_per_cluster] = clustering(Gene_expression_data, Image_data, gene_name_ge, GeneName, gene_network, Clinical_numeric, survival, target_cluster_count)
function [clusters, cluster_entropy, subject_count_per_cluster] = clustering(X1, X2, labels, target_cluster_count)

%MAKE DATA SAME SCALE - normalize by column
%%%Gene_expression_data = zscore(Gene_expression_data);
%%%Image_data = zscore(Image_data);



%put everything in the right order
%[Y I] = sort(labels,'descend'); %sort labels by class - I is indices
%labels_sorted = labels(I);
%W0 = [labels_sorted > 0 labels_sorted < 0];
%survival = survival(I,:);

%part1 = Gene_expression_data(I,:);
%part2 = Image_data(I,:);


%combine into data matrix
X = [X1.data X2.data];


%%snn
disp('creating adjacency matrix')
[A ranks distances] = snn(X,30);
figure
imagesc(ranks)
colorbar
title({'snn', X1.name, X2.name})

%%cluster
disp('clustering')
%calculate laplacian
L = diag(sum(A))-A;
[eigenvectors eigenvalues] = eig(A);
%number of clusters
k = target_cluster_count;
eigenvectors_part = eigenvectors(:,(1:k));

clusters = kmeans_pp(eigenvectors_part',k)';
cluster_labels = unique(clusters)
cluster_label_count = numel(cluster_labels);
clinical_labels = labels.values
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
       subject_has_clinical_label_cj = labels.numeric == clinical_label;
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



