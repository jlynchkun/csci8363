function pipeline_gene_image(verbose)

% these should be arguments? Clinical_header, Clinical_numeric,
% Clinical_string, Gene_name, Gene_expression_data
load('./data/Breast_cancer_image.mat');
load('./data/ppiMatrixTF.mat');
display(['clinical header size ' num2str(size(Clinical_header))])
display(['clinical numeric size ' num2str(size(Clinical_numeric))])
% Clinical_header gives labels for the columns of Clinical_numeric

display(['clinical string size ' num2str(size(Clinical_string))])
if verbose 
  disp(Clinical_header')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% avoid this kind of colomn or row: sum == 0
zero_sum_index = find(sum(Gene_expression_data,1) == 0);
Gene_expression_data(:,zero_sum_index) = Gene_expression_data(:,zero_sum_index) + eps;
zero_sum_index = find(sum(Image_data,1) == 0);
Image_data(:,zero_sum_index) = Image_data(:,zero_sum_index) + eps;

%clean up gene data
[~, all_gene_count] = size(gene_name_ge);
[Common_gene ai bi] = intersect(GeneName(:,2), gene_name_ge);
named_gene_count = numel(ai)
Gene_expression_data = Gene_expression_data(:,ai);
preprocess = true; %remove data with small absolute values and small variance
%preprocess genetic data - OPTIONAL - THIS WAS DONE IN ZHANG 2011
if preprocess
    mask1 = genelowvalfilter(Gene_expression_data','Percentile',60);
    mask2 = genevarfilter(Gene_expression_data','Percentile',30);    
    Gene_expression_data = Gene_expression_data(:,mask1 & mask2); %I THINK this applies both masks
    gene_network_size = size(ppiMatrixTF)
    ppiMatrixTF = ppiMatrixTF(mask1&mask2,mask1&mask2);
    gene_network_size = size(ppiMatrixTF)
end

%SORT BY survival/metastasis status
survival = survival/12; %survival from month to year scale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

metastasis_labels = [];
metastasis_label_index = 6;
% these are 0 = no metastasis, 1 = metastasis
% change to 1 = non-metastasis, -1 = metastasis
metastasis_labels.numeric = Clinical_numeric(:, metastasis_label_index);
metastasis_labels.numeric(metastasis_labels.numeric == 1) = -1;
metastasis_labels.numeric(metastasis_labels.numeric == 0) =  1;
% no string labels for column 6 of Clinical_string
metastasis_labels.values = unique(metastasis_labels.numeric)'
metastasis_labels.names = {'metastasis', 'non-metastasis'}

[gene_expression_subject_count gene_expression_feature_count] = size(Gene_expression_data);
display(['gene expression subject count: ' num2str(gene_expression_subject_count)])
display(['gene expression feature count: ' num2str(gene_expression_feature_count)])

[image_data_subject_count image_data_feature_count] = size(Image_data);
display(['image data subject count: ' num2str(image_data_subject_count)])
display(['image data feature count: ' num2str(image_data_feature_count)])

for label_i = 1:numel(metastasis_labels.values)
  label = metastasis_labels.values(label_i);
  subject_count_for_label_i = nnz(metastasis_labels.numeric == metastasis_labels.values(label_i));
  display([num2str(subject_count_for_label_i) ' subjects have label "' metastasis_labels.names{label_i} '"'])
end

[~, sorted_metastasis_label_indices] = sort(metastasis_labels.numeric);
X1 = [];
X1.data = zscore(Image_data(sorted_metastasis_label_indices, :));
X1.name = 'Image data';
X2 = [];
X2.data = zscore(Gene_expression_data(sorted_metastasis_label_indices, :));
X2.name = 'Gene expression data';
figure
subplot(1,2,1)
imagesc(X1.data)
title({X1.name, 'sorted'})
subplot(1,2,2)
imagesc(X2.data)
title({X2.name, 'sorted'})

parameters.cluster_count = 5;
pipeline(X1, X2, metastasis_labels, parameters)


end