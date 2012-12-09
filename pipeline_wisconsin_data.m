function pipeline_wisconsin_data()
% loads an array called 'wisconsin'
load('.\data\wisconsin breast cancer\wisconsin.mat')

% break the wisconsin data into separate parts
X1 = []
X1.data = wisconsin(:, 1:10);
X1.name = 'Wisconsin X1';
X2.data = wisconsin(:, 11:20);
X2.name = 'Wisconsin X2';
labels = [];
labels.numeric = wisconsin(:, end);

[subject_count column_count] = size(wisconsin);
feature_count = column_count - 1;
labels.values = unique(labels.numeric)'
labels.names = {'malignant', 'benign'}
display(['wisconsin subject count: ' num2str(subject_count)])
display(['wisconsin feature count: ' num2str(feature_count)])
display(['wisconsin labels: ' num2str(labels.values)])

parameters.cluster_count = 2;
pipeline(X1, X2, labels, parameters);

end