function pipeline_wisconsin_data()
% loads an array called 'wisconsin'
load('.\data\wisconsin breast cancer\wisconsin.mat')

% break the wisconsin data into separate parts
df.skip = true;
df.X1 = []
df.X1.data = wisconsin(:, 1:10);
df.X1.name = 'Wisconsin X1';
df.X2.data = wisconsin(:, 11:20);
df.X2.name = 'Wisconsin X2';
df.labels = [];
df.labels.numeric = wisconsin(:, end);
correlations = corrcoef(df.X2.data);
B = correlations(end-size(df.X2.data, 2)+1:end,1:size(df.X2.data, 2)); %B is just the lower quarter of this correlation matrix - WEIRD DISTRIBUTION?
B(B<.7)=0; %make sparse - DO THIS???
df.X2_self_associations = B;
figure
imagesc(B)
title('Wisconsin X2 self associations')
df.dataname = 'Wisconsin_data_dual_factorization';
df.parameters.gamma1=20;
df.parameters.gamma2=10;
df.parameters.lambda1=.0001 * 0.1;
df.parameters.lambda2=.01;
df.parameters.lambda3=.0001;
df.parameters.iterations = 1;

[subject_count column_count] = size(wisconsin);
feature_count = column_count - 1;
df.labels.values = unique(df.labels.numeric)'
df.labels.names = {'malignant', 'benign'}

% break the wisconsin data into separate parts
clustering.skip = true;
clustering.X1 = []
clustering.X1.data = wisconsin(:, 1:10);
clustering.X1.name = 'Wisconsin X1';
clustering.X2.data = wisconsin(:, 11:20);
clustering.X2.name = 'Wisconsin X2';
clustering.labels = [];
clustering.labels.numeric = wisconsin(:, end);

[subject_count column_count] = size(wisconsin);
feature_count = column_count - 1;
clustering.labels.values = unique(clustering.labels.numeric)'
clustering.labels.names = {'malignant', 'benign'}
display(['wisconsin subject count: ' num2str(subject_count)])
display(['wisconsin feature count: ' num2str(feature_count)])
display(['wisconsin labels: ' num2str(clustering.labels.values)])

tf.skip = false;
tf.labels.numeric = wisconsin(:, end);
tf.X1.data = zscore(wisconsin(:, 1:10));
tf.X1.name = 'Wisconsin X1';
tf.X2.data = zscore(wisconsin(:, 11:20));
tf.X2.name = 'Wisconsin X2';

clustering.cluster_count = 2;
pipeline(clustering, df, tf);

end