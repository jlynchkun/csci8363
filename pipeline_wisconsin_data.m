function pipeline_wisconsin_data(noise_feature_count)
% loads an array called 'wisconsin'
load('.\data\wisconsin breast cancer\wisconsin.mat')

%cell1_data = 
[subject_count ~] = size(wisconsin);

% break the wisconsin data into separate parts
df.skip = true;
df.X1 = []
df.X1.data = zscore(wisconsin(:, 1:10));
df.X1.name = 'Wisconsin X1';
df.X2.data = zscore([wisconsin(:, 11:20) randn(subject_count, noise_feature_count)]);
df.X2.name = ['Wisconsin X2 with ' num2str(noise_feature_count) ' noise features'];
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
df.parameters.K = 25;
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
clustering.skip = false;
clustering.X1 = []
clustering.X1.data = df.X1.data;
clustering.X1.name = df.X1.name;
clustering.X2.data = df.X2.data;
clustering.X2.name = df.X2.name;
clustering.labels = [];
clustering.labels.numeric = wisconsin(:, end);

[subject_count column_count] = size(wisconsin);
feature_count = column_count - 1;
clustering.labels.values = unique(clustering.labels.numeric)'
clustering.labels.names = {'malignant', 'benign'}
display(['wisconsin subject count: ' num2str(subject_count)])
display(['wisconsin feature count: ' num2str(feature_count)])
display(['wisconsin labels: ' num2str(clustering.labels.values)])

tf.skip = true;
tf.labels.numeric = wisconsin(:, end);
tf.X1.data = df.X1.data;
tf.X1.name = df.X1.name;
tf.X2.data = df.X2.data;
tf.X2.name = df.X2.name;
tf.parameters.max_iteration_count = 100;

clustering.cluster_count = 2;
pipeline(clustering, df, tf);

end