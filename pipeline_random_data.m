function pipeline_random_data(subject_count, feature_count, label_count)

%subject_count = 100;
%feature_count = 1000;
%label_count = 2;

raw_X1 = randn(subject_count, feature_count);
figure
raw_X2 = randn(subject_count, feature_count);

df.skip = false;
df.X1 = []
df.X1.data = zscore(raw_X1);
df.X1.name = 'Random X1';
df.X2.data = zscore(raw_X2);
df.X2.name = 'Random X2';
df.labels = [];
df.labels.numeric = randi(label_count, subject_count, 1);
correlations = corrcoef(df.X2.data);
B = correlations(end-size(df.X2.data, 2)+1:end,1:size(df.X2.data, 2)); %B is just the lower quarter of this correlation matrix - WEIRD DISTRIBUTION?
B(B<.7)=0; %make sparse - DO THIS???
df.X2_self_associations = B;
figure
imagesc(B)
title('Random X2 self associations')
df.dataname = 'Random_data_dual_factorization';
df.parameters.gamma1=20;
df.parameters.gamma2=10;
df.parameters.lambda1=.0001 * 0.1;
df.parameters.lambda2=.01;
df.parameters.lambda3=.0001;
df.parameters.iterations = 1;

df.labels.values = unique(df.labels.numeric)'
df.labels.names = {'malignant', 'benign'}

% break the wisconsin data into separate parts
clustering.skip = false;
clustering.X1 = []
clustering.X1.data = df.X1.data;
clustering.X1.name = 'Random X1';
clustering.X2.data = df.X2.data;
clustering.X2.name = 'Random X2';
clustering.labels = [];
clustering.labels.numeric = df.labels.numeric;

clustering.labels.values = unique(clustering.labels.numeric)'
clustering.labels.names = {'malignant', 'benign'}

tf.skip = false;
tf.labels.numeric = df.labels.numeric;
tf.X1.data = df.X1.data;
tf.X1.name = 'Random X1';
tf.X2.data = df.X2.data;
tf.X2.name = 'Random X2';

clustering.cluster_count = 2;
pipeline(clustering, df, tf);

end