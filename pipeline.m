% X1, X2, labels must have the same number of rows
% labels should be OBSERVATIONSx1
function tri = pipeline(clustering_data, df_data, tf_data)

figure
subplot(1,2,1)
svd_singular_value_plot(clustering_data.X1.data, clustering_data.X1.name)
subplot(1,2,2)
svd_singular_value_plot(clustering_data.X2.data, clustering_data.X2.name)

figure
subplot(1,2,1)
svd_explore(clustering_data.X1, clustering_data.labels)
subplot(1,2,2)
svd_explore(clustering_data.X2, clustering_data.labels)

% dual factorization
%  for each algorithm, output p-value for each of the k columns
%  Box plot for p-values
%  co-modules ???
if df_data.skip == true
  display('skipping dual factorization')
else
  dual_factorization(df_data.X1, df_data.X2, df_data.X2_self_associations, df_data.labels, df_data.parameters, df_data.dataname);
end

% tri-factorization
%   imagesc images for clusters?
%   svd explore on combined matrix
%tri_factorization();

% spectral clustering with shared nn
%   reduced dimensionality plot of clusters
%   purity measure for clusters
if clustering_data.skip == true
  display('skipping clustering')
else
  clustering(clustering_data.X1, clustering_data.X2, clustering_data.labels, clustering_data.cluster_count);
end

tri.F = [];
tri.S = [];
tri.G = [];
if tf_data.skip == true
  display('skipping tri factorization')
else
  [tri.F tri.S tri.G] = tri_factorization(tf_data.X1, tf_data.X2, tf_data.labels, tf_data.parameters);
end

end

function svd_explore(X, labels)
  [obs_count feature_count] = size(X.data)
  % zscore the columns
  X_data_norm = zscore(X.data);
  [u s v] = svd(X_data_norm, 'econ');
  dp = X_data_norm * v(:,1:2) * inv(s(1:2,1:2));
  obs_numbers = 1:obs_count;
  for i = 1:numel(labels.values)
    label_numeric_value = labels.values(i)
    obs_with_label = labels.numeric == label_numeric_value;
    scatter(dp(obs_with_label,1), dp(obs_with_label,2), 30, labels.numeric(obs_with_label), 'filled');
    %text(dp(obs_with_label,1), dp(obs_with_label,2), num2str(obs_numbers(obs_with_label)'))
    hold on
  end
  legend(labels.names, 'Location', 'NorthOutside')
  title(X.name)
end

function svd_remove_outliers(X, labels)
  [obs_count feature_count] = size(X.data)
  % zscore the columns
  X_data_norm = zscore(X.data);
  [u s v] = svd(X_data_norm, 'econ');
  dp = X_data_norm * v(:,1:3) * inv(s(1:3,1:3));
  g = gmdistribution.fit(dp(:, 2:3), 2, 'Replicates', 5);
  p = posterior(g, dp(:, 2:3));
  outliers = p(:, 1) < 0.001;
  scatter(dp(:,1), dp(:,2), 30, p(:, 1), 'filled');
  %hold on
  %ezcontour(@(x,y)pdf(g,[x y]), xlim(), ylim(), 100);
  %for i = 1:numel(labels.values)
    %label_numeric_value = labels.values(i)
    %obs_with_label = labels.numeric == label_numeric_value;
    %text(dp(obs_with_label,1), dp(obs_with_label,2), num2str(obs_numbers(obs_with_label)'))
    %hold on
  %end
  %legend(labels.names, 'Location', 'NorthOutside')
  title({'Outliers', X.name, [num2str(nnz(outliers)) ' outliers in ' num2str(obs_count) ' observations']})
  colorbar
end
