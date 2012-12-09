% X1, X2, labels must have the same number of rows
% labels should be OBSERVATIONSx1
function pipeline(X1, X2, labels, parameters)

figure
subplot(1,2,1)
svd_explore(X1, labels)
subplot(1,2,2)
svd_explore(X2, labels)

figure
subplot(1,2,1)
svd_singular_value_plot(X1, labels)
subplot(1,2,2)
svd_singular_value_plot(X2, labels)

% dual factorization
%  for each algorithm, output p-value for each of the k columns
%  Box plot for p-values
%  co-modules ???
%dual_factorization();

% tri-factorization
%   imagesc images for clusters?
%   svd explore on combined matrix
%tri_factorization();

% spectral clustering with shared nn
%   reduced dimensionality plot of clusters
%   purity measure for clusters
clustering(X1, X2, labels, parameters.cluster_count);

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

function svd_singular_value_plot(X, labels)
  X_data_norm = zscore(X.data);
  [u s v] = svd(X_data_norm, 'econ');
  %bar(diag(s))
  %hold on
  p = diag(s).^2;
  bar(p/sum(p))
  title({X.name, 'singular values'})
end