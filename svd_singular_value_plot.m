function k = svd_singular_value_plot(X_data, X_name)
  X_data_norm = zscore(X_data);
  [u s v] = svd(X_data_norm, 'econ');
  %bar(diag(s))
  %hold on
  p = diag(s).^2;
  sum_p = sum(p);
  k = find(cumsum(p) >= 0.9 * sum_p, 1, 'first');
  bar(p/sum_p)
  title({X_name, 'singular values', ['k = ' num2str(k)]})
end