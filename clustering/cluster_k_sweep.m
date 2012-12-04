function cluster_k_sweep(k_range, varargin)
  k_count = numel(k_range);
  minimum_entropy = zeros(k_count, 1);
  wmean_entropy = zeros(k_count, 1);
  for ki = 1:numel(k_range)
      k = k_range(ki);
      [clusters cluster_entropy subject_count_per_cluster] = main(varargin{:}, k);
      subject_count_per_cluster
      minimum_entropy(ki) = min(cluster_entropy);
      wmean_entropy(ki) = mean(cluster_entropy .* (subject_count_per_cluster / sum(subject_count_per_cluster)))
  end
  figure
  plot(k_range, minimum_entropy, k_range, wmean_entropy)
  title({'minimum cluster entropy', 'weighted mean of cluster entropy'})
  legend({'min', 'weighted mean'})
end