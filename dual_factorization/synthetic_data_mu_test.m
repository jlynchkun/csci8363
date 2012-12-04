function synthetic_data_mu_test()

  synthetic_data

  K = 1
  
  [n,m1] = size(s_X1_nn)
  [n,m2] = size(s_X2_nn)

  s_W = rand(n,K)
  s_H1 = rand(K,m1) %initialize H to random numbers between 0 and 1
  s_H2 = rand(K,m2)

  L1 = 0.1
  L2 = 0.1
  L3 = 0.1
  r1 = 0.1
  r2 = 0.1
 

  [s_W,s_H1,s_H2] = multiplicative_update(s_X1_nn,s_X2_nn,s_W,s_H1,s_H2,s_A,s_B,L1,L2,L3,r1,r2,K)

  figure
  subplot(1,3,1)
  imagesc(s_W)
  title('s W')
  subplot(1,3,2)
  imagesc(s_H1)
  title('s H1')
  subplot(1,3,3)
  imagesc(s_H2)
  title('s H2')
end