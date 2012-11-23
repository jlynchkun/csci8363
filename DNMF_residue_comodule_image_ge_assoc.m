function [Initial_W,W1,H11,H12,W2,H21,H22,W3,H31,H32,W_ge,H_ge,W_image,H_image] = DNMF_residue_comodule_image_ge_assoc(X1,X2,W0, H10, H20, A,B,L1,L2,r1,r2,r3,K,opt)
%
% INPUT:
% X1 (N,M1): input profile matrix
% X2 (N,M2): input profile matrix
% A (M1,M1): input adjacent matrix

% r1       : limit the growth of W
% r2       : limit the growth of H1 and H2
% L1       : weigh the must link constraints in A
% K        : Number of components
%
% avoid this kind of colomn or row: sum == 0
index = find(sum(X1,1) == 0);
X1(:,index) = X1(:,index) + eps;

index = find(sum(X2,1) == 0);
X2(:,index) = X2(:,index) + eps;

% edited by taehyun        
% index = find(sum(A,1) == 0);
% A(:,index) = A(:,index) + eps;
% 
% index = find(sum(B,1) == 0);
% B(:,index) = B(:,index) + eps;

% set the iteration number and initiate the output matrices
% nloop = 5; 
nloop = 1; 
verbose=1;

[n,m1] = size(X1);
[n,m2] = size(X2);

bestW=zeros(n,K);
bestH1=zeros(K,m1);
bestH2=zeros(K,m2);



%%   OPTION 1 - don't include B gene-image correlation matrix at all
disp('Option 1')
bestobj1=1000000000; %objective function - 
bestobj2=1000000000;
r1 = 15;
r2 = 15;
fid = fopen(['record_K' int2str(K) '_r1=' num2str(r1) '_r2=' num2str(r2) '.txt'],'wt+');
for iloop=1:nloop
    if verbose 
        fprintf(fid,' iteration %d\n',iloop); 
    end 

%edited by taehyun    
%     maxiter=500;
    maxiter=50;
    speak=1; 
    [Initial_W,W,H1,H2]=DNMF_mm_image_ge_wo_assoc(X1, X2, W0, H10, H20, A, B, L1, L2, r1, r2, r3, K, maxiter, speak, fid, opt);
    % compute residue
    newobj1 = sum(sum((X1-W*H1).^2));
    newobj2 = sum(sum((X2-W*H2).^2));
    
    if (newobj1<bestobj1)||(newobj2<bestobj2)        
        bestobj1 = newobj1;
        bestobj2 = newobj2;  
        bestW = W;
        bestH1 = H1;
        bestH2 = H2;
    end
end
fclose(fid);
%  compute the modules according to bestW, bestH1 and bestH2
W1 = bestW; H11 = bestH1; H12 = bestH2;








%%  OPTION 2, INCLUDE B = G1^T*G2
disp('Option 2')
bestobj1=1000000000;
bestobj2=1000000000;
opt = 6;
r1 = 1;
r2 = 1;
fid = fopen(['record_K_assoc' int2str(K) '_r1=' num2str(r1) '_r2=' num2str(r2) '.txt'],'wt+');
for iloop=1:nloop
    if verbose 
        fprintf(fid,' iteration %d\n',iloop); 
    end 

%edited by taehyun    
%     maxiter=500;
    maxiter=50;
    speak=1; 
    [Initial_W,W,H1,H2]=DNMF_mm_image_ge_use_assoc_corr(X1, X2, W0, H10, H20, A, B, L1, L2, r1, r2, r3, K, maxiter, speak, fid, opt);
    % compute residue
    newobj1 = sum(sum((X1-W*H1).^2));
    newobj2 = sum(sum((X2-W*H2).^2));
    
    if (newobj1<bestobj1)||(newobj2<bestobj2)        
        bestobj1 = newobj1;
        bestobj2 = newobj2;  
        bestW = W;
        bestH1 = H1;
        bestH2 = H2;
    end
end
fclose(fid);
%  compute the modules according to bestW, bestH1 and bestH2
W2 = bestW; H21 = bestH1; H22 = bestH2;

%% B = G1^T*G2 +++ ||B-B0||
disp('Option 3')
bestobj1=1000000000;
bestobj2=1000000000;
opt = 6;
r1 = 1;
r2 = 1;
fid = fopen(['record_K_assoc' int2str(K) '_r1=' num2str(r1) '_r2=' num2str(r2) '.txt'],'wt+');
for iloop=1:nloop
    if verbose 
        fprintf(fid,' iteration %d\n',iloop); 
    end 

%edited by taehyun    
%     maxiter=500;
    maxiter=50;
    speak=1; 
    [Initial_W,W,H1,H2]=DNMF_mm_image_ge_learn_assoc_corr(X1, X2, W0, H10, H20, A, B, L1, L2, r1, r2, r3, K, maxiter, speak, fid, opt);
    % compute residue
    newobj1 = sum(sum((X1-W*H1).^2));
    newobj2 = sum(sum((X2-W*H2).^2));
    
    if (newobj1<bestobj1)||(newobj2<bestobj2)        
        bestobj1 = newobj1;
        bestobj2 = newobj2;  
        bestW = W;
        bestH1 = H1;
        bestH2 = H2;
    end
end
fclose(fid);
%  compute the modules according to bestW, bestH1 and bestH2
W3 = bestW; H31 = bestH1; H32 = bestH2;

%% Only image
disp('Only Image Data')
maxiter = 50;
[W, H] = nnmf_by_taehyun(X1, W0, H10, maxiter);
W_ge = W; H_ge = H; 

%% Only GE
disp('Only Genetic Data')
[W, H] = nnmf_by_taehyun(X2, W0, H20, maxiter);
W_image = W; H_image = H; 


