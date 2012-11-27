function [W,H1,H2] = DNMF_residue_comodule(X1,X2,A,B,r1,r2,L1,L2,K)
%
% INPUT:
% X1 (N,M1): input profile matrix
% X2 (N,M2): input profile matrix
% A (M2,M2): input adjacent matrix
% B (M1,M2): input adjacent matrix

% r1       : limit the growth of W
% r2       : limit the growth of H1 and H2
% L1       : weigh the must link constraints in A
% L2       : weigh the must link constraints in B
% K        : Number of components
%
% avoid this kind of colomn or row: sum == 0
index = find(sum(X1,1) == 0);
X1(:,index) = X1(:,index) + eps;

index = find(sum(X2,1) == 0);
X2(:,index) = X2(:,index) + eps;

index = find(sum(A,1) == 0);
A(:,index) = A(:,index) + eps;

index = find(sum(B,1) == 0);
B(:,index) = B(:,index) + eps;

% set the iteration number and initiate the output matrices
nloop = 5; 
verbose=1;

[n,m1] = size(X1);
[n,m2] = size(X2);

bestW=zeros(n,K);
bestH1=zeros(K,m1);
bestH2=zeros(K,m2);

bestobj1=1000000000;
bestobj2=1000000000;
fid = fopen(['record_K' int2str(K) '_L1=' num2str(L1) '_L2=' num2str(L2) '_r1=' num2str(r1) '_r2=' num2str(r2) '.txt'],'wt+');
for iloop=1:nloop
    if verbose 
        fprintf(fid,' iteration %d\n',iloop); 
    end 
    
    maxiter=500; 
    speak=1; 
    [W,H1,H2]=DNMF_mm(X1, X2, A, B, r1, r2, L1, L2, K, maxiter, speak, fid);
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
W = bestW; H1 = bestH1; H2 = bestH2;

