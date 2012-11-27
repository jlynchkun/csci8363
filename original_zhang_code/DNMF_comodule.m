function [Co_module] = DNMF_comodule(W, H1, H2, tt, isdouble)
%
% INPUT
% W         : common basis matrix
% H1        : coefficient matrix
% H2        : coefficient matrix
% tt        : a given threshold for z-score.
% isdouble  : a flag
%
% OUTPUT
% Co_module : the index list of microRNAs and Genes in a Co-module.
%
% Compute the mean(meadia) and std in columns of W and rows in H1, H2 to determine
% the module member and output the Co-module based on W and H1, H2
% matrices.
%
m1 = size(H1,2);
m2 = size(H2,2);
n = size(W,1);
K = size(W,2);

MW = mean(W,1);     MH1 = mean(H1,2);   MH2 = mean(H2,2);
VW = std(W,0,1);    VH1 = std(H1,0,2);  VH2 = std(H2,0,2);

% Co-Module
if isdouble
    for i = 1:K
        c1 = find(H1(i,:) > MH1(i) + tt*VH1(i));
        c1(c1 > m1/2) = c1(c1 > m1/2) - m1/2;% tranform the double microRNA index into origin index
        c1 = unique(c1);
        
        c2 = find(H2(i,:) > MH2(i) + tt*VH2(i));
        c2(c2 > m2/2) = c2(c2 > m2/2) - m2/2; % tranform the double gene index into origin index
        c2 = unique(c2);
        
        Co_module{i,1}=c1; Co_module{i,2}=c2;
    end
else
    for i = 1:K
        c1 = find(H1(i,:) > MH1(i) + tt*VH1(i));
        c2 = find(H2(i,:) > MH2(i) + tt*VH2(i));
        Co_module{i,1}=c1; Co_module{i,2}=c2;
    end
end