%% load the real data: four matrices
%
% INPUT:
% OX1 (N,M1): input profile matrix
% OX2 (N,M2): input profile matrix
% OA (M2,M2): input adjacent matrix
% OB (M1,M2): input adjacent matrix

load('micro3_30.mat');%microRNA expression profile
OX1 = microdata3_30;
load('gene3_30.mat');%gene expression profile
OX2 = genedata3_30;
load('ppiMatrixTF.mat');% protein-protein interaction | TF interaction
OA = ppiMatrixTF;
load('geneMicroMatrix_v5.mat');%gene-MicroRNA interaction
OB = geneMicroMatrix;
clear microdata3_30 genedata3_30 ppiMatrixTF geneMicroMatrix