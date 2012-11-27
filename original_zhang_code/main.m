%% load the real input data
getdata;
%% preprocess the real input data
OA = sparse(OA);
OB = sparse(OB);
%% Due to the nonegativity of the X1 and X2.
% we 'double' the Original matrix--OX1 and OX2
if ~isempty(find(OX1 < 0)) | ~isempty(find(OX2 < 0))
    X1 = [max(OX1,0) max(-OX1,0)];
    X2 = [max(OX2,0) max(-OX2,0)];
    
    A =[OA OA;
        OA OA];
    B =[OB OB;
        OB OB];
    isdouble = 1;
else
    X1 = OX1;
    X2 = OX2;
    A  = OA;
    B  = OB;
    isdouble = 0;
end
clear OX1 OX2 OA OB;
%% get factors by applying SNMNMF algorithm
L1 = 0.0001; L2 = 0.01; r1 = 10; r2 = 10; K = 50;
[W,H1,H2] = DNMF_residue_comodule(X1,X2,A,B,r1,r2,L1,L2,K);
%% get comodules(index)
tt = 7;
[Co_module] = DNMF_comodule(W, H1, H2, tt, isdouble);
save sampleComodule.mat Co_module;