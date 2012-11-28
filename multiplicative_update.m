function  [W,H1,H2] = multiplicative_update(X1,X2,W,H1,H2,A,B,L1,L2,L3,r1,r2,K)
%multiplicative update from page i404 of Zhang et al 2011
%from paper, the following are equivalent
%gamma1 -> r1 = 20
%gamma2 -> r2 = 10
%lambda1 -> L1 = .0001
%lambda2 -> L2 = .01
%lambda3 -> L3 = ??? THIS IS FOR "learn B" parameter. 0 means don't learn
%B (i.e. don't push H1'*H2 to be like B)

errordiff = 1000000;

error = [];
while errordiff > 1000 %TODO: need a better (but princpled) ending point
    
    %this update is directly from Zhang's code, plus 
    W = W.*([H1 H2]*[X1 X2]')'./(W*([H1 H2]*[H1 H2]'+r1*eye(K))+eps);    %Update rule-1;  
%    HH1 = H1.*(W'*X1 + L2/2*H2*B')./((W'*W+r2*ones(K))*H1+eps);     %%Update H1 H2 simultaneously
    HH1 = H1.*(W'*X1 + L2/2*H2*B' + L3*2*H2*B')./((W'*W+r2*ones(K))*H1+eps + L3*2*H2*H2'*H1);     %%Update H1 H2 simultaneously
%    H2 = H2.*(W'*X2 + L1*H2*A + L2/2*H1*B)./((W'*W+r2*ones(K))*H2+eps); 
    H2 = H2.*(W'*X2 + L1*H2*A + L2/2*H1*B + 2*L3*H1*B)./((W'*W+r2*ones(K))*H2+eps + L3*2*H1*H1'*H2); 
    H1 = HH1;
    
    error = [error norm(X1-W*H1,'fro')^2 + norm(X2-W*H2,'fro')^2];

    plot(error)
    pause(.1) %to get the plot to show

    try
        errordiff = error(end-1)-error(end);
    end
end