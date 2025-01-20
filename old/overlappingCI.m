function [K,B,Y,U] = overlappingCI(H,R,C,Yb)
    % Sizes
    n = size(H,2);
    epsilon = 1e-4;
    % o = size(H,1);
    m = size(Yb{1,1},1);
    M = size(Yb,1);
    % Build SDP
    omega = sdpvar(M,1);
    Y = zeros(m,m);
    for i = 1:M
        Y = Y + omega(i,1)*Yb{i,1};
    end
    gamma = sdpvar(1,1); 
    U = sdpvar(m,m);
    B = sdpvar(n,n);
    cntr = [gamma>=0, U>=0, B>=0]; %, Y>= epsilon*eye(m)];
    cntr = [cntr, omega>=0, omega<=1, ones(1,M)*omega == 1];
    cntr = [cntr, [U eye(m); eye(m) Y+(C'/R)*C]>=0];
    cntr = [cntr, [B eye(n); eye(n) (H'/R)*(R-C*U*C')*(R\H)]>=0];
    cntr = [cntr, gamma>= trace(B)];
    optimize(cntr,gamma);
    Y = value(Y);
    U = value(U);
    B = value(B);
    % Compute K
    K = ((H'/R)*(R-C*U*C')*(R\H))\(((H'/R)*(R-C*U*C'))/R);
    %K = (((H'/(R+(C/Y_outer)*C'))*H)\H')/(R+(C/Y_outer)*C');
end