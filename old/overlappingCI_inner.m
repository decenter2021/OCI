function [K,B,Y_outer,Y_inner,Ybinner] = overlappingCI_inner(H,R,C,Yb)
    % Sizes
    n = size(H,2);
    o = size(H,1);
    % m = size(Yb{1,1},1);
    M = size(Yb,1);
    % Compute Yinner
    Ybinner = cell(M,1);
    for i = 1:M
        Ybinner{i,1} = inv(R)-((R\C)/(Yb{i,1}+(C'/R)*C))*(C'/R); %inv(R+(C/Yb{i,1})*C');
    end
    % Build SDP
    omega = sdpvar(M,1);
    Q1 = zeros(o,o);
    for i = 1:M
        Q1 = Q1 + omega(i,1)*Ybinner{i,1};
    end
    gamma = sdpvar(1,1); 
    Q2 = sdpvar(n,n);
    cntr = [gamma>=0, Q2>=0];
    cntr = [cntr, omega>=0, omega<=1, ones(1,M)*omega == 1];
    cntr = [cntr, [Q2 eye(n); eye(n) H'*Q1*H]>=0];
    cntr = [cntr, gamma>= trace(Q2)];
    optimize(cntr,gamma);
    Y_inner = value(Q1);
    Y_outer = invertBound(C,inv(inv(Y_inner)-R));
    B = value(Q2);
    % Compute K
    K = (H'*Y_inner*H)\(H'/Y_inner);
end