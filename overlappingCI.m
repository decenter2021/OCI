% Description: Implementation of Kahan-family-optimal solution to the OCI
% problem devised in [1].
% Inputs: 
%   - Matrices H,R,C
%   - Cell of bounds Yb
%   - Metric 'trace' (default) or 'det'
%   - verbose (optinal)
% Outputs: Matrices U,B,Y,K
% Assumptions and limitations: none
% Other m-files required: none
% MAT-files required: none
% Toolboxes required: YALMIP [2], SDPT3 [3,4]
% Authors: Leonardo Pedroso, Pedro Batista, W.P.M.H. (Maurice) Heemels
% Revision history:
%   - 20/01/2025 - First draft of full implementation (Leonardo Pedroso)
% References:
% [1] L. Pedroso, P. Batista, and W.P.M.H. Heemels, "Overlapping Covariance 
%     Intersection: Fusion with Partial Structural Knowledge of Correlation
%     from Multiple Sources", 2025. (submitted)
% [2] J. Lofberg, "YALMIP: A toolbox for modeling and optimization in 
%     MATLAB." In 2004 IEEE international conference on robotics and
%     automation, pp. 284-289, 2004.
% [3] K.C. Toh, M.J. Todd, and R.H. Tutuncu, SDPT3 — a Matlab software 
%     package for semidefinite programming, Optimization Methods and 
%     Software, 11 (1999), pp. 545–581.
% [4] R.H Tutuncu, K.C. Toh, and M.J. Todd, Solving 
%     semidefinite-quadratic-linear programs using SDPT3, Mathematical 
%     Programming Ser. B, 95 (2003), pp. 189–217.

function [K,B,Y,U] = overlappingCI(H,R,C,Yb,metric,verbose)
    % Optional arguments
    if nargin < 6 || isempty(verbose)
        verbose = 1;
    end
    % Sizes
    n = size(H,2);
    m = size(Yb{1,1},1);
    M = size(Yb,1);
    % Build SDP
    opts = sdpsettings('solver','mosek','verbose',verbose);
    omega = sdpvar(M,1);
    Y = zeros(m,m);
    for i = 1:M
        Y = Y + omega(i,1)*Yb{i,1};
    end
    U = sdpvar(n,n);
    B = sdpvar(n,n);
    cntr = [U>=0, B>=0];
    cntr = [cntr, omega>=0, omega<=1, ones(1,M)*omega == 1];
    cntr = [cntr, [U (H'/R)*C; ((H'/R)*C)' Y+(C'/R)*C]>=0];
    cntr = [cntr, [B eye(n); eye(n) (H'/R)*H-U]>=0];
    % Solve SDP
    if strcmp(metric, 'det')
        sol = optimize(cntr,-logdet((H'/R)*H-U),opts);
    elseif strcmp(metric, 'trace')
        sol = optimize(cntr,trace(B),opts);
    else
        warning("Unknown bound metric. Using trace as default.")
        sol = optimize(cntr,trace(B),opts);
    end   
    % Throw an error if not solved to optimality
    if sol.problem ~= 0
        error('OCI optimization failed (code %d): %s\n%s', ...
            sol.problem, sol.info,yalmiperror(sol.problem));
    end
    % Get solution
    Y = value(Y);
    U = value(U);
    if strcmp(metric, 'det')
        B = inv((H'/R)*H-U);
    else
        B = value(B);
    end
    % Compute K
    aux = (R\(R-C*pinv(Y+(C'/R)*C)*C'))/R;
    K = (H'*aux*H)\(H'*aux);
end