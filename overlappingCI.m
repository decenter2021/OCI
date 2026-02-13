% Description: Implementation of Kahan-family-optimal solution to the OCI
% problem devised in [1].
% Inputs: 
%   - Matrices H,R,C
%   - Cell of bounds Yb
%   - Metric 'trace' (default) or 'det'
%   - struct 'options' (optinal)
%       - 'verbose': solver verbose option (default: 1)
%       - 'normalization': 0: no normalization; 1: empirical normalization;
%          3: normalization that minimizes a cost function (default: 0)
%       - 'warning_on_numerical_problems': 0: throw error if terminates
%          with error 4 (happens if the orders of magnitude of the inputs 
%          are too different); 1: throw warning and output last iteration
%          of the solver as the solution (default: 0)
%       - 'cachesolvers': solver cachesolvers option (default: 0)
% Outputs: Matrices U,B,Y,K
% Assumptions and limitations: none
% Other m-files required: none
% MAT-files required: none
% Toolboxes required: YALMIP [2], SDPT3 [3,4]
% Authors: Leonardo Pedroso, Pedro Batista, W.P.M.H. (Maurice) Heemels
% Revision history:
%   - 20/01/2025 - First draft of full implementation (Leonardo Pedroso)
%   - 13/02/2026 - Robustness to numerical problem: Normalization and 
%                  Symmetrization  (Leonardo Pedroso)
% References:
% [1] L. Pedroso, P. Batista, and W.P.M.H. Heemels, "Overlapping Covariance 
%     Intersection: Fusion with Partial Structural Knowledge of Correlation
%     from Multiple Sources", 2026. (submitted)
% [2] J. Lofberg, "YALMIP: A toolbox for modeling and optimization in 
%     MATLAB." In 2004 IEEE international conference on robotics and
%     automation, pp. 284-289, 2004.
% [3] K.C. Toh, M.J. Todd, and R.H. Tutuncu, SDPT3 — a Matlab software 
%     package for semidefinite programming, Optimization Methods and 
%     Software, 11 (1999), pp. 545–581.
% [4] R.H Tutuncu, K.C. Toh, and M.J. Todd, Solving 
%     semidefinite-quadratic-linear programs using SDPT3, Mathematical 
%     Programming Ser. B, 95 (2003), pp. 189–217.

function [K,B,Y,U] = overlappingCI(H,R,C,Yb,metric,options)
    % Optional arguments
    if nargin < 6 || isempty(options)
        options.verbose = 1;
        options.normalization = 0;
        options.warning_on_numerical_problems = 0;
    end
    if ~isfield(options, 'verbose')
        options.verbose = 1;
    end
    if ~isfield(options, 'normalization')
        options.normalization = 0;
    end
    if ~isfield(options, 'warning_on_numerical_problems')
        options.warning_on_numerical_problems = 0;
    end
    if ~isfield(options,'cachesolvers')
        options.cachesolvers = 0;
    end
    % Sizes
    n = size(H,2);
    m = size(Yb{1,1},1);
    M = size(Yb,1);
    o = size(R,1);
    % Normalization constants
    if options.normalization == 2
        trY = 0;
        for i = 1:M
            trY = trY + trace(Yb{i,1});
        end
        trY = trY/(m*M);
        trR = trace(R)/o;
        normalization_cost_wrapper = @(alpha) normalization_cost(alpha,trY,trR,C,H);
        opts = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','none');
        [alpha,~,exitflag,output] = fminunc(normalization_cost_wrapper,...
                                            [mean([norm(C,'inf') norm(H,'inf')]);min(eig(R))],opts);
        if exitflag >= 0 
            alpha_HC = alpha(1);
            alpha_R = alpha(2);
        else
            warning('Optimized normalization failed: %s\nUsing default normalization.',output.message);
            alpha_HC = mean([norm(C,'inf') norm(H,'inf')]);
            alpha_R = min(eig(R));
        end
    elseif options.normalization == 1
        alpha_HC = mean([norm(C,'inf') norm(H,'inf')]);
        alpha_R = min(eig(R));
    else
        alpha_HC = 1;
        alpha_R = 1;
    end
    alpha_U = alpha_HC^2/alpha_R;
    alpha_B = 1/alpha_U;
    alpha_Y = alpha_U;    
    % Some robustness to numerical issues
    epsl_omega = 1e-6;
    epsl_U = 1e-10;
    Cn = C/alpha_HC; 
    Hn = H/alpha_HC; 
    Rn = (0.5/alpha_R)*(R+R');
    HRHn = ((H'/Rn)*H)/(alpha_HC^2);
    HRHn = 0.5*(HRHn + HRHn');
    CRCn = ((C'/Rn)*C)/(alpha_HC^2);
    CRCn = 0.5*(CRCn + CRCn');
    HRCn = ((H'/Rn)*C)/(alpha_HC^2);

    % Build SDP
    opts = sdpsettings('solver','mosek','verbose',options.verbose,'cachesolvers',options.cachesolvers);
    Un = sdpvar(n,n,'symmetric');
    Bn = sdpvar(n,n,'symmetric');
    Yn = sdpvar(m,m,'symmetric');
    omega = sdpvar(M,1);
    Yn_aux = zeros(m,m);
    for i = 1:M
        Yn_aux = Yn_aux + omega(i,1)*(0.5/alpha_Y)*(Yb{i,1}+Yb{i,1}');
    end

    % Constraints
    cntr = [];
    cntr = [cntr, HRHn - Un >= epsl_U*eye(n)];  % For numerical robustness of dual feasibility
    cntr = [cntr, Yn == Yn_aux];
    cntr = [cntr, omega >= epsl_omega, omega <= 1-epsl_omega, sum(omega)==1];
    cntr = [cntr, [Un HRCn; HRCn' Yn+CRCn] >= 0];
    cntr = [cntr, [Bn eye(n); eye(n) HRHn-Un] >= epsl_U*eye(2*n)];

    % Solve SDP
    if strcmp(metric, 'det')
        sol = optimize(cntr,-logdet(HRHn-Un),opts);
    elseif strcmp(metric, 'trace')
        sol = optimize(cntr,trace(Bn),opts);
    else
        warning("Unknown bound metric. Using trace as default.")
        sol = optimize(cntr,trace(Bn),opts);
    end   
    
    % Get solution
    Y = value(Yn)*alpha_Y;
    U = value(Un)*alpha_U;
    if strcmp(metric, 'det')
        aux = HRHn-value(Un);
        aux = 0.5*(aux+aux');
        [V,D] = eig(aux);
        B = (V*diag(1./max(D,epsl_inv))*V')*alpha_B;
    else
        B = value(Bn)*alpha_B;
    end

    % Compute K
    aux = (Rn\(Rn-Cn*pinv(value(Yn)+CRCn)*Cn'))/Rn;
    K = (1/alpha_HC)*((Hn'*aux*Hn)\(Hn'*aux));

    % Throw an error if not solved to optimality
    if sol.problem ~= 0
        if options.warning_on_numerical_problems && sol.problem == 4
            warning('OCI optimization did not meet termination criteria (code %d): %s\n%s', ...
                sol.problem, sol.info,yalmiperror(sol.problem));
        else
            error('OCI optimization failed (code %d): %s\n%s', ...
                sol.problem, sol.info,yalmiperror(sol.problem));
        end 
    end
end

function [f,g] = normalization_cost(alpha,trY,trR,C,H)
    f = 3*(trY*(alpha(2)/(alpha(1)^2))-1)^2+((trR/alpha(2))-1)^2 + ((norm(C,'inf')/alpha(1))-1)^2+((norm(H,'inf')/alpha(1))-1)^2;
    g = 2*[-6*(trY*(alpha(2)/(alpha(1)^2))-1)*trY*(alpha(2)/(alpha(1)^3))-(norm(C,'inf')/alpha(1)-1)*norm(C,'inf')/(alpha(1)^2) - (norm(H,'inf')/alpha(1)-1)*norm(H,'inf')/(alpha(1)^2);...
                       3*(trY*(alpha(2)/(alpha(1)^2))-1)*trY/(alpha(1)^2)-((trR/alpha(2))-1)*trR/(alpha(2)^2)];
end
