function x1 = plotCov2(invP)
    n = size(invP,1);
    % (for plotting we assume that X and Q are positive defnite)
    % Plotting parameters 
    Npoints = 100;
    theta = 0:(2*pi)/(Npoints-1):2*pi;
    % Elipse X1
    %[S1,D1,~] = eig(P); % X = SDS^-1
    R1 = chol(invP); % X = R'*R 
    x1 = zeros(n,Npoints);
    for i = 1:Npoints
        x1(:,i) = [cos(theta(i));sin(theta(i))];
        x1(:,i) = (R1)\x1(:,i);
    end
end