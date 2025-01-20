function Y = invertBound(W,Xinv)
    % -- Description --
    % Given a bound WEW' <= X on E
    % We want to find a bound Y on E^-1, i.e., 
    % E^-1 >= Y
    Cpar = null(null(W)');
    Y = Cpar*(W*Cpar)'*Xinv*(W*Cpar)*Cpar';
    
    % % --- Old ---
    % % Sizes
    % o = size(W,1);
    % m = size(W,2);
    % % Assemble basis chnage matrix
    % Cperp = null(W);
    % Cpar = null(Cperp');
    % S = [Cpar Cperp];
    % T = (S\[W;zeros(m-o,m)])*S;
    % Y = (T(:,1:o)*Xinv)*T(:,1:o)';
end