function [R, E] = gauss_jordan(A)
    % GAUSS_JORDAN performs Gauss-Jordan elimination and returns the RREF and 
    % the product of elementary matrices that reduce A to its RREF.
    % Input:
    %   A - The input matrix to be reduced.
    % Output:
    %   R - The reduced row echelon form (RREF) of matrix A.
    %   E - The product of elementary matrices that reduce A to R.

    % Initialize the size of the matrix
    [m, n] = size(A);
    
    % Initialize the matrix of elementary matrices as an identity matrix
    E = eye(m);
    
    % Perform Gauss-Jordan elimination
    for i = 1:min(m, n)
        % Find pivot row and column
        [~, pivot_row] = max(abs(A(i:m, i)));
        pivot_row = pivot_row + i - 1;
        
        % Swap rows in A and in the elementary matrix
        if pivot_row ~= i
            A([i, pivot_row], :) = A([pivot_row, i], :);
            E([i, pivot_row], :) = E([pivot_row, i], :);
        end
        
        % Scale the pivot row to make the pivot element 1
        if A(i, i) ~= 0
            scale_factor = A(i, i);
            A(i, :) = A(i, :) / scale_factor;
            E(i, :) = E(i, :) / scale_factor;
        end
        
        % Eliminate all entries in column i (except for the pivot itself)
        for j = 1:m
            if j ~= i
                factor = A(j, i);
                A(j, :) = A(j, :) - factor * A(i, :);
                E(j, :) = E(j, :) - factor * E(i, :);
            end
        end
    end
    
    % The resulting matrix A is now in reduced row echelon form
    R = A;
end