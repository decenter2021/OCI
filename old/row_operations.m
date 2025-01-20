
function [rref_mat, row_ops, elem_matrices] = row_operations(A)
% Input: A is the augmented matrix
% Output: rref_mat is the RREF of A
% row_ops is a cell array containing the row operations performed to obtain the RREF form
% elem_matrices is a cell array containing the corresponding elementary matrices for each row operation
% Initialization
[r, c] = size(A);
rref_mat = A;
row_ops = cell(r-1, 1);
elem_matrices = cell(r-1, 1);
% Perform Gaussian elimination to obtain upper triangular matrix
for i = 1:r-1
    % Find the first non-zero entry in the ith column
    pivot_row = i;
    while (pivot_row <= r) && (rref_mat(pivot_row, i) == 0)
        pivot_row = pivot_row + 1;
    end
    % If there is no non-zero entry in the ith column, move to the next column
    if pivot_row > r
        continue;
    end
    % Swap rows if necessary to bring pivot to the ith row
    if pivot_row > i
        rref_mat([i pivot_row], :) = rref_mat([pivot_row i], :);
        row_ops{i} = sprintf('R%d <-> R%d', i, pivot_row);
        elem_matrices{i} = eye(r);
        elem_matrices{i}([i pivot_row], [i pivot_row]) = elem_matrices{i}([pivot_row i], [pivot_row i]);
    end
    % Scale the ith row to make the pivot equal to 1
    pivot = rref_mat(i, i);
    rref_mat(i, :) = rref_mat(i, :) / pivot;
    row_ops{i} = sprintf('R%d -> (1/%d) R%d', i, pivot, i);
    elem_matrices{i} = eye(r);
    elem_matrices{i}(i, i) = 1/pivot;
    % Eliminate the entries below the pivot in the ith column
    for j = i+1:r
        if rref_mat(j, i) ~= 0
            factor = rref_mat(j, i);
            rref_mat(j, :) = rref_mat(j, :) - factor*rref_mat(i, :);
            row_ops{j-1} = sprintf('R%d -> R%d - (%d) R%d', j, j, factor, i);
            elem_matrices{j-1} = eye(r);
            elem_matrices{j-1}(j, i) = -factor;
        end
    end
end
% Perform back-substitution to obtain reduced row echelon form
for i = r:-1:2
    for j = i-1:-1:1
        if rref_mat(j, i) ~= 0
            factor = rref_mat(j, i);
            rref_mat(j, :) = rref_mat(j, :) - factor*rref_mat(i, :);
            row_ops{j} = sprintf('R%d -> R%d - (%d) R%d', j, j, factor, i);
            elem_matrices{j} = eye(r);
            elem_matrices{j}(j, i) = -factor;
        end
    end
end
% Round the entries to 6 decimal places to avoid rounding errors
rref_mat = round(rref_mat, 6);
% Print the row operations and corresponding elementary matrices
for i = 1:length(row_ops)
    fprintf('R%d: %s\n',i,row_ops{i});
end

end
% Compute RREF and row operations 