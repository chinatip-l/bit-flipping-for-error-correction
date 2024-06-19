function write_alist(H, filename)
    % Ensure the matrix H is binary (mod 2)
    H = mod(H, 2);

    % Get the size of the matrix
    [m, n] = size(H);

    % Calculate column weights
    col_weights = sum(H, 1);

    % Calculate row weights
    row_weights = sum(H, 2);

    % Maximum column weight
    maxcolwt = max(col_weights);

    % Maximum row weight
    maxrowwt = max(row_weights);

    % Open the file to write
    fileID = fopen(filename, 'w');

    % Write the header information
    fprintf(fileID, '%d %d\n', n, m);
    fprintf(fileID, '%d %d\n', maxcolwt, maxrowwt);

    % Write the column weights
    fprintf(fileID, '%d ', col_weights);
    fprintf(fileID, '\n');

    % Write the row weights
    fprintf(fileID, '%d ', row_weights);
    fprintf(fileID, '\n');

    % Write the column indices
    for col = 1:n
        nz_indices = find(H(:, col))';
        fprintf(fileID, '%d ', nz_indices);
        fprintf(fileID, '\n');
    end

    % Write the row indices
    for row = 1:m
        nz_indices = find(H(row, :));
        fprintf(fileID, '%d ', nz_indices);
        fprintf(fileID, '\n');
    end

    % Close the file
    fclose(fileID);
end

% Parameters
K = 252;  % Number of information bits
N = 504;  % Total number of bits (assuming rate 1/2)
M = N - K;  % Number of parity-check bits

% Generate a random LDPC parity-check matrix H (example construction)
% For real applications, use a structured construction method like PEG or QC.
H = sparse(M, N);
for i = 1:M
    col_indices = randperm(N, 3);  % Ensure low density, e.g., 3 ones per row
    H(i, col_indices) = 1;
end

% Write the AList file
write_alist(H, 'ldpc_matrix.alist');
