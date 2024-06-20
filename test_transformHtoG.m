% [N, M, maxVNd, maxCNd, VNd, CNd, VNlink, CNlink, H] = f_readPCM_2024b('N15_K7_M8.txt');
[N, M, maxVNd, maxCNd, VNd, CNd, VNlink, CNlink, H] = f_readPCM_2024b('N96_K48_M48.txt');
% [N, M, maxVNd, maxCNd, VNd, CNd, VNlink, CNlink, H] = f_readPCM_2024b('N504_K252_M252.txt');
[G,sys]=transformHtoG(H)

%%
function [G,H_sys] = transformHtoG(H)
    % Ensure H is a binary matrix
    if ~all(H(:) == 0 | H(:) == 1)
        error('H matrix should be binary.');
    end
    
    % Get dimensions of H
    [m, n] = size(H)
    r=n-m

    rank(H)
    m

    
    % Perform Gaussian elimination to get H into systematic form [I | B]
    H_systematic = H;
    
    for i = 1:m
        % Make the diagonal element 1
        if H_systematic(i, i) == 0
            for j = i+1:m
                if H_systematic(j, i) == 1
                    H_systematic([i j], :) = H_systematic([j i], :);
                    break;
                end
            end
        end
        
        % If no swap was possible, continue to the next row
        if H_systematic(i, i) == 0
            continue;
        end
        
        % Make the elements below the diagonal 0
        for j = i+1:m
            if H_systematic(j, i) == 1
                H_systematic(j, :) = xor(H_systematic(j, :) , H_systematic(i, :));
            end
        end
    end
    
    % Make the elements above the diagonal 0
    for i = m:-1:1
        if H_systematic(i, i) == 1
            for j = 1:i-1
                if H_systematic(j, i) == 1
                    H_systematic(j, :) = mod(H_systematic(j, :)+ H_systematic(i, :),2);
                end
            end
        end
    end
    H_systematic(:, m+1:end)'
    H_sys=H_systematic
    % Extract P matrix
    G = [H_systematic(:, m+1:end)' eye(r)]
end