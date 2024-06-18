G=[1 1 0 1 0 0;
   0 1 1 0 1 0;
   1 0 1 0 0 1]
H=[1 0 0 1 0 1;
   0 1 0 1 1 0;
   0 0 1 0 1 1]
c=[0 1 1]
v=(1-mod(c*G,2)*2)

EbN0dB=[0:1:7]

sigma = sqrt(1./ (2*R*(10.^(EbN0dB/10))));
for idx=1:size(EbN0dB,2)
    EbN0dB(idx)

    y = v + sigma(idx) * randn(1,size(v,2))
    
    % y=[ -0.2 0.7 0.1 -1 0.5 0.3]
    % abs_y=abs(y)
    v
    % x_guess=sign(y)
    
    % syn_n=(x_guess<0)*H'
    
    r=GDBF_SingleBitFlip_Decoder(y,H,100)
end


function decoded = GDBF_SingleBitFlip_Decoder(received, H, max_iterations)
    % GDBF_SingleBitFlip_Decoder - Gradient Descent Bit Flipping Decoder for LDPC codes
    % Single-bit flip version
    % Inputs:
    % received       - received word (vector)
    % H              - parity-check matrix
    % max_iterations - maximum number of iterations
    
    % Initialize decoded vector with hard decision
    decoded = sign(received);
    decoded(decoded == 0) = 1;
    
    fprintf('Initial decoded vector: %s\n', mat2str(decoded));

    for iteration = 1:max_iterations
        fprintf('Iteration %d:\n', iteration);
        
        % Compute objective function
        [obj_function_value, correlation_term, parity_term] = compute_objective_function(decoded, received, H);
        fprintf('  Objective function value: %.4f\n', obj_function_value);
        fprintf('  Correlation term: %.4f\n', correlation_term);
        fprintf('  Parity term: %.4f\n', parity_term);
        
        if check_parity(decoded, H)
            fprintf('  Valid codeword found: %s\n', mat2str(decoded));
            return; % Decoding successful
        end

        % Single-bit mode
        fprintf('  Single-bit mode:\n');
        min_inversion_value = Inf;
        flip_position = NaN;
        for k = 1:length(decoded)
            inv_value = inversion_function(decoded, received, H, k);
            fprintf('    Inversion function value for bit %d: %.4f\n', k, inv_value);
            if inv_value < min_inversion_value
                min_inversion_value = inv_value;
                flip_position = k;
            end
        end
        if ~isnan(flip_position)
            decoded(flip_position) = -decoded(flip_position);
            fprintf('    Bit %d flipped to %d\n', flip_position, decoded(flip_position));
        end
    end
    fprintf('Maximum iterations reached. Final decoded vector: %s\n', mat2str(decoded));
end

function sign_val = sign(x)
    % Custom sign function to handle zero values
    sign_val = ones(size(x));
    sign_val(x < 0) = -1;
end

function [obj_val, correlation_term, parity_term] = compute_objective_function(decoded, received, H)
    % Compute the objective function
    % Correlation term: sum(decoded .* received)
    correlation_term = sum(decoded .* received);
    
    % Parity term: sum(arrayfun(@(i) prod(decoded(H(i, :) == 1)), 1:size(H, 1)))
    parity_term = 0;
    for i = 1:size(H, 1)
        parity_product = 1;
        for j = find(H(i, :))
            parity_product = parity_product * decoded(j);
        end
        parity_term = parity_term + parity_product;
    end
    
    % Objective function: correlation_term + parity_term
    obj_val = correlation_term + parity_term;
end

function inv_val = inversion_function(decoded, received, H, k)
    % Compute the inversion function
    % inv_val = decoded(k) * received(k) + sum(arrayfun(@(i) prod(decoded(H(i, :) == 1)), find(H(:, k) == 1)))
    inv_val = decoded(k) * received(k);
    for i = find(H(:, k))'
        parity_product = 1;
        for j = find(H(i, :))
            parity_product = parity_product * decoded(j);
        end
        inv_val = inv_val + parity_product;
    end
end

function is_valid = check_parity(decoded, H)
    % Check if all parity-check equations are satisfied
    is_valid = all(arrayfun(@(i) prod(decoded(H(i, :) == 1)) == 1, 1:size(H, 1)));
end
