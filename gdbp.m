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
    
    r=GDBF_Decoder(y,H,100,1,1)
end


function decoded = GDBF_Decoder(received, H, max_iterations, theta1, theta2)
    % GDBF_Decoder - Gradient Descent Bit Flipping Decoder for LDPC codes
    % Inputs:
    % received       - received word (vector)
    % H              - parity-check matrix
    % max_iterations - maximum number of iterations
    % theta1         - initial threshold for multi-bit flipping
    % theta2         - threshold for escape process
    
    % Initialize decoded vector with hard decision
    decoded = sign(received);
    decoded(decoded == 0) = 1;
    mode_flag = 0; % Start in multi-bit mode
    theta = theta1; % Initial threshold

    fprintf('Initial decoded vector: %s\n', mat2str(decoded));

    for iteration = 1:max_iterations
        fprintf('Iteration %d:\n', iteration);
        
        % Compute objective function
        obj_function_value = compute_objective_function(decoded, received, H);
        fprintf('  Objective function value: %.4f\n', obj_function_value);
        
        if check_parity(decoded, H)
            fprintf('  Valid codeword found: %s\n', mat2str(decoded));
            return; % Decoding successful
        end

        if mode_flag == 0
            % Multi-bit mode
            fprintf('  Multi-bit mode:\n');
            for k = 1:length(decoded)
                inv_value = inversion_function(decoded, received, H, k);
                fprintf('    Inversion function value for bit %d: %.4f\n', k, inv_value);
                if inv_value < theta
                    decoded(k) = -decoded(k);
                    fprintf('    Bit %d flipped to %d\n', k, decoded(k));
                end
            end
            new_obj_function_value = compute_objective_function(decoded, received, H);
            fprintf('  New objective function value: %.4f\n', new_obj_function_value);
            if obj_function_value > new_obj_function_value
                mode_flag = 1; % Switch to single-bit mode
                fprintf('  Switching to single-bit mode.\n');
            end
        else
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
    end
    fprintf('Maximum iterations reached. Final decoded vector: %s\n', mat2str(decoded));
end

function sign_val = sign(x)
    % Custom sign function to handle zero values
    sign_val = ones(size(x));
    sign_val(x < 0) = -1;
end

function obj_val = compute_objective_function(decoded, received, H)
    % Compute the objective function
    correlation_term = sum(decoded .* received);
    parity_term = sum(arrayfun(@(i) prod(decoded(H(i, :) == 1)), 1:size(H, 1)));
    obj_val = correlation_term + parity_term;
end

function inv_val = inversion_function(decoded, received, H, k)
    % Compute the inversion function
    inv_val = decoded(k) * received(k) + sum(arrayfun(@(i) prod(decoded(H(i, :) == 1)), find(H(:, k) == 1)));
end

function is_valid = check_parity(decoded, H)
    % Check if all parity-check equations are satisfied
    is_valid = all(arrayfun(@(i) prod(decoded(H(i, :) == 1)) == 1, 1:size(H, 1)));
end
