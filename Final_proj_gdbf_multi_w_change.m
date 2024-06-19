clear;

%%
% Load PCM and related parameters
% [N, M, maxVNd, maxCNd, VNd, CNd, VNlink, CNlink, H] = f_readPCM_2024b('k12.txt');
% [N, M, maxVNd, maxCNd, VNd, CNd, VNlink, CNlink, H] = f_readPCM_2024b('ldpc_matrix.alist');
% [N, M, maxVNd, maxCNd, VNd, CNd, VNlink, CNlink, H] = f_readPCM_2024b('N15_K7_M8.txt');
% [N, M, maxVNd, maxCNd, VNd, CNd, VNlink, CNlink, H] = f_readPCM_2024b('N96_K48_M48.txt');
[N, M, maxVNd, maxCNd, VNd, CNd, VNlink, CNlink, H] = f_readPCM_2024b('N504_K252_M252.txt');
% [N_bf, M_bf, maxVNd_bf, maxCNd_bf, VNd_bf, CNd_bf, VNlink_bf, CNlink_bf, H_bf] = f_readPCM_2024b('N15_K7_M15.txt');
K = N-M;
% Generator matrix
% (15,7)
Ht=H';
G_1=[1 0 0 0 1 0 1 1 1 0 0 0 0 0 0;
   1 1 0 0 1 1 1 0 0 1 0 0 0 0 0;
   0 1 1 0 0 1 1 1 0 0 1 0 0 0 0;
   1 0 1 1 1 0 0 0 0 0 0 1 0 0 0;
   0 1 0 1 1 1 0 0 0 0 0 0 1 0 0;
   0 0 1 0 1 1 1 0 0 0 0 0 0 1 0;
   0 0 0 1 0 1 1 1 0 0 0 0 0 0 1]
G=transformHtoG(H);
% b=[ binary_null_space(H) eye(K)]
% Coderate
R = K/N; 
I_lim=(2:2:100)
I_max=max(I_lim)

%%
EbN0dB = [0:1:7];
sigma = sqrt(1./ (2*R*(10.^(EbN0dB/10))));
minimum_codeword_error_num = 500;
soft_minimum_codeword_error_num = 500;

coded_bit_error_num = zeros(size(EbN0dB));
codeword_error_num = zeros(size(EbN0dB));

hbf_coded_bit_error_num= zeros(size(EbN0dB));
hbf_codeword_error_num = zeros(size(EbN0dB));

wbf_ber= zeros(size(EbN0dB,2),size(I_lim,2));
wbf_ber(:,:)=-1;
wbf_bler = zeros(size(EbN0dB,2),size(I_lim,2));
wbf_bler(:,:)=-1;

total_codeword_num = zeros(size(EbN0dB));
alpha=0.01
theta = -0.5; % Initial threshold
theta_increment=0.1

%%
% allmsgs_dec = [0:2^K-1];
% allmsgs_bin = dec2bin(allmsgs_dec,K);
% allmsgs_str = str2num(allmsgs_bin(:));
% msg_set = reshape(allmsgs_str,2^K,[]);
% codeword_set= mod(msg_set*G,2);

%%

%%
tic

for idx=1:length(EbN0dB)
    while any(wbf_ber(idx,:) < minimum_codeword_error_num)
        total_codeword_num(idx) = total_codeword_num(idx) + 1;
        % r=wbf_codeword_error_num(idx,:)
        % K-bit source data generation
        Tx_data = randi([0 1],1,K); 

        % Encoding
        Tx_codeword = mod(Tx_data * G,2); 

        % BPSK modulation
        Tx_codeword_BPSK = 1 - 2 * Tx_codeword; 

        % AWGN channel
        Rx_received = Tx_codeword_BPSK + sigma(idx) * randn(1,N); %AWGN channel

        % BPSK demodulation
        Rx_BPSK_demod = (Rx_received < 0); 
        
        Rx_hbf=Rx_BPSK_demod;

        Rx_wbf=Rx_received;
       

        % Computing syndrome
        % Syn = mod(Rx_BPSK_demod*H.',2);
        % Soft_syn=mod(Soft_Rx*H.',2);

        % Error correction BPSK
        % Syn_bit_sum = sum(Syn);
        % if Syn_bit_sum~=0
        %     Hamming_distance = vecnorm(ones(2^K,1)*Rx_BPSK_demod-codeword_set,2,2);
        %     [minHammD,minI] = min(Hamming_distance);
        %     decoded_codeword = codeword_set(minI,:);
        % else
        %     decoded_codeword = Rx_BPSK_demod;
        % end

        % IMWBF algo
        % calculate imwbf

        guess_rx=sign(Rx_wbf);
        mode_flag = 0; % Start in multi-bit mode
        
        prev_flip_position=0


        for i=1:I_max
            % En=zeros(1,N)
            % calculate En from guess rx
          
            Syn_wbf = mod((guess_rx<0)*H.',2);

            fprintf('Iteration %d:\n', i);
            
            % Compute objective function
            [obj_function_value, correlation_term, parity_term] = compute_objective_function(guess_rx, Rx_wbf, H);
            fprintf('  Objective function value: %.4f\n', obj_function_value);
            fprintf('  Correlation term: %.4f\n', correlation_term);
            fprintf('  Parity term: %.4f\n', parity_term);
            


            if sum(Syn_wbf)>0
                
                
                if mode_flag == 0
                    % Multi-bit mode
                    fprintf('  Multi-bit mode:\n');
                    % for k = 1:length(decoded)
                    %     inv_value = inversion_function(decoded, received, H, k);
                    %     fprintf('    Inversion function value for bit %d: %.4f\n', k, inv_value);
                    %     if inv_value < theta
                    %         decoded(k) = -decoded(k);
                    %         fprintf('    Bit %d flipped to %d\n', k, decoded(k));
                    %     end
                    % end
                    inv_vect=arrayfun(@(k) inversion_function(guess_rx,Rx_wbf,H,k),1:numel(guess_rx));
                    
                    below_thr=-(inv_vect<theta);
                    below_thr(below_thr==0)=1;
                    guess_rx=guess_rx.*below_thr;

                    for z = 1:length(find(below_thr==-1))
                        fprintf('    Bit %d flipped to %d\n', z, guess_rx(z));
                    end


                    new_obj_function_value = compute_objective_function(guess_rx, Rx_wbf, H);
                    fprintf('  New objective function value: %.4f\n', new_obj_function_value);
                    if obj_function_value >= new_obj_function_value && i>1
                        mode_flag = 1; % Switch to single-bit mode
                        fprintf('  Switching to single-bit mode.\n');
                    end
                else
                    % fprintf('  Single-bit mode:\n');
                    inv_vect=arrayfun(@(k) inversion_function(guess_rx,Rx_wbf,H,k),1:numel(guess_rx));
                    [min_val,min_idx]=min(inv_vect);
                    
                    % guess_rx(min_idx) = -guess_rx(min_idx);
                    

                    % % check if better
                    % new_obj_function_value = compute_objective_function(guess_rx, Rx_wbf, H);
                    % fprintf('  New objective function value: %.4f\n', new_obj_function_value);
                    % if obj_function_value >= new_obj_function_value && i>1
                    %     mode_flag = 1; % Switch to single-bit mode
                    %     fprintf('  Switching to single-bit mode.\n');
                    % end
                    
                    % Check for oscillation
                    if min_idx == prev_flip_position
                        oscillation_count = oscillation_count + 1;
                    else
                        oscillation_count = 0;
                        step_size=1;
                    end
                    
                    % Escape mechanism
                    if oscillation_count > 2  % Adjust this threshold based on empirical observations
                        % theta = theta + theta_increment;
                        step_size = step_size + 1;
                        [~, sorted_indices] = sort(guess_rx);
                        flip_positions = sorted_indices(1:step_size);
                        % flip_positions = find(inv_vect < theta);
                        % if isempty(flip_positions)
                        %     flip_positions = min_idx;  % Ensure at least one bit is flipped
                        % end
                        oscillation_count = 0;  % Reset oscillation counter after escape attempt
                    else
                        flip_positions = min_idx;
                    end
                    
                    
                    % Flip the selected bits
                    flip_positions
                    guess_rx(flip_positions) = -guess_rx(flip_positions);
                    for z = 1:length(flip_positions)
                        fprintf('    Bit %d flipped to %d\n', flip_positions(z), guess_rx(flip_positions(z)));
                    end
                    
                    % Track the previous flip position
                    prev_flip_position = flip_positions;



                    
                    res_wbf_tmp=guess_rx;
                end
                
                col=find(I_lim==i);
                % not reach any lim
                if isempty(col)

                else
                % reach some limit
                    res_wbf_tmp=guess_rx<0;
                    wbf_Num_error_bit_temp = sum(Tx_codeword ~= res_wbf_tmp);
                    if wbf_Num_error_bit_temp > 0
                        wbf_ber(idx,col) = wbf_ber(idx,col) + wbf_Num_error_bit_temp
                        wbf_bler(idx,col) = wbf_bler(idx,col) + 1;
                    end
                end
            else
                
                % wbf_codeword_error_num(idx,col)
                res_wbf_tmp=guess_rx<0;
                fprintf("Got %d\n",i)
                res_wbf_tmp
                Tx_codeword
                wbf_Num_error_bit_temp = sum(Tx_codeword ~= res_wbf_tmp);
                if wbf_Num_error_bit_temp > 0
                    for col=find(I_lim>=i)
                        wbf_ber(idx,col) = wbf_ber(idx,col) + wbf_Num_error_bit_temp;
                        wbf_bler(idx,col) = wbf_bler(idx,col) + 1;
                    end
                end
                break
            end
        end
        zzzzzz=1

        % Num_error_bit_temp = sum(Tx_codeword ~= decoded_codeword);
        % if Num_error_bit_temp > 0
        %     codeword_error_num(idx) = codeword_error_num(idx) + Num_error_bit_temp;
        %     coded_bit_error_num(idx) = coded_bit_error_num(idx) + 1;
        % end

        % wbf_Num_error_bit_temp = sum(Tx_codeword ~= res_wbf);
        % if wbf_Num_error_bit_temp > 0
        %     wbf_codeword_error_num(idx) = wbf_coded_bit_error_num(idx) + wbf_Num_error_bit_temp
        %     wbf_coded_bit_error_num(idx) = wbf_coded_bit_error_num(idx) + 1;
        % end

    end
end

toc

%% Compute simulated error rates
% hard_BER_sim = codeword_error_num./(total_codeword_num*N);
% hard_FER_sim = coded_bit_error_num./(total_codeword_num); % FER is aka BLER

% WBF_BER_sim = wbf_ber./repmat((total_codeword_num*N)',1,size(I_lim,2))
WBF_BER_sim = wbf_ber./repmat((total_codeword_num*N)',1,size(I_lim,2))
WBF_FER_sim = wbf_bler./repmat((total_codeword_num)',1,size(I_lim,2)) % FER is aka BLER

%% BER of uncoded BPSK for comparison
uncodedSNR_EbN0dB = EbN0dB;
uncodedSNR_EbN0 = 10.^(uncodedSNR_EbN0dB/10);
BPSK_BER_ana = 0.5*erfc(sqrt(uncodedSNR_EbN0)) ;

%%
figure;
    semilogy(EbN0dB,BPSK_BER_ana,'-','color',[0.2,0.2,0.2],'DisplayName','Uncoded BPSK BER');
    % plot(EbN0dB,BPSK_BER_ana,'-','color',[0.2,0.2,0.2],'DisplayName','Uncoded BPSK BER');
hold on;
    % semilogy(EbN0dB,hard_BER_sim,'--*','color',[0.2,0.2,0.8],'DisplayName','BER (Hard Decision)');

hold on;
    % semilogy(EbN0dB,hard_FER_sim,'-*','color',[0.2,0.2,0.8],'DisplayName','BLER (Hard Decision)');

for i=1:size(I_lim,2)
    colour=hsv2rgb([i/size(I_lim,2)*0.7,1,0.8])
    % hold on;
        % semilogy(EbN0dB,WBF_BER_sim(:,i),'-o','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);
    hold on;
        semilogy(EbN0dB,WBF_BER_sim(:,i),'-o','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);
        % plot(EbN0dB,WBF_BER_sim(:,i),'-o','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);
    hold on;
        semilogy(EbN0dB,WBF_FER_sim(:,i),'--*','color',colour,'DisplayName',['BLER (WBF Algo) I=' num2str(I_lim(i))]);
        % plot(EbN0dB,WBF_FER_sim(:,i),'-*','color',colour,'DisplayName',['BLER (WBF Algo) I=' num2str(I_lim(i))]);

end

hold off;
legend;
title(['Error-Rate Performance of (' num2str(N) ',' num2str(K) ') GDBF Code; AWGN Channel']);
xlabel('Eb/N0 (dB)');
ylabel('Error Rates');

toc

%%
function G = transformHtoG(H)
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
                    H_systematic(j, :) = mod(H_systematic(j, :) + H_systematic(i, :), 2);
                end
            end
        end
    end
    
    % Extract P matrix
    G = [H_systematic(:, m+1:end)' eye(r)]
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

function B = binary_null_space(H)
    % Ensure the matrix H is binary (mod 2)
    H = mod(H, 2);
    
    % Perform Gaussian elimination modulo 2
    [rows, cols] = size(H);
    augmented_matrix = [H eye(rows)]; % Augment with identity matrix

    % Gaussian elimination
    for i = 1:rows
        % Find pivot
        for j = i:rows
            if augmented_matrix(j, i) == 1
                % Swap rows
                temp = augmented_matrix(i, :);
                augmented_matrix(i, :) = augmented_matrix(j, :);
                augmented_matrix(j, :) = temp;
                break;
            end
        end
        
        % Make all rows below this one 0 in the current column
        for k = i+1:rows
            if augmented_matrix(k, i) == 1
                augmented_matrix(k, :) = mod(augmented_matrix(k, :) + augmented_matrix(i, :), 2);
            end
        end
    end

    % Back substitution to reduce to row echelon form
    for i = rows:-1:1
        % Make all rows above this one 0 in the current column
        for k = i-1:-1:1
            if augmented_matrix(k, i) == 1
                augmented_matrix(k, :) = mod(augmented_matrix(k, :) + augmented_matrix(i, :), 2);
            end
        end
    end

    % Extract the null space basis vectors
    B = augmented_matrix(:, cols+1:end);

    % Display result
    disp('The null space basis vectors in binary form are:');
    disp(B);
end