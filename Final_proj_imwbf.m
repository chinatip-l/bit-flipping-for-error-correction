clear;

%%
% Load PCM and related parameters
% [N, M, maxVNd, maxCNd, VNd, CNd, VNlink, CNlink, H] = f_readPCM_2024b('N15_K7_M8.txt');
% [N, M, maxVNd, maxCNd, VNd, CNd, VNlink, CNlink, H] = f_readPCM_2024b('N96_K48_M48.txt');
[N, M, maxVNd, maxCNd, VNd, CNd, VNlink, CNlink, H] = f_readPCM_2024b('N504_K252_M252.txt');
% [N_bf, M_bf, maxVNd_bf, maxCNd_bf, VNd_bf, CNd_bf, VNlink_bf, CNlink_bf, H_bf] = f_readPCM_2024b('N15_K7_M15.txt');
K = N-M;
% Generator matrix
% (15,7)
Ht=H'
G_1=[1 0 0 0 1 0 1 1 1 0 0 0 0 0 0;
   1 1 0 0 1 1 1 0 0 1 0 0 0 0 0;
   0 1 1 0 0 1 1 1 0 0 1 0 0 0 0;
   1 0 1 1 1 0 0 0 0 0 0 1 0 0 0;
   0 1 0 1 1 1 0 0 0 0 0 0 1 0 0;
   0 0 1 0 1 1 1 0 0 0 0 0 0 1 0;
   0 0 0 1 0 1 1 1 0 0 0 0 0 0 1]
G=transformHtoG(H)
% Coderate
R = K/N; 
I_lim=(1:2:50)
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
        abs_y=abs(Rx_wbf);

        % original weight
        ck_weight=zeros(1,M);
        for i=1:M
            cur_check_en_bit=H(i,:);
            cur_check_weight=abs_y.*cur_check_en_bit;
            cur_check_weight(cur_check_weight == 0) = inf;
            [minw,minI]=min(cur_check_weight);
            ck_weight(i)=minw;
        end

        ck_weight_im=zeros(M,maxCNd);
        for i=1:M
            for j=1:CNd(i)
                cur_check_en_bit_excld=H(i,:);
                cur_check_en_bit_excld(CNlink(i,j))=0;
                cur_check_weight=abs_y.*cur_check_en_bit_excld;
                cur_check_weight(cur_check_weight == 0) = inf;
                [minw,minI]=min(cur_check_weight);
                ck_weight_im(i,j)=minw;
            end
        end

        ck_weight_im=ck_weight_im';

        

        guess_rx_wbf=Rx_wbf<0;
        res_wbf=guess_rx_wbf;

        for i=1:I_max
            % En=zeros(1,N)
            % calculate En from guess rx
          
            Syn_wbf = mod(guess_rx_wbf*H.',2);

            if sum(Syn_wbf)>0
                
                % original WBF
                % WBF=-((1-2*Syn_wbf).*ck_weight)*H;
                % IMWBF
                % WBF=-(((1-2*Syn_wbf).*ck_weight_im*H)-(alpha*abs_y));

                % this one is the modified one
                % shold be normalised to 
                % rz=((1-2*Syn_wbf)'.*ck_weight_im)
                rz=(1-2*Syn_wbf).*ck_weight_im;
                rz=sum(rz,1);

                WBF=-((rz*H)-(alpha*abs_y));

                
                [max_v,max_i]=max(WBF);
                guess_rx_wbf(max_i)=~guess_rx_wbf(max_i);
                res_wbf=guess_rx_wbf;
                
                col=find(I_lim==i);
                % not reach any lim
                if isempty(col)

                else
                % reach some limit
                    res_wbf_tmp=guess_rx_wbf;
                    wbf_Num_error_bit_temp = sum(Tx_codeword ~= res_wbf_tmp);
                    if wbf_Num_error_bit_temp > 0
                        wbf_ber(idx,col) = wbf_ber(idx,col) + wbf_Num_error_bit_temp
                        wbf_bler(idx,col) = wbf_bler(idx,col) + 1;
                    end
                end
            else
                res_wbf=guess_rx_wbf;
                % fprintf("Got %d",res_wbf)
                % wbf_codeword_error_num(idx,col)
                res_wbf_tmp=guess_rx_wbf;
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
    colour=hsv2rgb([i/size(I_lim,2),1,0.8])
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
title(['Error-Rate Performance of (' num2str(N) ',' num2str(K) ') IMWBF Code; AWGN Channel']);
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