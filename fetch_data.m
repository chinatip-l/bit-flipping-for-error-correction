% Define the folder where your .mat files are located
folder_path = './'; % Replace with your actual folder path
iter_cnt=2:2:50
% Define the algorithm names, bit pairs, and theta values
algos = {'imwbf','gdbf_multi','gdbf_multi_escape_paper','gdbf_multi_escape_improve' };  % Replace with your algorithm names
nk_values = [96 48 ;
    504 252];         % Replace with your k values
theta_values = 0:-0.1:-1.5;     % Replace with your theta values

EbN0dB=0:1:7
% Initialize a structure to hold the combined data
combined_data = struct();
NK_set={}
for i=1:length(nk_values)
    nk=nk_values(i,:)
    n = nk(1);
    k = nk(2);
    NK = sprintf('N%dK%d', n, k);
    NK_set=[NK_set;NK]
end

% Loop through each algorithm, n, k, and theta combination
for algo_idx = 1:length(algos)
    for i=1:length(nk_values)
        nk_name=NK_set{i}
        nk=nk_values(i,:)
        n = nk(1);
        k = nk(2);

        % Initialize matrices to store WBF_BER_sim and WBF_FER_sim
        WBF_BER_sim_all = zeros( length(theta_values),8,25);
        WBF_FER_sim_all = zeros( length(theta_values),8,25);
        % BER_ANA =

        for theta_idx = 1:length(theta_values)
            theta = theta_values(theta_idx);
            if strcmp(algos{algo_idx},'imwbf') && theta~=0
                theta=-theta
            end
            % Construct the file name pattern
            file_pattern = sprintf('%s_N%dK%d_theta_%.2f.mat', algos{algo_idx}, n, k, theta);

            % Check if the file exists
            file_path = fullfile(folder_path, file_pattern);
            if exist(file_path, 'file') == 2
                % Load the data from the file
                data = load(file_path);

                % Extract WBF_BER_sim and WBF_FER_sim
                WBF_BER_sim = data.WBF_BER_sim;
                WBF_FER_sim = data.WBF_FER_sim;
                BER_ANA_sim = data.BPSK_BER_ana;

                % Store in the combined matrices
                WBF_BER_sim_all( theta_idx,:,:) = WBF_BER_sim(:,:);
                WBF_FER_sim_all( theta_idx,:,:) = WBF_FER_sim(:,:);
            else
                % Handle case where file doesn't exist (optional)
                disp(['File not found: ' file_pattern]);
            end
        end

        % Store the combined matrices in the structure


        combined_data.(algos{algo_idx}).(nk_name).BER = WBF_BER_sim_all;
        combined_data.(algos{algo_idx}).(nk_name).FER = WBF_FER_sim_all;
        combined_data.('BER').BER = BER_ANA_sim;
    end
end

% figure;
%     semilogy(EbN0dB,combined_data.('BER').BER,'-','color',[0.2,0.2,0.2],'DisplayName','Uncoded BPSK BER');
%     % plot(EbN0dB,BPSK_BER_ana,'-','color',[0.2,0.2,0.2],'DisplayName','Uncoded BPSK BER');
% I_lim=10:10:50
% for i=1:numel(I_lim)
%     colour=hsv2rgb([i/size(I_lim,2)*0.7,1,0.8])
%     % hold on;
%         % semilogy(EbN0dB,WBF_BER_sim(:,i),'-o','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);
%     hold on;
%         semilogy(EbN0dB,combined_data.(algos{algo_idx}).(NK).FER,'-o','color',colour,'DisplayName',['BER I=' num2str(I_lim(i))]);
%         % plot(EbN0dB,WBF_BER_sim(:,i),'-o','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);
%     hold on;
%         semilogy(EbN0dB,WBF_FER_sim(:,i),'--*','color',colour,'DisplayName',['BLER I=' num2str(I_lim(i))]);
%         % plot(EbN0dB,WBF_FER_sim(:,i),'-*','color',colour,'DisplayName',['BLER (WBF Algo) I=' num2str(I_lim(i))]);
%
% end
% hold off;


% Now you have combined_data structure containing all the data
% Example usage:
% To access WBF_BER_sim for algo1 with N=96, K=48:
% combined_data.algo1_N96K48.WBF_BER_sim
% To access WBF_FER_sim for algo2 with N=504, K=252:
% combined_data.algo2_N504K252.WBF_FER_sim

% Optional: Save combined_data to a .mat file
save('combined_data.mat', 'combined_data');
%%
figure;
semilogy(EbN0dB,combined_data.('BER').BER,'-','color',[0.2,0.2,0.2],'DisplayName','Uncoded BPSK BER');
% plot(EbN0dB,BPSK_BER_ana,'-','color',[0.2,0.2,0.2],'DisplayName','Uncoded BPSK BER');
title("IMWBF ")
I_lim=10:10:50
for k=1:length(algos)
    algo=algos{k}
for i=1:numel(I_lim)
    colour=hsv2rgb([i/size(I_lim,2)*0.7,1,0.8])
    % hold on;
    % semilogy(EbN0dB,WBF_BER_sim(:,i),'-o','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);

    hold on;
    semilogy(EbN0dB,combined_data.(algo).(NK_set{1}).BER(8,:,I_lim(i)/2),'-o','color',colour,'DisplayName',[NK_set{1} ' BER I=' num2str(I_lim(i))]);
    % semilogy(EbN0dB,combined_data.imwbf.(NK_set{2}).BER(8,:,I_lim(i)/2),'-^','color',colour,'DisplayName',[NK_set{2} ' BER I=' num2str(I_lim(i))]);
    % plot(EbN0dB,WBF_BER_sim(:,i),'-o','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);
    hold on;
    semilogy(EbN0dB,combined_data.(algo).(NK_set{1}).FER(8,:,I_lim(i)/2),'--o','color',colour,'DisplayName',[NK_set{1} ' BLER I=' num2str(I_lim(i))]);
    % semilogy(EbN0dB,combined_data.imwbf.(NK_set{2}).FER(8,:,I_lim(i)/2),'--^','color',colour,'DisplayName',[NK_set{2} ' BLER I=' num2str(I_lim(i))]);
    % plot(EbN0dB,WBF_FER_sim(:,i),'-*','color',colour,'DisplayName',['BLER (WBF Algo) I=' num2str(I_lim(i))]);

end
for i=1:numel(I_lim)
    colour=hsv2rgb([i/size(I_lim,2)*0.7,1,0.8])
    % hold on;
    % semilogy(EbN0dB,WBF_BER_sim(:,i),'-o','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);

    hold on;
    % semilogy(EbN0dB,combined_data.imwbf.(NK_set{1}).BER(8,:,I_lim(i)/2),'-o','color',colour,'DisplayName',[NK_set{1} ' BER I=' num2str(I_lim(i))]);
    semilogy(EbN0dB,combined_data.(algo).(NK_set{2}).BER(8,:,I_lim(i)/2),'-^','color',colour,'DisplayName',[NK_set{2} ' BER I=' num2str(I_lim(i))]);
    % plot(EbN0dB,WBF_BER_sim(:,i),'-o','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);
    hold on;
    % semilogy(EbN0dB,combined_data.imwbf.(NK_set{1}).FER(8,:,I_lim(i)/2),'--o','color',colour,'DisplayName',[NK_set{1} ' BLER I=' num2str(I_lim(i))]);
    semilogy(EbN0dB,combined_data.(algo).(NK_set{2}).FER(8,:,I_lim(i)/2),'--^','color',colour,'DisplayName',[NK_set{2} ' BLER I=' num2str(I_lim(i))]);
    % plot(EbN0dB,WBF_FER_sim(:,i),'-*','color',colour,'DisplayName',['BLER (WBF Algo) I=' num2str(I_lim(i))]);

end
end
hold off;


