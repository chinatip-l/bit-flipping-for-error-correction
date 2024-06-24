%%
figure;
idx=1
theta=1:1:16
% N=96
% K=48
N=504
K=252
theta_sel=(theta(idx)-1)/10
semilogy(EbN0dB,combined_data.('BER').BER,'-','color',[0.2,0.2,0.2],'DisplayName','Uncoded BPSK BER');
% plot(EbN0dB,BPSK_BER_ana,'-','color',[0.2,0.2,0.2],'DisplayName','Uncoded BPSK BER');
title(['Performance of (' num2str(N) ',' num2str(K) ') IMWBF; theta=' num2str(theta_sel)]);

I_lim=2:2:12
t_algo='gdbf_single'
% t_algo='imwbf'
% t_algo='gdbf_multi'
% t_algo='gdbf_multi_escape_paper'
% t_algo='gdbf_multi_escape_improve'
% NK_sel='N96K48'
NK_sel='N504K252'

for i=1:numel(I_lim)
    colour=hsv2rgb([i/size(I_lim,2)*0.7,1,0.8])
    % hold on;
    % semilogy(EbN0dB,WBF_BER_sim(:,i),'-o','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);

    hold on;
    semilogy(EbN0dB,combined_data.(t_algo).(NK_sel).BER(idx,:,I_lim(i)/2),'-o','color',colour,'DisplayName',[' BER I=' num2str(I_lim(i))]);
    % semilogy(EbN0dB,combined_data.imwbf.(NK_set{2}).BER(8,:,I_lim(i)/2),'-^','color',colour,'DisplayName',[NK_set{2} ' BER I=' num2str(I_lim(i))]);
    % plot(EbN0dB,WBF_BER_sim(:,i),'-o','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);
    hold on;
    semilogy(EbN0dB,combined_data.(t_algo).(NK_sel).FER(idx,:,I_lim(i)/2),'--o','color',colour,'DisplayName',[' BLER I=' num2str(I_lim(i))]);
    % semilogy(EbN0dB,combined_data.imwbf.(NK_set{2}).FER(8,:,I_lim(i)/2),'--^','color',colour,'DisplayName',[NK_set{2} ' BLER I=' num2str(I_lim(i))]);
    % plot(EbN0dB,WBF_FER_sim(:,i),'-*','color',colour,'DisplayName',['BLER (WBF Algo) I=' num2str(I_lim(i))]);

end
% for i=1:numel(I_lim)
%     colour=hsv2rgb([i/size(I_lim,2)*0.7,1,0.8])
%     % hold on;
%     % semilogy(EbN0dB,WBF_BER_sim(:,i),'-o','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);
%
%     hold on;
%     % semilogy(EbN0dB,combined_data.imwbf.(NK_set{1}).BER(8,:,I_lim(i)/2),'-o','color',colour,'DisplayName',[NK_set{1} ' BER I=' num2str(I_lim(i))]);
%     semilogy(EbN0dB,combined_data.(algo).(NK_set{2}).BER(8,:,I_lim(i)/2),'-^','color',colour,'DisplayName',[NK_set{2} ' BER I=' num2str(I_lim(i))]);
%     % plot(EbN0dB,WBF_BER_sim(:,i),'-o','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);
%     hold on;
%     % semilogy(EbN0dB,combined_data.imwbf.(NK_set{1}).FER(8,:,I_lim(i)/2),'--o','color',colour,'DisplayName',[NK_set{1} ' BLER I=' num2str(I_lim(i))]);
%     semilogy(EbN0dB,combined_data.(algo).(NK_set{2}).FER(8,:,I_lim(i)/2),'--^','color',colour,'DisplayName',[NK_set{2} ' BLER I=' num2str(I_lim(i))]);
%     % plot(EbN0dB,WBF_FER_sim(:,i),'-*','color',colour,'DisplayName',['BLER (WBF Algo) I=' num2str(I_lim(i))]);
%
% end
legend;
xlabel('Eb/N0 (dB)');
ylabel('Error Rates');
hold off;


