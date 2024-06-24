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
end
legend;
title(['Error-Rate Performance of (' num2str(N) ',' num2str(K) ') GDBF Multiple bit with escape (paper) Code; AWGN Channel']);
xlabel('Eb/N0 (dB)');
ylabel('Error Rates');
hold off;


