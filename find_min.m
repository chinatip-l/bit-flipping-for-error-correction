A=combined_data.imwbf.N96K48.BER
A=single_96.WBF_BER_sim
[min_value, min_idx] = min(A(A~=0));
[idx1, idx2, idx3] = ind2sub(size(A), min_idx);

% Display results
fprintf('Minimum value: %ld\n', min_value);
fprintf('Indices: (%d, %d, %d)\n', idx1, idx2, idx3);

%%
%%
EbN0dB = [0:1:7];
uncodedSNR_EbN0dB = EbN0dB;
uncodedSNR_EbN0 = 10.^(uncodedSNR_EbN0dB/10);
BPSK_BER_ana = 0.5*erfc(sqrt(uncodedSNR_EbN0)) ;
load("combined_data.mat")
single_504=load("gdbf_single_N504K252.mat")
single_96=load("gdbf_single_N96K48.mat")
figure;
idx=4
theta=1:1:16
theta_sel=-(theta(idx)-1)/10
semilogy(EbN0dB,BPSK_BER_ana,'-','color',[0.2,0.2,0.2],'DisplayName','Uncoded BPSK BER');
% plot(EbN0dB,BPSK_BER_ana,'-','color',[0.2,0.2,0.2],'DisplayName','Uncoded BPSK BER');
title(['Performance Comparison' ]);

I_lim=10:10:50
% t_algo='gdbf_single'
% t_algo='gdbf_multi'
% t_algo='gdbf_multi_escape_paper'
t_algo='gdbf_multi_escape_improve'
NK_9='N96K48'
NK_5='N504K252'
hold on;
tmp_algo="imwbf"
tmp_sz=NK_5
theta=1
i_lim=10
colour=hsv2rgb([1/10*0.9,1,0.7])
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).BER(theta,:,i_lim),'-^','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).FER(theta,:,i_lim),'--*','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);
tmp_algo="imwbf"
tmp_sz=NK_9
theta=1
i_lim=5
colour=hsv2rgb([2/10*0.9,1,0.7])
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).BER(theta,:,i_lim),'-^','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).FER(theta,:,i_lim),'--*','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);

tmp_algo="gdbf_single"
tmp_sz=NK_5
theta=8
i_lim=10
colour=hsv2rgb([3/10*0.9,1,0.7])
semilogy(EbN0dB,single_504.WBF_BER_sim(:,i_lim),'-^','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);
semilogy(EbN0dB,single_504.WBF_FER_sim(:,i_lim),'--*','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);

tmp_algo="gdbf_single"
tmp_sz=NK_9
theta=1
i_lim=5
colour=hsv2rgb([4/10*0.9,1,0.7])
semilogy(EbN0dB,single_96.WBF_BER_sim(:,i_lim),'-^','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);
semilogy(EbN0dB,single_96.WBF_FER_sim(:,i_lim),'--*','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);


tmp_algo="gdbf_multi"
tmp_sz=NK_5
theta=6
i_lim=6
colour=hsv2rgb([5/10*0.9,1,0.7])
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).BER(theta,:,i_lim),'-^','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).FER(theta,:,i_lim),'--*','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);

tmp_algo="gdbf_multi"
tmp_sz=NK_9
theta=7
i_lim=6
colour=hsv2rgb([6/10*0.9,1,0.7])
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).BER(theta,:,i_lim),'-^','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).FER(theta,:,i_lim),'--*','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);

tmp_algo="gdbf_multi_escape_paper"
tmp_sz=NK_5
theta=5
i_lim=12
colour=hsv2rgb([7/10*0.9,1,0.7])
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).BER(theta,:,i_lim),'-^','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).FER(theta,:,i_lim),'--*','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);

tmp_algo="gdbf_multi_escape_paper"
tmp_sz=NK_9
theta=9
i_lim=4
colour=hsv2rgb([8/10*0.9,1,0.7])
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).BER(theta,:,i_lim),'-^','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).FER(theta,:,i_lim),'--*','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);

tmp_algo="gdbf_multi_escape_improve"
tmp_sz=NK_5
theta=4
i_lim=18
colour=hsv2rgb([9/10*0.9,1,0.7])
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).BER(theta,:,i_lim),'-^','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).FER(theta,:,i_lim),'--*','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);

tmp_algo="gdbf_multi_escape_improve"
tmp_sz=NK_9
theta=8
i_lim=14
colour=hsv2rgb([10/10*0.9,1,0.7])
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).BER(theta,:,i_lim),'-^','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);
semilogy(EbN0dB,combined_data.(tmp_algo).(tmp_sz).FER(theta,:,i_lim),'--*','color',colour,'DisplayName',[sprintf("%s %s",tmp_algo,tmp_sz)]);

% for i=1:numel(I_lim)
%     colour=hsv2rgb([i/size(I_lim,2)*0.9,1,0.7])
%     % hold on;
%     % semilogy(EbN0dB,WBF_BER_sim(:,i),'-^','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);
% 
%     hold on;
%     % semilogy(EbN0dB,WBF_BER_sim(:,i),'-^','color',colour,'DisplayName',[' BER I=' num2str(I_lim(i))]);
%     % semilogy(EbN0dB,combined_data.imwbf.(NK_set{2}).BER(8,:,I_lim(i)/2),'--*','color',colour,'DisplayName',[NK_set{2} ' BER I=' num2str(I_lim(i))]);
%     % plot(EbN0dB,WBF_BER_sim(:,i),'-^','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);
%     hold on;
%     % semilogy(EbN0dB,WBF_FER_sim(:,i),'--*','color',colour,'DisplayName',[' BLER I=' num2str(I_lim(i))]);
%     % semilogy(EbN0dB,combined_data.imwbf.(NK_set{2}).FER(8,:,I_lim(i)/2),'---*','color',colour,'DisplayName',[NK_set{2} ' BLER I=' num2str(I_lim(i))]);
%     % plot(EbN0dB,WBF_FER_sim(:,i),'-*','color',colour,'DisplayName',['BLER (WBF Algo) I=' num2str(I_lim(i))]);
% 
% end
% for i=1:numel(I_lim)
%     colour=hsv2rgb([i/size(I_lim,2)*0.9,1,0.7])
%     % hold on;
%     % semilogy(EbN0dB,WBF_BER_sim(:,i),'-^','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);
%
%     hold on;
%     % semilogy(EbN0dB,combined_data.imwbf.(NK_set{1}).BER(8,:,I_lim(i)/2),'-^','color',colour,'DisplayName',[NK_set{1} ' BER I=' num2str(I_lim(i))]);
%     semilogy(EbN0dB,combined_data.(algo).(NK_set{2}).BER(8,:,I_lim(i)/2),'--*','color',colour,'DisplayName',[NK_set{2} ' BER I=' num2str(I_lim(i))]);
%     % plot(EbN0dB,WBF_BER_sim(:,i),'-^','color',colour,'DisplayName',['BER (WBF Algo) I=' num2str(I_lim(i))]);
%     hold on;
%     % semilogy(EbN0dB,combined_data.imwbf.(NK_set{1}).FER(8,:,I_lim(i)/2),'--*','color',colour,'DisplayName',[NK_set{1} ' BLER I=' num2str(I_lim(i))]);
%     semilogy(EbN0dB,combined_data.(algo).(NK_set{2}).FER(8,:,I_lim(i)/2),'---*','color',colour,'DisplayName',[NK_set{2} ' BLER I=' num2str(I_lim(i))]);
%     % plot(EbN0dB,WBF_FER_sim(:,i),'-*','color',colour,'DisplayName',['BLER (WBF Algo) I=' num2str(I_lim(i))]);
%
% end
legend;
xlabel('Eb/N0 (dB)');
ylabel('Error Rates');
hold off;


