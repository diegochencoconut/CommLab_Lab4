clear; clc; close all
load sample.mat

%% Estimation and Correction of Carrier Frequency Offset
% Estimation Using the Short Training
STS_received = rx_obtained_original(1:length(STS_signal));
STS_for_alpha_estimation = STS_received(end-32+1:end);

% Estimate alpha_ST
alpha_ST = 0;
for ii = 1:16
    alpha_ST = alpha_ST + conj(STS_for_alpha_estimation(ii))*STS_for_alpha_estimation(ii+16);
end

alpha_ST = phase(alpha_ST);
alpha_ST = alpha_ST/16;

% Estimation and Correction Using Long Training
LTS_received = rx_obtained_original(length(STS_signal)+1+CP_size*2:length(STS_signal)+length(LTS_signal));

% m = 0:127;
% m = m';
% LTS_alphaST_corrected = LTS_received .* exp(-1i * m * alpha_ST);

alpha_LT = 0;
% for ii = 1:64
%     alpha_LT = alpha_LT + conj(LTS_alphaST_corrected(ii))*LTS_alphaST_corrected(ii+64);
% end
% 
% alpha_LT = phase(alpha_LT);
% alpha_LT = alpha_LT/64;

for ii = 1:64
    alpha_LT = alpha_LT + conj(LTS_received(ii))*LTS_received(ii+64);
end

alpha_LT = phase(alpha_LT);
alpha_LT = alpha_LT/64;

% LTS_alphaLT_corrected = LTS_alphaST_corrected .* exp(-1i * m * alpha_LT);
% LTS_alphaLT_corrected = LTS_received .* exp(-1i * m * alpha_LT);

% alpha = alpha_ST + alpha_LT;
deltaf_ST = alpha_ST*freq_sample/(2*pi*16);
deltaf_LT = alpha_LT*freq_sample/(2*pi*64);

time = (0:length(rx_obtained_original)-1) / freq_sample;
time = time';
rx_cfo_corrected = rx_obtained_original .* exp(-1i*2*pi*(deltaf_ST+deltaf_LT).*time);
% rx_cfo_corrected = rx_obtained_original;

% %% Data correction
% 
% n_correction = 0:(length(rx_obtained_original) - 321);
% n_correction = n_correction' + 128;
% 
% data_received = rx_obtained_original(321:end);
% data_cfo_corrected = data_received .* exp(-1i * n_correction * alpha);
% 
% rx_cfo_corrected = [rx_obtained_original(1:length(preamble_signal)); data_cfo_corrected];

%% Channel Estimation
LTS_cfo_corrected = rx_cfo_corrected(length(STS_signal)+1+CP_size*2:length(STS_signal)+length(LTS_signal));
LTS_for_channel_estimation = 0.5*(LTS_cfo_corrected(1:length(LTS_signalf)) + LTS_cfo_corrected(length(LTS_signalf)+1:end));
LTS_cfo_correctedf = fft(LTS_for_channel_estimation);
channel_estimated = fftshift(LTS_cfo_correctedf) .* fftshift(LTS_signalf);
channel_estimated_inv = 1 ./ channel_estimated;
channel_estimated_inv(channel_estimated == 0) = 0;

%% Modify power
% rx_power = rms(rx_cfo_corrected)^2;
% tx_power = rms(preamble_signal)^2;
% rx_obtained_CFO = rx_cfo_corrected * sqrt((tx_power) / (rx_power));
rx_obtained_CFO = rx_cfo_corrected;

%% Process OFDM symbols
rx_recovered_original = zeros(data_subcarrier_num*symbolnum, 1);
rx_remove_preamble_original = rx_obtained_original(length(preamble_signal)+1:end);
for ii = 1:symbolnum
    rx_currentsymbol_original = rx_remove_preamble_original((ii - 1)*(FFT_size + CP_size)+1:ii*(FFT_size + CP_size));

    % remove CP
    rx_cpremoved_original = rx_currentsymbol_original(CP_size+1:end);

    % recover and remove pilot tones
    rx_recovered_symbol_original = fft(rx_cpremoved_original);
    rx_recovered_symbol_original = fftshift(rx_recovered_symbol_original);
    rx_dataonly_original = [rx_recovered_symbol_original( 7:11); ...
                   rx_recovered_symbol_original(13:25); ...
                   rx_recovered_symbol_original(27:32); ...
                   rx_recovered_symbol_original(34:39); ...
                   rx_recovered_symbol_original(41:53); ...
                   rx_recovered_symbol_original(55:59)];
    rx_recovered_original((ii-1)*data_subcarrier_num+1:ii*data_subcarrier_num) = rx_dataonly_original;
end

rx_recovered_CFO = zeros(data_subcarrier_num*symbolnum, 1);
rx_remove_preamble = rx_obtained_CFO(length(preamble_signal)+1:end);
for ii = 1:symbolnum
    rx_currentsymbol = rx_remove_preamble((ii - 1)*(FFT_size + CP_size)+1:ii*(FFT_size + CP_size));

    % remove CP
    rx_cpremoved = rx_currentsymbol(CP_size+1:end);

    % recover and remove pilot tones
    rx_recovered_symbol = fft(rx_cpremoved);
    rx_recovered_symbol = fftshift(rx_recovered_symbol);
    rx_dataonly = [rx_recovered_symbol( 7:11); ...
                   rx_recovered_symbol(13:25); ...
                   rx_recovered_symbol(27:32); ...
                   rx_recovered_symbol(34:39); ...
                   rx_recovered_symbol(41:53); ...
                   rx_recovered_symbol(55:59)];
    rx_recovered_CFO((ii-1)*data_subcarrier_num+1:ii*data_subcarrier_num) = rx_dataonly;
end

rx_recovered_wchannel_estimation = zeros(data_subcarrier_num*symbolnum, 1);
rx_remove_preamble = rx_obtained_CFO(length(preamble_signal)+1:end);
for ii = 1:symbolnum
    rx_currentsymbol = rx_remove_preamble((ii - 1)*(FFT_size + CP_size)+1:ii*(FFT_size + CP_size));

    % remove CP
    rx_cpremoved = rx_currentsymbol(CP_size+1:end);

    % recover and remove pilot tones
    rx_recovered_symbol = fft(rx_cpremoved);
    rx_recovered_symbol = fftshift(rx_recovered_symbol);
    rx_channel_fixed = rx_recovered_symbol .* channel_estimated_inv;
    rx_dataonly = [rx_channel_fixed( 7:11); ...
                   rx_channel_fixed(13:25); ...
                   rx_channel_fixed(27:32); ...
                   rx_channel_fixed(34:39); ...
                   rx_channel_fixed(41:53); ...
                   rx_channel_fixed(55:59)];
    rx_recovered_wchannel_estimation((ii-1)*data_subcarrier_num+1:ii*data_subcarrier_num) = rx_dataonly;
end

rx_recovered_wequalization = zeros(data_subcarrier_num*symbolnum, 1);
rx_remove_preamble = rx_obtained_CFO(length(preamble_signal)+1:end);
for ii = 1:symbolnum
    rx_currentsymbol = rx_remove_preamble((ii - 1)*(FFT_size + CP_size)+1:ii*(FFT_size + CP_size));

    % remove CP
    rx_cpremoved = rx_currentsymbol(CP_size+1:end);

    % recover and remove pilot tones
    rx_recovered_symbol = fft(rx_cpremoved);
    rx_recovered_symbol = fftshift(rx_recovered_symbol);
    rx_channel_fixed = rx_recovered_symbol .* channel_estimated_inv;
    rx_dataonly = [rx_channel_fixed( 7:11); ...
                   rx_channel_fixed(13:25); ...
                   rx_channel_fixed(27:32); ...
                   rx_channel_fixed(34:39); ...
                   rx_channel_fixed(41:53); ...
                   rx_channel_fixed(55:59)];

    % Find offset with pilot tones
    rx_pilotonly = [rx_channel_fixed(12), rx_channel_fixed(26), rx_channel_fixed(40), rx_channel_fixed(54)];
    rx_current_offset = 0.25*sum(original_modulated_pilot .* rx_pilotonly, 'all');
    
    % Remove offset
    rx_recovered_wequalization((ii-1)*data_subcarrier_num+1:ii*data_subcarrier_num) = rx_dataonly ./ rx_current_offset;
end

%% Demodulate the received symbols
% Original
rx_demoded_original = qamdemod(rx_recovered_original, md_order);
rx_bit_original = int2bit(rx_demoded_original, log2(md_order));

% After CFO
rx_demoded_CFO = qamdemod(rx_recovered_CFO, md_order);
rx_bit_CFO = int2bit(rx_demoded_CFO, log2(md_order));

% After Channel Estimation
rx_demoded_wchannel_estimation = qamdemod(rx_recovered_wchannel_estimation, md_order);
rx_bit_wchannel_estimation = int2bit(rx_demoded_wchannel_estimation, log2(md_order));

% After Equalization
rx_demoded_wequalization = qamdemod(rx_recovered_wequalization, md_order);
rx_bit_wequalization = int2bit(rx_demoded_wequalization, log2(md_order));

%% Compare to the transmitted bit and obtain the bit error rate
% dataIn = int2bit(dataSymbolIn, log2(md_order));

error_bit_original(loop) = biterr(rx_bit_original, dataIn);
fprintf("Original BER: %f\n", error_bit_original(loop) / length(dataIn));

error_bit_CFO(loop) = biterr(rx_bit_CFO, dataIn);
fprintf("BER after CFO: %f\n", error_bit_CFO(loop) / length(dataIn));

error_bit_wchannel_estimation(loop) = biterr(rx_bit_wchannel_estimation, dataIn);
fprintf("BER after channel estimation: %f\n", error_bit_wchannel_estimation(loop) / length(dataIn));

error_bit_wequalization(loop) = biterr(rx_bit_wequalization, dataIn);
fprintf("BER after equalization: %f\n", error_bit_wequalization(loop) / length(dataIn));

%% Final plot

ideal_modulation = qammod(0:md_order-1, md_order);
color = [0.0000 0.4470 0.7410;
         0.8500 0.3250 0.0980;
         0.6350 0.0780 0.1840;
         0.4940 0.1840 0.5560];

figure;

set(gcf, "Position", [300,150,1120,280]);
subplot(1, 4, 1);
rx_recovered_original_0 = rx_recovered_original(dataMod == qammod(0, md_order));
rx_recovered_original_1 = rx_recovered_original(dataMod == qammod(1, md_order));
rx_recovered_original_2 = rx_recovered_original(dataMod == qammod(2, md_order));
rx_recovered_original_3 = rx_recovered_original(dataMod == qammod(3, md_order));
scatter(real(rx_recovered_original_0), imag(rx_recovered_original_0), MarkerEdgeColor=color(1, :), Marker='.'); hold on
scatter(real(rx_recovered_original_1), imag(rx_recovered_original_1), MarkerEdgeColor=color(2, :), Marker='.');
scatter(real(rx_recovered_original_2), imag(rx_recovered_original_2), MarkerEdgeColor=color(3, :), Marker='.'); 
scatter(real(rx_recovered_original_3), imag(rx_recovered_original_3), MarkerEdgeColor=color(4, :), Marker='.'); 
scatter(real(ideal_modulation(1)), imag(ideal_modulation(1)), MarkerEdgeColor=color(1, :));
scatter(real(ideal_modulation(2)), imag(ideal_modulation(2)), MarkerEdgeColor=color(2, :));
scatter(real(ideal_modulation(3)), imag(ideal_modulation(3)), MarkerEdgeColor=color(3, :));
scatter(real(ideal_modulation(4)), imag(ideal_modulation(4)), MarkerEdgeColor=color(4, :));
xline(0, 'black');
yline(0, 'black');
xlim([-2, 2]);
ylim([-2, 2]);
xlabel("Real");
ylabel("Imag");
title("Original Signal");

subplot(1, 4, 2);
rx_recovered_0 = rx_recovered_CFO(dataMod == qammod(0, md_order));
rx_recovered_1 = rx_recovered_CFO(dataMod == qammod(1, md_order));
rx_recovered_2 = rx_recovered_CFO(dataMod == qammod(2, md_order));
rx_recovered_3 = rx_recovered_CFO(dataMod == qammod(3, md_order));
scatter(real(rx_recovered_0), imag(rx_recovered_0), MarkerEdgeColor=color(1, :), Marker='.'); hold on
scatter(real(rx_recovered_1), imag(rx_recovered_1), MarkerEdgeColor=color(2, :), Marker='.');
scatter(real(rx_recovered_2), imag(rx_recovered_2), MarkerEdgeColor=color(3, :), Marker='.'); 
scatter(real(rx_recovered_3), imag(rx_recovered_3), MarkerEdgeColor=color(4, :), Marker='.'); 
scatter(real(ideal_modulation(1)), imag(ideal_modulation(1)), MarkerEdgeColor=color(1, :));
scatter(real(ideal_modulation(2)), imag(ideal_modulation(2)), MarkerEdgeColor=color(2, :));
scatter(real(ideal_modulation(3)), imag(ideal_modulation(3)), MarkerEdgeColor=color(3, :));
scatter(real(ideal_modulation(4)), imag(ideal_modulation(4)), MarkerEdgeColor=color(4, :));
xline(0, 'black');
yline(0, 'black');
xlim([-2, 2]);
ylim([-2, 2]);
xlabel("Real");
ylabel("Imag");
title("After CFO");
sgtitle("Signal Constellation");

subplot(1, 4, 3);
rx_recovered_0 = rx_recovered_wchannel_estimation(dataMod == qammod(0, md_order));
rx_recovered_1 = rx_recovered_wchannel_estimation(dataMod == qammod(1, md_order));
rx_recovered_2 = rx_recovered_wchannel_estimation(dataMod == qammod(2, md_order));
rx_recovered_3 = rx_recovered_wchannel_estimation(dataMod == qammod(3, md_order));
scatter(real(rx_recovered_0), imag(rx_recovered_0), MarkerEdgeColor=color(1, :), Marker='.'); hold on
scatter(real(rx_recovered_1), imag(rx_recovered_1), MarkerEdgeColor=color(2, :), Marker='.');
scatter(real(rx_recovered_2), imag(rx_recovered_2), MarkerEdgeColor=color(3, :), Marker='.'); 
scatter(real(rx_recovered_3), imag(rx_recovered_3), MarkerEdgeColor=color(4, :), Marker='.'); 
scatter(real(ideal_modulation(1)), imag(ideal_modulation(1)), MarkerEdgeColor=color(1, :));
scatter(real(ideal_modulation(2)), imag(ideal_modulation(2)), MarkerEdgeColor=color(2, :));
scatter(real(ideal_modulation(3)), imag(ideal_modulation(3)), MarkerEdgeColor=color(3, :));
scatter(real(ideal_modulation(4)), imag(ideal_modulation(4)), MarkerEdgeColor=color(4, :));
xline(0, 'black');
yline(0, 'black');
xlim([-2, 2]);
ylim([-2, 2]);
xlabel("Real");
ylabel("Imag");
title("After channel estimation");
sgtitle("Signal Constellation");

subplot(1, 4, 4);
rx_recovered_0 = rx_recovered_wequalization(dataMod == qammod(0, md_order));
rx_recovered_1 = rx_recovered_wequalization(dataMod == qammod(1, md_order));
rx_recovered_2 = rx_recovered_wequalization(dataMod == qammod(2, md_order));
rx_recovered_3 = rx_recovered_wequalization(dataMod == qammod(3, md_order));
scatter(real(rx_recovered_0), imag(rx_recovered_0), MarkerEdgeColor=color(1, :), Marker='.'); hold on
scatter(real(rx_recovered_1), imag(rx_recovered_1), MarkerEdgeColor=color(2, :), Marker='.');
scatter(real(rx_recovered_2), imag(rx_recovered_2), MarkerEdgeColor=color(3, :), Marker='.'); 
scatter(real(rx_recovered_3), imag(rx_recovered_3), MarkerEdgeColor=color(4, :), Marker='.'); 
scatter(real(ideal_modulation(1)), imag(ideal_modulation(1)), MarkerEdgeColor=color(1, :));
scatter(real(ideal_modulation(2)), imag(ideal_modulation(2)), MarkerEdgeColor=color(2, :));
scatter(real(ideal_modulation(3)), imag(ideal_modulation(3)), MarkerEdgeColor=color(3, :));
scatter(real(ideal_modulation(4)), imag(ideal_modulation(4)), MarkerEdgeColor=color(4, :));
xline(0, 'black');
yline(0, 'black');
xlim([-2, 2]);
ylim([-2, 2]);
xlabel("Real");
ylabel("Imag");
title("After channel estimation");
sgtitle("Signal Constellation");