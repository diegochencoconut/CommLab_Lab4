clear; clc; close all;

%% Define parameters
md_order = 4;
symbolnum = 100;
FFT_size = 64;
data_subcarrier_num = 48;
N = symbolnum * data_subcarrier_num;      % For each symbol, it will have 48 subcarriers
CP_size = 16;
freq_sample = 10e6;
expect_sig_length = 0.1;
MAXLOOP = 3;

%% USRP parameters
MasterClockRate = 100e6;
freq_carrier = 892e6;      % For VERT900: 824MHz ~ 960MHz
freq_upsampled = MasterClockRate;
interpolate_factor = freq_upsampled / freq_sample;

%% A. Generate STS
STS_position = [-24, -20, -16, -12, -8, -4, 0, 4, 8, 12, 16, 20, 24]' + 33;
STS_value = [1, -1, 1, -1, -1, 1, 0, -1, -1, 1, 1, 1, 1]';
STS_signalf = zeros(64, 1);
STS_signalf(STS_position) = sqrt(13/6) * (1+1i) * STS_value;
STS_signalf = [STS_signalf(33:end); STS_signalf(1:32)];
STS_signal = ifft(STS_signalf);
STS_signal = [STS_signal(33:64); STS_signal; STS_signal];

%% B. Generate LTS
LTS_signalf = [ 1  1 -1 -1 ...
                1  1 -1  1 ...
               -1  1  1  1 ...
                1  1  1 -1 ...
               -1  1  1 -1 ...
                1 -1  1  1 ...
                1  1  0  1 ...
               -1 -1  1  1 ...
               -1  1 -1  1 ...
               -1 -1 -1 -1 ...
               -1  1  1 -1 ...
               -1  1 -1  1 ...
               -1  1  1  1  1]';
LTS_signalf = [zeros(6, 1); LTS_signalf; zeros(5, 1)];
LTS_signalf = [LTS_signalf(33:end); LTS_signalf(1:32)];
LTS_signal = ifft(LTS_signalf);
LTS_signal = [LTS_signal(33:64); LTS_signal; LTS_signal];

%% Generate random bits
dataIn = randi([0, 1], log2(md_order)*N, 1);
dataSymbolIn = bit2int(dataIn, log2(md_order));

%% Modulate the bits with QAM
dataMod = qammod(dataSymbolIn, md_order);
original_modulated_pilot = [1; 1; 1; -1];
txwithcp = zeros((CP_size+FFT_size) * symbolnum, 1);
for ii = 1:symbolnum
    % adding pilot tones
    modulated_pilot = original_modulated_pilot;
    retracted_data = dataMod((ii-1)*data_subcarrier_num+1:ii*data_subcarrier_num);
    subcarrier = [zeros(6,1); retracted_data(1:5);   modulated_pilot(1); ...
                              retracted_data(6:18);  modulated_pilot(2); ...
                              retracted_data(19:24); 0; ...
                              retracted_data(25:30); modulated_pilot(3); ...
                              retracted_data(31:43); modulated_pilot(4); ...
                              retracted_data(44:48); zeros(5,1)];
    
    % adding pilot tones
    subcarrier = fftshift(subcarrier);
    current_symbol = ifft(subcarrier);
    current_symbol_withcp = [current_symbol(end-CP_size+1:end); current_symbol];
    txwithcp((ii-1)*(CP_size+FFT_size)+1:ii*(CP_size+FFT_size)) = current_symbol_withcp;
end
 
preamble_signal = [STS_signal; LTS_signal];
preamble_power = rms(preamble_signal)^2;
signal_power = rms(txwithcp)^2;
preamble_signal = preamble_signal * sqrt((signal_power) / (preamble_power));

real_signal = [preamble_signal; txwithcp];
total_signal = [zeros(30, 1); real_signal; zeros(30, 1)];
transmitted_copies = floor(expect_sig_length * freq_sample / length(total_signal));

transmitted_signal = zeros(transmitted_copies * length(total_signal), 1);
for ii = 1:transmitted_copies
    transmitted_signal((ii-1)*length(total_signal) + 1:ii*length(total_signal)) = total_signal;
end
total_length = expect_sig_length*freq_sample;
num_compensate_zeros = total_length - length(transmitted_signal);

if mod(num_compensate_zeros, 2) == 1
    transmitted_signal = [zeros(floor(num_compensate_zeros/2), 1); transmitted_signal; zeros(floor(num_compensate_zeros/2)+1, 1)];
else
    transmitted_signal = [zeros(num_compensate_zeros/2, 1); transmitted_signal; zeros(num_compensate_zeros/2, 1)];
end

%% Transmit data by USRP
radio_tx = comm.SDRuTransmitter(...
    "Platform",                         "N200/N210/USRP2", ...
    "CenterFrequency",                  freq_carrier, ...
    "Gain",                             4, ...
    "InterpolationFactor",              interpolate_factor);

radio_rx = comm.SDRuReceiver(...
                'Platform',             "N200/N210/USRP2", ...
                'CenterFrequency',      freq_carrier, ...
                'Gain',                 1, ...
                'DecimationFactor',     interpolate_factor, ...
                'SamplesPerFrame',      length(transmitted_signal)*5, ...
                'OutputDataType',       'double');
release(radio_tx);
release(radio_rx);

ber_list = ones(1, MAXLOOP);
for loop = 1:MAXLOOP
    fprintf("Trial %d\n", loop);
    disp('Transmitting......')
    tic
    underrun = radio_tx(transmitted_signal);
    while (underrun ~= 0)
        % pause(0.1);
        underrun = radio_tx(transmitted_signal);
    end
    tx_time = toc;
    fprintf("Transmission time: %f\n", tx_time)
    
    %% Receive data by USRP
    disp('Receiving......')
    tic
    [rx_data, ~, overflow, rx_time_stamp] = radio_rx();   
    while overflow ~= 0
        % pause(0.1);
        [rx_data, ~, overflow, rx_time_stamp] = radio_rx();   
    end
    rx_time = toc;
    fprintf("Receive time: %f\n", rx_time)

    tic
    %% Synchronization
    [result, lags] = xcorr(rx_data, STS_signal);
    [resulttest, lagstest] = xcorr(transmitted_signal, preamble_signal);
    suspected_startpoint = lags(abs(result) == max(abs(result))) + 1;
    suspected_startpoint = suspected_startpoint(1);
    
    %% Obtaining data
    try
        rx_obtained_original = rx_data(suspected_startpoint:suspected_startpoint+length(real_signal)-1);
    catch
        disp("The detection is too late!");
        continue;
    end
    
    %% Estimation and Correction of Carrier Frequency Offset
    % Estimation Using the Short Training
    STS_received = rx_obtained_original(1:length(STS_signal));
    STS_for_alpha_estimation = STS_received(end-32+1:end);
    
    alpha_ST = 0;
    for ii = 1:16
        alpha_ST = alpha_ST + conj(STS_for_alpha_estimation(ii))*STS_for_alpha_estimation(ii+16);
    end
    
    alpha_ST = phase(alpha_ST);
    alpha_ST = alpha_ST/16;
    
    % Estimation and Correction Using Long Training
    LTS_received = rx_obtained_original(length(STS_signal)+1+CP_size*2:length(STS_signal)+length(LTS_signal));
        
    alpha_LT = 0;  
    for ii = 1:64
        alpha_LT = alpha_LT + conj(LTS_received(ii))*LTS_received(ii+64);
    end
    
    alpha_LT = phase(alpha_LT);
    alpha_LT = alpha_LT/64;
    
    deltaf_ST = alpha_ST*freq_sample/(2*pi*16);
    deltaf_LT = alpha_LT*freq_sample/(2*pi*64);
    
    time = (0:length(rx_obtained_original)-1) / freq_sample;
    time = time';
    rx_cfo_corrected = rx_obtained_original .* exp(-1i*2*pi*(deltaf_ST+deltaf_LT).*time);
    rx_obtained_CFO = rx_cfo_corrected;

    %% Channel Estimation
    LTS_cfo_corrected = rx_cfo_corrected(length(STS_signal)+1+CP_size*2:length(STS_signal)+length(LTS_signal));
    LTS_for_channel_estimation = 0.5*(LTS_cfo_corrected(1:length(LTS_signalf)) + LTS_cfo_corrected(length(LTS_signalf)+1:end));
    LTS_cfo_correctedf = fft(LTS_for_channel_estimation);
    channel_estimated = fftshift(LTS_cfo_correctedf) .* fftshift(LTS_signalf);
    channel_estimated_inv = 1 ./ channel_estimated;
    channel_estimated_inv(channel_estimated == 0) = 0;

    %% Process OFDM symbols
    rx_recovered = zeros(data_subcarrier_num*symbolnum, 1);
    rx_remove_preamble = rx_cfo_corrected(length(preamble_signal)+1:end);
    for ii = 1:symbolnum
        rx_currentsymbol = rx_remove_preamble((ii - 1)*(FFT_size + CP_size)+1:ii*(FFT_size + CP_size));
    
        % remove CP
        rx_cpremoved = rx_currentsymbol(CP_size+1:end);
    
        % recover symbols
        rx_recovered_symbol = fft(rx_cpremoved);
        rx_recovered_symbol = fftshift(rx_recovered_symbol);

        % remove channel effect
        rx_channel_fixed = rx_recovered_symbol .* channel_estimated_inv;

        % obtain subcarriers
        rx_dataonly = [rx_channel_fixed( 7:11); ...
                       rx_channel_fixed(13:25); ...
                       rx_channel_fixed(27:32); ...
                       rx_channel_fixed(34:39); ...
                       rx_channel_fixed(41:53); ...
                       rx_channel_fixed(55:59)];
    
        % find offset with pilot tones
        rx_pilotonly = [rx_channel_fixed(12), rx_channel_fixed(26), rx_channel_fixed(40), rx_channel_fixed(54)];
        rx_current_offset = 0.25*sum(original_modulated_pilot .* rx_pilotonly, 'all');
        
        % remove offset
        rx_recovered((ii-1)*data_subcarrier_num+1:ii*data_subcarrier_num) = rx_dataonly ./ rx_current_offset;
    end
    
    %% Demodulate the received symbols
    rx_demoded = qamdemod(rx_recovered, md_order);
    rx_bit = int2bit(rx_demoded, log2(md_order));
    process_time = toc;
    fprintf("Process time: %f\n", process_time);
    fprintf("\n");

    %% Compare to the transmitted bit and obtain the bit error rate
    dataIn = int2bit(dataSymbolIn, log2(md_order));
    error_bit = biterr(rx_bit, dataIn);
    fprintf("Bit error rate: %f\n", error_bit / length(dataIn));
    ber_list(loop) = error_bit / length(dataIn);

    % if error_bit / length(dataIn) <= 0.1
    %     break;
    % end
end

%% Release USRP
release(radio_tx);
release(radio_rx);
    
%% Final plot
% (a)
figure;

set(gcf, "Position", [300,300,560,400]);
subplot(2,1,1);
rxtime_us = (0:length(rx_obtained_original)-1)/freq_sample*1e6;
plot(rxtime_us, abs(rx_obtained_original)); hold on
plot(rxtime_us(1:length(STS_signal)), abs(rx_obtained_original(1:length(STS_signal))));
plot(rxtime_us(length(STS_signal):length(STS_signal)+length(LTS_signal)), abs(rx_obtained_original(length(STS_signal):length(STS_signal)+length(LTS_signal))));
title("Obtained Rx signal (Power unmodified)");
legend("data", "STS", "LTS", 'Location','best');

subplot(2,1,2);
txtime_us = (0:length(real_signal)-1)/freq_sample*1e6;
plot(txtime_us, abs(real_signal)); hold on
plot(txtime_us(1:length(STS_signal)), abs(real_signal(1:length(STS_signal))));
plot(txtime_us(length(STS_signal):length(STS_signal)+length(LTS_signal)), abs(real_signal(length(STS_signal):length(STS_signal)+length(LTS_signal))));
title("Tx signal");

% (b)
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
ideal_modulation = qammod(0:md_order-1, md_order);
color = [0.0000 0.4470 0.7410;
         0.8500 0.3250 0.0980;
         0.6350 0.0780 0.1840;
         0.4940 0.1840 0.5560];

figure;

set(gcf, "Position", [200,150,750,700]);
subplot(2, 2, 1);
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

subplot(2, 2, 2);
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

subplot(2, 2, 3);
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

subplot(2, 2, 4);
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
title("After equalization");
sgtitle("Signal Constellation");