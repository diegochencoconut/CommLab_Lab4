clear; clc; close all;

%% Define parameters
md_order = 16;
symbolnum = 100;
FFT_size = 64;
data_subcarrier_num = 48;
N = symbolnum * data_subcarrier_num;      % For each symbol, it will have 48 subcarriers
CP_size = 16;
freq_sample = 10e6;
expect_sig_length = 0.1;
MAXLOOP = 100;
tx_gain = [0.5, 2, 5, 10, 15, 20, 25, 30];
rx_gain = 1;

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
    
%% Transmit data by USRP with varying power
avg_ber = zeros(1,length(tx_gain));
SNR = zeros(1,length(tx_gain));
noise = zeros(length(real_signal), 1);

for power_times = 1:length(tx_gain)
    radio_tx = comm.SDRuTransmitter(...
        "Platform",                         "N200/N210/USRP2", ...
        "CenterFrequency",                  freq_carrier, ...
        "Gain",                             tx_gain(power_times), ...
        "InterpolationFactor",              interpolate_factor);
    
    radio_rx = comm.SDRuReceiver(...
                    'Platform',             "N200/N210/USRP2", ...
                    'CenterFrequency',      freq_carrier, ...
                    'Gain',                 rx_gain, ...
                    'DecimationFactor',     interpolate_factor, ...
                    'SamplesPerFrame',      length(transmitted_signal)*5, ...
                    'OutputDataType',       'double');
    release(radio_tx);
    release(radio_rx);
    
    ber_list = ones(1, MAXLOOP);

    %% Get noise power
    signal_empty = zeros(length(transmitted_signal), 1);
    fprintf("Noise measuring\n");
    disp('Transmitting......')
    tic
    underrun = radio_tx(signal_empty);
    while (underrun ~= 0)
        % pause(0.1);
        underrun = radio_tx(signal_empty);
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

    noise_power = rms(rx_data)^2;
    release(radio_tx);
    release(radio_rx);
    %% Transmit data
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
            loop = loop - 1;
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
    avg_ber(power_times) = sum(ber_list(2:end)) / (MAXLOOP - 1);
    rx_signal_power = rms(rx_obtained_original)^2;
    SNR(power_times) = 10 * log10((rx_signal_power-noise_power) / noise_power);
end

%% Release USRP
release(radio_tx);
release(radio_rx);
    
%% Final plot
% (d) 
figure;
set(gcf, 'Position', [300, 300, 600, 280]);
scatter(SNR(1:length(tx_gain)), avg_ber(1:length(tx_gain)), 'o', 'filled');
set(gca,'yscale','log')
xlabel("SNR");
ylabel("Average BER");
title("BER v.s. SNR");
xlim([min(SNR(1:length(tx_gain)))-5, max(SNR(1:length(tx_gain)))+5]);

