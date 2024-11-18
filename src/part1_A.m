clear; clc; close all;

%% Define parameters
md_order = 4;
symbolnum = 5;
FFT_size = 64;
data_subcarrier_num = 48;
N = symbolnum * data_subcarrier_num;      % For each symbol, it will have 48 subcarriers
CP_size = 16;
freq_sample = 10e6;
expect_sig_length = 0.2;
MAXLOOP = 100;

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

% figure;
% set(gcf, "Position", [300,300,560,200]);
% time = 0:1/freq_sample:(length(STS_signal)-1)/freq_sample;
% time = time';
% plot(time*1e6, abs(STS_signal));
% title("STS signal");
% xlabel("Time[\mus]");
% ylabel("Magnitude");

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

% figure;
% set(gcf, "Position", [300,300,560,200]);
% time = 0:1/freq_sample:(length(LTS_signal)-1)/freq_sample;
% time = time';
% plot(time*1e6, abs(LTS_signal));
% title("LTS signal");
% xlabel("Time[\mus]");
% ylabel("Magnitude");

%% Generate random bits
dataIn = randi([0, 1], log2(md_order)*N, 1);
dataSymbolIn = bit2int(dataIn, log2(md_order));
% dataSymbolIn = [0	0	0	0	1	1	0	0	0	0	0	3	2	0	2	0	0	2	3	1	3	2	2	1	2	1	3	2	2	2	3	0	0	1	1	0	0	3	3	0	0	1	3	0	3	3	3	1	0	0	2	2	1	1	2	1	3	0	3	3	2	2	1	1	2	0	2	0	0	0	1	1	2	1	0	2	3	3	1	1	3	0	1	3	2	3	1	0	2	3	3	3	0	2	0	3	2	1	3	0	3	2	2	2	0	3	0	1	2	0	3	0	2	1	1	3	1	2	2	0	3	1	1	2	2	2	1	2	1	1	3	0	3	3	2	2	0	3	3	3	2	0	3	3	1	2	2	0	1	0	0	1	0	1	1	2	3	1	2	1	3	1	3	2	0	0	1	2	0	0	1	1	0	0	1	1	3	0	1	2	1	2	3	1	0	3	1	3	0	1	0	3	0	1	0	2	3	0	3	3	1	0	1	3	2	2	3	2	1	3	3	3	1	0	2	1	0	0	1	0	0	2	3	3	0	2	3	2	2	0	3	3	1	1	3	3	0	2	0	1]';

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
                'SamplesPerFrame',      length(transmitted_signal), ...
                'OutputDataType',       'double');
release(radio_tx);
release(radio_rx);

for loop = 1:MAXLOOP
    fprintf("Trial %d\n", loop);
    disp('Transmitting......')
    underrun = radio_tx(transmitted_signal);
    
    %% Receive data by USRP
    disp('Receiving......')
    [rx_data, ~, overflow, rx_time_stamp] = radio_rx();   
    
    %% Synchronization
    [result, lags] = xcorr(rx_data, preamble_signal);
    [resulttest, lagstest] = xcorr(transmitted_signal, preamble_signal);
    suspected_startpoint = lags(abs(result) == max(abs(result))) + 1;
    suspected_startpoint = suspected_startpoint(1);
    % figure;
    % set(gcf, "Position", [300,300,560,400]);
    % subplot(2, 1, 1)
    % plot(lags, abs(result));
    % xlabel("Lags");
    % ylabel("Magnitude");
    % title("Result of match filter");
    % xlim([min(lags), max(lags)])
    % subplot(2, 1, 2)
    % plot(lagstest, abs(resulttest));
    % xlim([min(lags), max(lags)])
    
    %% Obtaining data
    try
        rx_obtained_original = rx_data(suspected_startpoint:suspected_startpoint+length(real_signal)-1);
    catch
        disp("The detection is too late!");
        continue;
    end
    
    preamble_signal = [STS_signal; LTS_signal];
    rx_power = rms(rx_obtained_original)^2;
    tx_power = rms(real_signal)^2;
    rx_obtained = rx_obtained_original * sqrt((tx_power) / (rx_power));

    %% Process OFDM symbols
    rx_recovered = zeros(data_subcarrier_num*symbolnum, 1);
    rx_remove_preamble = rx_obtained(length(preamble_signal)+1:end);
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
        rx_recovered((ii-1)*data_subcarrier_num+1:ii*data_subcarrier_num) = rx_dataonly;
    end
    
    %% Demodulate the received symbols
    rx_demoded = qamdemod(rx_recovered, md_order);
    rx_bit = int2bit(rx_demoded, log2(md_order));
    
    %% Compare to the transmitted bit and obtain the bit error rate
    dataIn = int2bit(dataSymbolIn, log2(md_order));
    error_bit = biterr(rx_bit, dataIn);
    fprintf("Bit error rate: %f\n", error_bit / length(dataIn));

    if error_bit / length(dataIn) <= 0.1
        break;
    end
end

%% Release USRP
release(radio_tx);
release(radio_rx);
    
%% Final plot
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

