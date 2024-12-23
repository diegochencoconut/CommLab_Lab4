clear; clc; close all;

%% Define parameters
md_order = 4;
symbolnum = 100;
FFT_size = 64;
data_subcarrier_num = 48;
N = symbolnum * data_subcarrier_num;      % For each symbol, it will have 48 subcarriers
CP_size = 16;
freq_sample = 10e6;
expect_sig_length = 0.2;
MAXLOOP = 100;
error_bit = zeros(MAXLOOP, 1);

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

    %% Estimation Using the Short Training
    STS_received = rx_obtained_original(1:length(STS_signal));
    STS_for_alpha_estimation = STS_received(end-32+1:end);
    
    % Estimate alpha_ST
    alpha_ST = 0;
    for ii = 1:16
        alpha_ST = alpha_ST + conj(STS_for_alpha_estimation(ii))*STS_for_alpha_estimation(ii+16);
    end
    
    alpha_ST = phase(alpha_ST);
    alpha_ST = alpha_ST/16;
    
    %% Estimation and Correction Using Long Training
    LTS_received = rx_obtained_original(length(STS_signal)+1+CP_size*2:length(STS_signal)+length(LTS_signal));   
    alpha_LT = 0;
    
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

    %% Power equalization
    rx_power = rms(rx_cfo_corrected)^2;
    tx_power = rms(real_signal)^2;
    rx_obtained = rx_cfo_corrected * sqrt((tx_power) / (rx_power));

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
    
    %% Demodulate the received symbols
    rx_demoded = qamdemod(rx_recovered, md_order);
    rx_bit = int2bit(rx_demoded, log2(md_order));
    
    %% Compare to the transmitted bit and obtain the bit error rate
    dataIn = int2bit(dataSymbolIn, log2(md_order));
    error_bit(loop) = biterr(rx_bit, dataIn);
    fprintf("Bit error rate: %f\n", error_bit(loop) / length(dataIn));
    
    if error_bit(loop) / length(dataIn) <= 0.1
        break;
    end
end

%% Release USRP
release(radio_tx);
release(radio_rx);
    
%% Final plot

ideal_modulation = qammod(0:md_order-1, md_order);

figure;

set(gcf, "Position", [300,300,560,280]);
subplot(1, 2, 1);
rx_recovered_original_0 = rx_recovered_original(dataMod == qammod(0, md_order));
rx_recovered_original_1 = rx_recovered_original(dataMod == qammod(1, md_order));
rx_recovered_original_2 = rx_recovered_original(dataMod == qammod(2, md_order));
rx_recovered_original_3 = rx_recovered_original(dataMod == qammod(3, md_order));
scatter(real(rx_recovered_original_0), imag(rx_recovered_original_0), 'Marker', '.'); hold on
scatter(real(rx_recovered_original_1), imag(rx_recovered_original_1), 'Marker', '.');
scatter(real(rx_recovered_original_2), imag(rx_recovered_original_2), 'Marker', '.'); 
scatter(real(rx_recovered_original_3), imag(rx_recovered_original_3), 'Marker', '.'); 
scatter(real(ideal_modulation), imag(ideal_modulation));
xlim([-2, 2]);
ylim([-2, 2]);
legend('0', '1', '2', '3', 'ideal');

title("Before equalization");
subplot(1, 2, 2);
rx_recovered_0 = rx_recovered(dataMod == qammod(0, md_order));
rx_recovered_1 = rx_recovered(dataMod == qammod(1, md_order));
rx_recovered_2 = rx_recovered(dataMod == qammod(2, md_order));
rx_recovered_3 = rx_recovered(dataMod == qammod(3, md_order));
scatter(real(rx_recovered_0), imag(rx_recovered_0), 'Marker', '.'); hold on
scatter(real(rx_recovered_1), imag(rx_recovered_1), 'Marker', '.');
scatter(real(rx_recovered_2), imag(rx_recovered_2), 'Marker', '.'); 
scatter(real(rx_recovered_3), imag(rx_recovered_3), 'Marker', '.'); 
scatter(real(ideal_modulation), imag(ideal_modulation));
xlim([-2, 2]);
ylim([-2, 2]);

title("After equalization");
sgtitle("Signal Constellation");