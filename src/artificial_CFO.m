clear; clc; close all;

%% Define parameters
md_order = 4;
symbolnum = 5;
FFT_size = 64;
data_subcarrier_num = 48;
N = symbolnum * data_subcarrier_num;      % For each symbol, it will have 48 subcarriers
CP_size = 16;
freq_sample = 10e6;
channel_SNR = 0;
freq_offset_num = 0;                          % deltaf = 2*pi

transmitted_copies = 1;

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

figure;
set(gcf, "Position", [300,300,560,200]);
time = 0:1/freq_sample:(length(STS_signal)-1)/freq_sample;
time = time';
plot(time*1e6, abs(STS_signal));
title("STS signal");
xlabel("Time[\mus]");
ylabel("Magnitude");

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

figure;
set(gcf, "Position", [300,300,560,200]);
time = 0:1/freq_sample:(length(LTS_signal)-1)/freq_sample;
time = time';
plot(time*1e6, abs(LTS_signal));
title("LTS signal");
xlabel("Time[\mus]");
ylabel("Magnitude");

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

total_signal = [preamble_signal; txwithcp];

transmitted_signal = zeros(transmitted_copies * length(total_signal), 1);
for i = 1:transmitted_copies
    transmitted_signal((i-1)*length(total_signal) + 1:i*length(total_signal)) = total_signal;
end
transmitted_signal = [zeros(30, 1); transmitted_signal];

%% Upsampling from 20MHz sampling rate to 5GHz sampling rate
tx_upsampled = upsample(transmitted_signal, freq_upsampled / freq_sample);

% interpolation filter
L = freq_upsampled / freq_sample;
tx_upsampledf = fft(tx_upsampled);
omega = 0:2*pi/length(tx_upsampledf):2*pi-2*pi/length(tx_upsampledf);
omega = omega';
interpolated_signalf = tx_upsampledf .*((omega <= pi/L) + (omega >= 2*pi - pi/L));
interpolated_signal = ifft(interpolated_signalf);

%% Perform digital upconversion
time = 0:1/freq_upsampled:1/freq_upsampled*(length(interpolated_signal)-1);
time = time';
interpolated_signal_real = real(interpolated_signal);
interpolated_signal_imag = imag(interpolated_signal);
tx_upconversion = real(interpolated_signal .* exp(1i*2*pi*freq_carrier*time));

%% Add Gaussian noise for a target SNR of 20dB
rx_noised = awgn(tx_upconversion, channel_SNR - 10*log10(250), 'measured');

%% Perform digital downconversion
% rx_downconversion = rx_noised .* exp(-1i*2*pi*freq_carrier*time);
rx_downconversion_real = rx_noised .* 2 .* cos(2*pi*freq_carrier*time);
rx_downconversion_imag = rx_noised .* (-2) .* sin(2*pi*freq_carrier*time);
rx_downconversion = rx_downconversion_real + 1i * rx_downconversion_imag;

rx_downconversionf = fft(rx_downconversion);
decimated_signalf = rx_downconversionf .* ((omega <= pi/L) + (omega >= 2*pi - pi/L));
decimated_signal = ifft(decimated_signalf) * L;

%% Downsampling
rx_received = downsample(decimated_signal, freq_upsampled / freq_sample);

%% Artificial frequency offset
rx_receivedf = fft(rx_received);
rx_receivedf = circshift(rx_receivedf, freq_offset_num);
rx_receivedwCFO = ifft(rx_receivedf);

%% Synchronization
[result, lags] = xcorr(rx_receivedwCFO, preamble_signal);
[resultorigin, lagsorigin] = xcorr(rx_received, preamble_signal);
[resulttest, lagstest] = xcorr(transmitted_signal, preamble_signal);
suspected_startpoint = lags(abs(result) == max(abs(result))) + 1;
figure;
set(gcf, "Position", [300,300,560,400]);
subplot(2, 1, 1)
plot(lagsorigin, abs(resultorigin)); hold on
plot(lags, abs(result)); 
xlabel("Lags");
ylabel("Magnitude");
title("Result of match filter");
xlim([min(lags), max(lags)])
subplot(2, 1, 2)
plot(lagstest, abs(resulttest));
xlim([min(lags), max(lags)])

%% Obtaining data
rx_obtained = rx_receivedwCFO(suspected_startpoint:suspected_startpoint+length(total_signal)-1);

%% Process OFDM symbols
rx_recovered = [];
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
    rx_recovered = [rx_recovered; rx_dataonly];
end

%% Demodulate the received symbols
rx_demoded = qamdemod(rx_recovered, md_order);
rx_bit = int2bit(rx_demoded, log2(md_order));

%% Compare to the transmitted bit and obtain the bit error rate
dataIn = int2bit(dataSymbolIn, log2(md_order));
error_bit = biterr(rx_bit, dataIn);
fprintf("Bit error rate: %f\n", error_bit / length(dataIn));

%% 
fprintf("Real frequency offset: %f\n", freq_offset_num/length(rx_receivedf))