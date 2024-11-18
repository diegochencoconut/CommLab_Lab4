clear; clc; close all;

%% Define parameters
md_order = 4;
symbolnum = 1;
FFT_size = 64;
data_subcarrier_num = 48;
N = symbolnum * data_subcarrier_num;      % For each symbol, it will have 48 subcarriers
CP_size = 16;
freq_sample = 10e6;

transmitted_copies = 40;

%% USRP parameters
MasterClockRate = 100e6;
freq_carrier = 900e6;      % For VERT900: 824MHz ~ 960MHz
freq_upsampled = MasterClockRate;
interpolate_factor = freq_upsampled / freq_sample;

%% make signal
sig_time = 0.5;
data = randi([0, 1], sig_time * freq_sample, 1);

transmitted_signal = pskmod(data, 2);
transmitted_signal = [zeros(100, 1); transmitted_signal; zeros(100, 1)];

rx_signal_list = zeros(length(transmitted_signal), transmitted_copies);

radio_tx = comm.SDRuTransmitter(...
    "Platform",                         "N200/N210/USRP2", ...
    "CenterFrequency",                  freq_carrier, ...
    "Gain",                             15, ...
    "InterpolationFactor",              interpolate_factor);

radio_rx = comm.SDRuReceiver(...
                'Platform',             "N200/N210/USRP2", ...
                'CenterFrequency',      freq_carrier, ...
                'Gain',                 10, ...
                'DecimationFactor',     interpolate_factor, ...
                'SamplesPerFrame',      length(transmitted_signal), ...
                'OutputDataType',       'double');
release(radio_tx);
release(radio_rx);

for ii = 1:transmitted_copies
    disp(ii);

    %% Transmit data by USRP
    disp('Transmitting......')
    underrun = radio_tx(transmitted_signal);
    % if underrun == 0
    %     disp('Transmission successful!')
    % else
    %     error('Transmission failed')
    % end
    fprintf('\n');

    %% Receive data by USRP

    disp('Receiving......')
    [rx_data, ~, overflow, rx_time_stamp] = radio_rx();
    % if overflow==0
    %     disp('Reception successful')
    % else
    %     disp('Reception failed')
    % end

    rx_signal_list(:, ii) = rx_data;
end

%% Save figures
time = (0:length(transmitted_signal)-1)/freq_sample;
time = time';
for ii = 1:transmitted_copies
    disp(ii);
    current_fig = figure('Visible','off');
    subplot(2,1,1);
    plot(time, real(transmitted_signal));
    subplot(2,1,2);
    plot(time, real(rx_signal_list(:, ii)));
    filename = sprintf("img/trial%d", ii);
    saveas(current_fig, filename, "png");
end
% 
% %% Release USRP
% release(radio_tx);
% release(radio_rx);
% % rx_data = transmitted_signal;
% 
% % Recover power
% preamble_power = rms(transmitted_signal)^2;
% signal_power = rms(rx_data)^2;
% % rx_powermodified = rx_data * sqrt((preamble_power) / (signal_power));
% rx_powermodified = rx_data;
% 
% %%
% % figure;
% % subplot(2,1,1)
% % plot(time, abs(transmitted_signal));
% % subplot(2,1,2)
% % plot(time, abs(rx_data));
% % 
% % figure;
% % n = length(rx_data);
% % freq = (-n/2:n/2-1)*(freq_sample/n);
% % subplot(2,1,1)
% % plot(freq, real(fftshift(fft(transmitted_signal))));
% % subplot(2,1,2)
% % plot(freq, real(fftshift(fft(rx_data))));