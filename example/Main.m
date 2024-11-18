clear; clc; close all;
md_order = 4 ;
% Set the rolloff factor used in the RRC filter 
rolloff = 0 ;  
[prmQPSKTransmitter,prmQPSKReceiver,qpskRx,hTx,radio_Tx,radio_Rx] = Tx_Rx_init(rolloff,md_order);

cleanupTx = onCleanup(@()release(hTx));
cleanupRadio = onCleanup(@()release(radio_Tx));
cleanupRadio = onCleanup(@()release(radio_Rx));
cleanupRx = onCleanup(@()release(qpskRx));

currentTime = 0;
underrun = uint32(0);
count = 0 ; 

data = 0 ; 
coarseCompSignal = 0 ;
flag = 0 ; 
buffer = zeros(3,1000); 

buffer_error = 0 ; 
buffer_symbol = 0; 
%while currentTime < prmQPSKTransmitter.StopTime
% Bit generation, modulation and transmission filtering 
while flag == 0     
    buffer_Tx = data ;
    data = hTx();
    buffer_Rx = coarseCompSignal ;
    % Data transmission
    tunderrun = radio_Tx(data);
    underrun = underrun + tunderrun;
    %% Receive data 
    [rcvdSignal, ~, toverflow] = step(radio_Rx);
    if ~toverflow % Avoid overflow frames for continuous-mode receiver synchronization
        %% Decode data 
        [AGCSignal,RCRxSignal, coarseCompSignal, timingRecSignal, fineCompSignal, symframe,BER] = qpskRx(rcvdSignal); % Receiver
    end
    %overflow = toverflow + overflow;
    % Update simulation time
    %% BER here is the acuumulated bit error rate 
    %% We want to find the BER of the current frame 

    symbol = BER(3) - buffer_symbol ; 
    if symbol > 0 
        if ((  BER(2) - buffer_error  )/symbol   < 0.1) && (BER(2) >0)
            flag = 1 ; 
        end 
    end 
    currentTime=currentTime+prmQPSKTransmitter.USRPFrameTime;
    buffer_error = BER(2) ; 
    buffer_symbol = BER(3) ;
    count = count +1 ;
    buffer(:,count) = BER ;
    disp(count);
    disp(BER);

end

cleanupTx = onCleanup(@()release(hTx));
cleanupRadio = onCleanup(@()release(radio_Tx));
cleanupRadio = onCleanup(@()release(radio_Rx));
cleanupRx = onCleanup(@()release(qpskRx));


sampleRate = 1000000;
L = length(data);
Fs = sampleRate ; 

figure ;
plot(Fs/L*(-L/2:L/2-1),abs(fftshift(fft(real(data)))))
title("Spectrum of transmitted signal")
xlabel("frequency")


figure ;
scatter(real(RCRxSignal),imag(RCRxSignal))
title("QPSK received signal after RRC filter")
xlabel("Real part ")
ylabel("imaginary part ")


figure ; 
scatter(real(symframe),imag(symframe))
title("QPSK received signal after processing")
xlabel("Real part ")
ylabel("imaginary part ")