clear; clf;

%Set the symbol period (must be identical to transmit file)
Symbol_period = 20;

% Build identical random noise to transmit file
header_size = 1000;
seed = 562019;
rng(seed);
constant_noise = sign(randn(header_size,1)) + 1i*sign(randn(header_size,1));

% Build identical signal for checking bit error against later
% Comment out if not using randomly generated signal
signal_size = 10000;
signal = sign(randn(signal_size,1)) + 1i*sign(randn(signal_size,1));

% Convolve constant noise with pulse, effectively upsampling constant_noise
pulse = ones(Symbol_period, 1);
header = zeros(Symbol_period*length(constant_noise),1);
header(1:Symbol_period:end) = constant_noise; % every 20 points there is a data value
header_tx = conv(pulse, header);              % makes that data value everywhere at mag of pulse



% Open the file containing the received samples
f2 = fopen('rx_MNJ2.dat', 'rb');
% read data from the file
rxfile = fread(f2, 'float32');
% close the file
fclose(f2);

%Assign y to be the data, ignoring the first 250 values due to glitching
%from hardware
y = zeros(length(rxfile)/2,1);
y = rxfile(1:2:end) + 1i*rxfile(2:2:end);
y = y(250:end);

% Open the file containing the transmitted samples (what you actually sent)

f1 = fopen('txtest.dat', 'rb');
txfile = fread(f1, 'float32');
fclose(f1);
txfile = txfile*100;        
xreal = txfile(1:2:end);    % pulling out the real values from every other index
ximag = txfile(2:2:end);    % pulling out the imag values from the other set of every other index
x = xreal + 1i*ximag;       % add them together to create the 2 channels
x = x(100001:(length(x)-100000));   % remove the 100000 padding zeroes from both sides


% Cross correlate the received signal with our reconstructed header to find
% where the header starts
[Ryx,lags] = xcorr(y,header_tx);

figure(1)
clf(1)
hold on
plot(lags,abs(Ryx))
hold off


% Finds the index and magnitude of the spike of the cross-correlation, this
% marks the start of the header section
[mm, ii] = max(abs(Ryx));

% The start of the actual data should be the index of the spike plus the
% length of the header
start = lags(ii)+1 + header_size*Symbol_period;
stop = start + signal_size*Symbol_period; %% start point plus length of data (10,000) * symbol period
y_data = y(start: stop); % should be only our data


% plot the real and imaginary parts of what we received
figure(2)
clf(2)
subplot(2,2,1);
hold on
plot(real(y_data));
title('Real portion of received Data');
ylabel('Magnitude of real (V)')
hold off

subplot(2,2,2);
hold on
plot(imag(y_data));
title('Imaginary portion of received Data');
xlabel('Index')
ylabel('Magnitude of imag (V)')
hold off

subplot(2,2,3);
hold on
plot(real(signal));
title('Real portion of transmitted Data');
xlabel('Index')
ylabel('Magnitude of imag (V)')
ylim([-1.5 1.5])
hold off

subplot(2,2,4);
hold on
plot(imag(signal));
title('Imaginary portion of transmitted Data');
xlabel('Index')
ylabel('Magnitude of imag (V)')
ylim([-1.5 1.5])
hold off


%%FREQUENCY OFFSET CORRECTION
% Take the fft of the received data
fft_x = fft(y_data.^4); % to the 4th because QPSK
x_axis = linspace(0, 2*pi*(length(y_data)-1)/length(y_data), length(y_data)); % remapping the x-axis to [0 2pi]

% Plot fft, notice the spike near frequency 0 --> this is offset
figure(3)
clf(3)
hold on
plot(x_axis,abs(fft_x))
xlabel('Frequency (rad/sample)')
ylabel('Magnitude')
hold off


% Pull out the value and index of peak
[max_val, max_index] = max(abs(fft_x));

% The frequency at which the spike occurs, is 4 time the frequency offset
freq_offset = x_axis(max_index)./4;

%The angle of the magnitude of the spike is 4 times the phase offset
theta_hat = angle(fft_x(max_index))./4 + 5*pi/4; % + 5*pi/4 dials the constellation into 
                                                 % the correct quadrants 

                                                 
% Creating a time vector to use in correction calculation
times = linspace(0,length(y_data)-1,length(y_data))'; 

% Using complex exponentials, and the offset values we found to calculate
% corrected data
x_hat = y_data.*(exp(-1i*(freq_offset.*times + theta_hat)));

% Downsample corrected data and separate into real and imaginary components
down_x = downsample(x_hat(10:end), 20);
down_x_real = sign(real(down_x));
down_x_imag = sign(imag(down_x));

%Interleave the values to make one long data vector
message_vec = ones(2*length(down_x),1);
message_vec(1:2:end) = down_x_real';
message_vec(2:2:end) = down_x_imag';
message_vec(find(message_vec == -1)) = 0;  % converting -1s to 0 to make binary

bin2str(message_vec)




% Constellation plot of corrected data
figure(4)
clf(4)
hold on
plot(real(down_x), imag(down_x), 'o');
xlabel('Real')
ylabel('Imaginary')
title('Frequency Offset Corrected Data Constellation')
hold off

% To check for accuracy, take our actual transmitted data and separate in
% the same way 
signal_real = real(signal);
signal_imag = imag(signal);

% Check for incongruities, add them, divide by total to calculate bit error
wrong_real = (signal_real ~= down_x_real);
wrong_imag = (signal_imag ~= down_x_imag); 
total_wrong = sum(wrong_real) + sum(wrong_imag);
bit_error = total_wrong/(length(signal)*2)



figure(5)
clf(5)

plot(signal_real(200:300),'.-')
hold on
plot(100* real(down_x(200:300)),'.-')
xlabel('Index')
ylabel('Magnitude(V)')
title('Comparison of tx vs rx in sub section')
ylim([-1.25 1.5])
legend('Transmitted data','Received data *100')



