% tx_samples_from_file --freq 250e3 --rate 200e3 --type float --gain 70 --file tx2.dat
%Lower receiver gain to like 20.


% Create header noise
seed = 562019;
rng(seed);
constant_noise = sign(randn(1000,1)) + 1i*sign(randn(1000,1));

% Create randomly generated message signal, comment out if using txt file
% to write
N = 10000;
signal = sign(randn(N,1)) + 1i*sign(randn(N,1));
bits = [constant_noise; signal];

Symbol_period = 20;

% create a generic pulse of unit height
% with width equal to symbol period
pulse = ones(Symbol_period, 1);

% spread out the values in "bits" by Symbol_period
% first create a vector to store the values
x = zeros(Symbol_period*length(bits),1);

% assign every Symbol_period-th sample to equal a value from bits
x(1:Symbol_period:end) = bits;

% now convolve the single generic pulse with the spread-out bits
x_tx = conv(pulse, x);
x_tx = x_tx(1:end);

% to visualize, make a stem plot
figure(1)
clf(1)
subplot(2,1,1)
hold on
stem(real(x_tx));
title('real')
hold off
subplot(2,1,2)
hold on
stem(imag(x_tx))
title('imag')
hold off

% zero pad the beginning with 100000 samples to ensure that any glitch that
% happens when we start transmitting doesn't effect the data

x_tx = [zeros(100000, 1); x_tx(1:end-19);zeros(100000, 1)];


% here we write the data into a format that the USRP can understand
% specifically, we use float32 numbers with real followed by imaginary
% values

% first create a vector to store the interleaved real and imaginary values

tmp = zeros(length(x_tx)*2, 1);

% then assign the real part of x_tx to every other sample and the imaginary
% part to the remaining samples. In this example, the imaginary parts are
% all zero since our original signal is purely real, but we still have to
% write the zero values 

tmp(1:2:end) = real(x_tx);
tmp(2:2:end) = imag(x_tx);

% open a file to write in binary format 
f1 = fopen('txtest2.dat', 'wb');
% write the values as a float32
fwrite(f1, tmp/100, 'float32');
% close the file
fclose(f1)