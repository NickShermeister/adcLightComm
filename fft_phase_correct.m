% load('cleanData.mat');

load('paige_plotteddata.mat')
 
bits = y;
 
% Step 1: Estimate channel magnitude
magnitude_estimate = rms(abs(bits));
 
% Step 2: Divide received signal by estimate of channel magnitude
y_hat = bits./magnitude_estimate;
 
% Step 3: Initialize variables
x_hat = zeros(length(y_hat), 1);

fft_x = fft(y_hat.^4);
x_axis = linspace(0, 2*pi*(length(bits)-1)/length(bits), length(bits));

[max_val, max_index] = max(abs(fft_x));

freq_offset = -1*x_axis(max_index)./4;

theta_hat = -1*angle(fft_x(max_index))./4;

for k = 1:length(y_hat)
    x_hat(k) = y_hat(k).*exp(1i*(freq_offset * (k-1) + theta_hat + 3*pi/4));
end

hold on

down_x = downsample(x_hat(10:end), 20);

plot(real(down_x), imag(down_x), 'o');

