% Phase lock loop: Shove the "fake" stuff back into the "real" stuff

clear;
clf;
% Open the file containing the received samples
f2 = fopen('rxq.dat', 'rb');
% read data from the file
rxfile = fread(f2, 'float32');
% close the file
fclose(f2);

symbol_period = 20;

f1 = fopen('txq.dat', 'rb');
txfile = fread(f1, 'float32');
fclose(f1);
txfile = txfile*100;
xreal = txfile(1:2:end);
ximag = txfile(2:2:end);


y = zeros(length(rxfile)/2,1);
y = rxfile(1:2:end)+j*rxfile(2:2:end);
y = y(250:end);

magnitude_estimate = rms(abs(y));
y = y./magnitude_estimate;

%basic filtering
tempreal = y(1:2:200);
tempimag = y(2:2:200);

maxreal = max(abs(tempreal));
maximag = max(tempimag);

start_constant = 0; % 200000;
end_constant = 0; %200000;

start = -1;
runningsum = zeros(1, 500);
for z = 1:2:length(y)
    runningsum(mod((z-1)/2, 500) + 1) = abs(y(z));
    if (sum(runningsum)/length(runningsum)) > maxreal*3
       start = z - length(runningsum)
       break
    end
end

ends = length(y) - end_constant;
runningsum = zeros(1, 100);
for z = length(y):-2:150
%     abs(y(z))
    runningsum(mod(round((z-1)/2), 100) + 1) = abs(y(z));
    if (sum(runningsum)/length(runningsum)) > maxreal*3
        ends = z + length(runningsum)
        break
    end
end


y = y((start - start_constant):(ends+ end_constant));


xreal = xreal(100000:(length(xreal)-100000));
ximag = ximag(100000:(length(ximag)-100000));

% size_corr = 1000;
% C = xcorr((y(1:size_corr)), (xreal(1:size_corr)));
% [a, b] = max(C)




% take FFT of square?
% negative offset of radians?
%

subplot(3,2,1);

stem(real(y(1:2:end)));
title('realy');
ylabel('realy');


subplot(3,2,2);
stem(imag(y(2:2:end)));
title('imagy');

subplot(3, 2, 3);

stem(real(xreal));
title('realx');



subplot(3, 2, 4);

stem((ximag));
title('imagx');

subplot(3, 2, 5);

y1 = y(1:round(length(y)/2));
% y1 = y;
y2 = y(round(length(y)/2):end);
% y2 = y;

fastforT1 = abs(fft(y1.^4));
x_axis1 = linspace(0, 2*pi*(length(y1)-1)/length(y1), length(y1)); %Unclear why we want this.

[max_val1, max_index1] = max(fastforT1);

freq_offset1 = -1*x_axis1(max_index1)./4; %Why divided by 2? ( Prob because we raised to two -Paige)

theta_hat1 = -1*angle(fastforT1(max_index1))./4; %still 2

for k = 1:length(y1)
   x_hat1(k) = y1(k).*exp(1i*(freq_offset1 * (k-1) + theta_hat1)); 
end

plot(x_hat1(round(symbol_period/2):symbol_period:end), '.') %Probably clipping, need to check amplitude
title('half1');


% plot(abs((C)/max(C)));
% title('Correlation');


subplot(3,2,6);

fastforT2 = abs(fft(y2.^4));
x_axis2 = linspace(0, 2*pi*(length(y2)-1)/length(y2), length(y2)); %Unclear why we want this.

[max_val2, max_index2] = max(fastforT2);

freq_offset2 = -1*x_axis2(max_index2)./4; %Why divided by 2? ( Prob because we raised to two -Paige)

theta_hat2 = -1*angle(fastforT2(max_index2))./4; %still 2

for k = 1:length(y2)
   x_hat2(k) = y2(k).*exp(1i*(freq_offset2 * (k-1) + theta_hat2)); 
end

plot(x_hat2(round(symbol_period/2) + mod(length(y1), symbol_period):symbol_period:end), '.') %Probably clipping, need to check amplitude
title('half2');


% plot(x_hat, '.')
% Every 20 to measure in the middle of a pulse
% beginning != end
%Send small noise-like data
%Recive it
%cross correlate it 