% Phase lock loop: Shove the "fake" stuff back into the "real" stuff

clear;
% Open the file containing the received samples
f2 = fopen('rx.dat', 'rb');
% read data from the file
rxfile = fread(f2, 'float32');
% close the file
fclose(f2);

f1 = fopen('tx2.dat', 'rb');
txfile = fread(f1, 'float32');
fclose(f1);
txfile = txfile*100;
xreal = txfile(1:2:end);
ximag = 1i*txfile(2:2:end);


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

start = -1;
runningsum = zeros(1, 100);
for z = 1:2:length(y)
    runningsum(mod((z-1)/2, 10) + 1) = abs(y(z));
    if (sum(runningsum)/length(runningsum)) > maxreal*3
       start = z - length(runningsum)
       break
    end
end

ends = -1;
runningsum = zeros(1, 100);
for z = length(y):-2:start
     runningsum(mod((z-1)/2, 10) + 1) = abs(y(z));
    if (sum(runningsum)/length(runningsum)) > maxreal*3
        ends = z + length(runningsum)
        break
    end
end

start_constant = 0; % 200000;
end_constant = 0; %200000;
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

stem(real(y(2:2:end)));
title('realy');
ylabel('realy');
subplot(3,2,2);

stem(imag(y(1:2:end)));
title('imagy');
subplot(3, 2, 3);

stem(real(xreal));
title('realx');



subplot(3, 2, 4);

stem((ximag));
title('imagx');

subplot(3, 2, 5);


% plot(abs((C)/max(C)));
% title('Correlation');


subplot(3,2,6);

fastforT = abs(fft(y.^2));
x_axis = linspace(0, 2*pi*(length(y)-1)/length(y), length(y)); %Unclear why we want this.

[max_val, max_index] = max(fastforT);

freq_offset = -1*x_axis(max_index)./2; %Why divided by 2? Prob because we raised to two

theta_hat = -1*angle(fastforT(max_index))./2; %still 2

for k = 1:length(y)
   x_hat(k) = y(k).*exp(1i*(freq_offset * (k-1) + theta_hat)); 
end

plot(x_hat, '.') %Probably clipping, need to check amplitude
%Send small noise-like data
%Recive it
%cross correlate it 