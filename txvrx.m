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
x = txfile;


y = zeros(length(rxfile)/2,1);
y = rxfile(1:2:end)+j*rxfile(2:2:end);
y = y(250:end);

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


x = x(200000:(length(x)-200000));

size_corr = 1000;
C = xcorr((y(1:size_corr)), (x(1:size_corr)));



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

stem(imag(x));
title('imagx');

subplot(3, 2, 4);

stem(real(x));
title('realx');

subplot(3, 2, 5);


plot(abs((C)/max(C)));
title('Correlation');

%Send small noise-like data
%Recive it
%cross correlate it 