clear;
% Open the file containing the received samples
f2 = fopen('rx.dat', 'rb');
% read data from the file
rxfile = fread(f2, 'float32');
% close the file
fclose(f2);

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


f1 = fopen('tx2.dat', 'rb');
txfile = fread(f1, 'float32');
fclose(f1);
x = txfile;



subplot(3,2,1);
title("realy");
stem(real(y(z:2:end)));
subplot(3,2,2);
title("imagy");
stem(imag(y((z-1):2:end)));
subplot(3, 2, 3);
title("realx");
stem(imag(x));
subplot(3, 2, 4);
title("imagx");
stem(real(x));
subplot(3, 2, 5);

C = xcorr(y(1:2:end), x(1:2:end));
plot((C)/max(C));
ylabel('Correlation');