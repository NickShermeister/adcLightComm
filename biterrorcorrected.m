% Phase lock loop: Shove the "fake" stuff back into the "real" stuff

clear;
clf;


N = 1000;
% make N random bits of values +- 1
seed = 562019;
rng(seed);
constant_bits = sign(randn(N,1)) + 1i*sign(randn(N,1));

%Set the symbol period
symbol_period = 20;


% Open the file containing the received samples
f2 = fopen('rxhello.dat', 'rb');
% read data from the file
rxfile = fread(f2, 'float32');
% close the file
fclose(f2);

%Assign y to be the data, ignoring the first 250 values due to weird noise.
y = zeros(length(rxfile)/2,1);
y = rxfile(1:2:end)+1i*rxfile(2:2:end);
y = y(250:end);

%open file and assign real/imag parts
f1 = fopen('txhello.dat', 'rb');
txfile = fread(f1, 'float32');
fclose(f1);
txfile = txfile*100;
xreal = txfile(1:2:end);
ximag = txfile(2:2:end);

%Trim down the data to the relevent section
y = movingAvg(y);
y = crossCorr(y, constant_bits);

%Supposed to be able to normalize the values.
magnitude_estimate = rms(abs(y));
y = y./magnitude_estimate;

%Hard set: change the tx data based on the known constants in the file
%making the tx data.
xreal = xreal(100000:(length(xreal)-100000));
ximag = ximag(100000:(length(ximag)-100000));

%Plot the received and transmitted data
subplot(3,2,1);

stem(real(y(1:2:end)));
title('Real portion of received Data');
ylabel('realy');

subplot(3,2,2);
stem(imag(y(2:2:end)));
title('Imaginary portion of received Data');

subplot(3, 2, 3);
stem(real(xreal));
title('Real portion of transmitted Data');

subplot(3, 2, 4);
stem((ximag));
title('Imaginary portion of transmitted Data');

% %normalize y
% temp_max = max(abs(y));
% y = y./temp_max;

%Divide the data in half because of the noisiness and because the clocks
%will drift over time.
y1 = y(1:round(length(y)/2));
y1 = y;
y2 = y(round(length(y)/2):end);

fastforT1 = abs(fft(y1.^4));
fastforT1(1) = 0;
x_axis1 = linspace(0, 2*pi*(length(y1)-1)/length(y1), length(y1)); %Unclear why we want this.
[max_val1, max_index1] = max(fastforT1);
freq_offset1 = -1*x_axis1(max_index1)./4; %Why divided by 2? ( Prob because we raised to two -Paige)
theta_hat1 = -1*angle(fastforT1(max_index1))./4; %still 2
% theta_hat1 = theta_hat1 + 1*pi/4;

for k = 1:length(y1)
   x_hat1(k) = y1(k).*exp(1i*(freq_offset1 * (k-1) + theta_hat1)); 
%    x_hat1(k) = y1(k).*exp(1i*(theta_hat1)); 

end

subplot(3, 2, 5);
% x_hat1 = y;

f1 = fopen('paige.dat', 'wb');
% write the values as a float32

fwrite(f1, x_hat1, 'float32');
% close the file
fclose(f1);

down_x = downsample(x_hat1(10:end), 20);
% plot(x_hat1(round(symbol_period/2):symbol_period:end), '.') %Probably clipping, need to check amplitude
plot(real(down_x), imag(down_x), '.');
title('half1');
axis square


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
axis square

% y = x_hat1';
% y = round(y);

%Trim off zeros at the front and set y's length equal to the sent data's
%length (possibly a bad presumption, but it should align technically)
test_loc = 1;
while(y(test_loc) == 0)
   test_loc = test_loc + 1;
end
%TODO: fix to make it actually the length of ximag
y = y((test_loc):(length(ximag) + test_loc-2 ));
complex = (1i*ximag + xreal);
complex = complex(2:end);


%HAND TUNE THESE
y_sum1 = (0 < real(y));
y_sum2 = (0 > real(y));
y_sum3 = (0 < imag(y));
y_sum4 = (0 > imag(y));
y_total1 = y_sum1*1 + y_sum2*2;
y_total2 = y_sum3*1 + y_sum4 * 2;

x_sum1 = (0 < real(complex));
x_sum2 = (0 > real(complex));
x_sum3 = (0 < imag(complex));
x_sum4 = (0 > imag(complex));
x_total1 = x_sum1*1 + x_sum2*2;
x_total2 = x_sum3*1 + x_sum4*2;

diffTest1 = (y_total1((symbol_period/2):symbol_period:end) == x_total1((symbol_period/2):symbol_period:end));
diffTest2 = (y_total2((symbol_period/2):symbol_period:end) == x_total2((symbol_period/2):symbol_period:end));

accuracy1 = sum(diffTest1/length(diffTest1))
accuracy2 = sum(diffTest2/length(diffTest2))

%Change the Rx basd on the values:
% final = y_total((symbol_period/2):symbol_period:end);
% for x = 1:length(final)
%    if final(x) == 4
%        final(x) = (-1 + 1i);
%    elseif final(x) == 1
%        final(x) = (1 + 1i);
%    elseif final(x) == 3
%        final(x) = (-1 - 1i);
%    elseif final(x) == 2
%        final(x) = (1 - 1i);
%    else
%        final(x) = 0;
%    end
% end
% final = final(1:(size(final)-1));
% 
% 
% temp_sum = sum(final==0);
% 
% temp_fin = zeros((length(final) - temp_sum-1)*2, 1);
% curr_loc = 1;
% for x = 1:1:length(final)
%    if x ~= 0
%        temp_fin(curr_loc) = real(final(x));
%        temp_fin(curr_loc+1) = imag(final(x));
%        curr_loc = curr_loc + 2;
%    end
% end

% final = temp_fin;


% Every 20 to measure in the middle of a pulse 
%Send small noise-like data
%Recive it
%cross correlate it 

function z = movingAvg(y)
    %basic filtering
tempreal = y(1:2:200);
tempimag = y(2:2:200);
maxreal = max(abs(tempreal));
maximag = max(tempimag);
start_constant = 0; % 200000;
end_constant = 0; %200000;
start = 1;
runningsum = zeros(1, 1000);
for z = 1:2:length(y)
    runningsum(mod((z-1)/2, 500) + 1) = abs(y(z));
    if (sum(runningsum)/length(runningsum)) > maxreal*3
       start = z - length(runningsum)
       break
    end
end

ends = length(y) - end_constant;
runningsum = zeros(1, 100);
for z = length(y):-2:(length(runningsum)+2)
%     abs(y(z))
    runningsum(mod(round((z-1)/2), 100) + 1) = abs(y(z));
    if (sum(runningsum)/length(runningsum)) > maxreal*3
        ends = z + length(runningsum)
        break
    end
end

z = y((start - start_constant):(ends+ end_constant));
end

function z = crossCorr(y, bitsIn)
   C = xcorr(y, bitsIn, 200);
   C = C/max(C);
   [maxval, maxindex] = max(C)
   
   z = y(maxindex:end);
end