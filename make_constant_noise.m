% tx_samples_from_file --freq 250e3 --rate 200e3 --type float --gain 70 --file tx2.dat
%Lower receiver gain to like 20.
clear;

fileName = 'constant_noise.mat';


N = 1000;
% make N random bits of values +- 1
seed = 562019;
rng(seed);
a = sign(randn(N,1));
% rng(seed + 1);
% b = sign(randn(N,1));
constant_bits = a - 1i*a;


% open a file to write in binary format 
save(fileName);