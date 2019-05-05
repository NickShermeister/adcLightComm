function binaryvec = str2bin(message, n)

% message is a string
% n sets the binary base
if nargin<2 % if not specified, n defaults to 8 
  n = 8;
end

bits = dec2bin(message,n);
binaryvec = reshape(bits.'-'0',1,[]);

end