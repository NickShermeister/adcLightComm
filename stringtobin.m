function binaryvec = stringtobin(message)

bits = dec2bin(message,8);
binaryvec = reshape(bits.'-'0',1,[]);

end