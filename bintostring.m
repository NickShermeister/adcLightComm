function str = bintostring(bin_vec)
    str = char(bin2dec(reshape(char(bin_vec+'0'), 8,[]).')');
end