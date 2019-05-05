function str = bin2str(bin_vec,n)
% bin_vec is a vector of binary digits
% n sets the binary base
    if nargin<2 % if not specified, n defaults to 8
        n = 8;
    end
    str = char(bin2dec(reshape(char(bin_vec+'0'), n,[]).')');
end