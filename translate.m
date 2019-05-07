function words = translate(signalVec)

binVec = ones(length(signalVec),1);
binVec(find(signalVec == -1)) = 0;


words = bin2str(binVec);

end
