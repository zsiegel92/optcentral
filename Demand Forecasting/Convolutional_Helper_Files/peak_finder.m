function  [I,pks]  = peak_finder(range, domain )
[pks,I]=findpeaks(range);
[pks,inds]=sort(pks);
pks = pks(end:-1:1);
I = domain(I(inds(end:-1:1)));
end

