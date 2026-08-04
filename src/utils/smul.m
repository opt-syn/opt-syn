function so = smul(s, d)
%smul scalar multiplication times lmim object

[i, j] = find(s);

so = lmim(s);
for k = 1:length(i)    
    so(i(k), j(k)) = so(i(k), j(k)) +  s(i(k), j(k)) * d;
end

end

