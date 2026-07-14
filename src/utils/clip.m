function xclip =  clip(x, lower, upper)

    c = size(lower, 1);
    x = reshape(x, c, []);
    sl = x <= lower;
    su = x >= upper;

    xclip = x;
    for i = 1:c
        xclip(sl(i, :)) = lower(i);
        xclip(su(i, :)) = upper(i);
    end
    
    xclip = reshape(xclip, [], 1);
end