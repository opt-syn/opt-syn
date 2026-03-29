function c =  clip(x, lower, upper)
    sl = x <= lower;
    su = x >= upper;

    c = x;
    c(sl) = lower;
    c(su) = upper;

end