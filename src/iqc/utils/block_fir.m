function [A, B] = block_fir(n)
%block_FIR jordan block for FIR system

if n==0
    A = [];
    B = [];
else
    A = [0, zeros(1, n-1);
        eye(n-1), zeros(n-1, 1)];
    B = [1; zeros(n-1, 1)];
end

end

