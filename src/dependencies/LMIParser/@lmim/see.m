function see(pi)
disp(' ')
for j=1:length(pi);
    p=pi(j);
    %long display of array of lmim objects
    disp('Matrix A')
    disp(p.A)
    disp('Matrix B')
    disp(p.B)
    disp('Matrix C')
    disp(p.C)
    mydisp(p)
    disp('-------------------------------------------------------------------------------')
end
end