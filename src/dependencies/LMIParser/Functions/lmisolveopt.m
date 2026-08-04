function opt=lmisolveopt;
reltol=1e-5;
maxiter=300;
feasradius=-1e7; %flexible bound method
L=20;
verbose=0;
opt=[reltol maxiter feasradius L verbose];
end
