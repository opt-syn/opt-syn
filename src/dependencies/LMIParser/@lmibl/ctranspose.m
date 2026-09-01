function s=ctranspose(p)
% function s=ctranspose(p)
%
% Transpose lmiblk object p
% X to X'.

if ~isscalar(p) %as an array of lmiblk objects
    error('Transposition of array of lmiblk objects not defined.')
end
s=p;
s.di=[p.di(2) p.di(1)];
s.tr=not(p.tr);