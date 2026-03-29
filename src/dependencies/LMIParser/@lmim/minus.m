function s=minus(s1,s2)
% function s=minus(s1,s2)
%
% Differnce s1-s2 of two lmim objects s1 and s2.
s1=lmim(s1);
s2=-lmim(s2);
s=s1+s2;
end
