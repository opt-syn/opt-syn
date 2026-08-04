function [n_out] = sys_compare(sys_1, sys_2)
%SYS_COMPARE use the elementwise norm to compare two systems
%only use for numerics, not for system responses
n_out = [norm(sys_1.A - sys_2.A), norm(sys_1.A - sys_2.A);
    norm(sys_1.C - sys_2.C), norm(sys_1.D - sys_2.D)];
end

