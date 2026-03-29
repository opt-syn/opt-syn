function [model_sys_opt] = wrap_optima_model(regulator)
%WRAP_OPTIMA_MODEL Add an extra dead state for optimality
%this is vestigial. 



s = size(regulator.Gam, 1);

model_sys = regulator.model_sys;

%add the optima state (to be stabilized to 0)
[AM, BM, CM, DM] = ssdata(model_sys);
AMO = blkdiag(0, AM);
BMO = [ones(1, s), zeros(1, 2*s); BM];
CMO = [zeros(2*s, 1), CM];
DMO = DM;
model_sys_opt = ss(AMO, BMO, CMO, DMO, 1);
model_sys_opt.StateName{1} = 'opt';
model_sys_opt.StateName(2:end) = model_sys.StateName;
model_sys_opt.InputName = model_sys.InputName;
model_sys_opt.OutputName = model_sys.OutputName;

end

