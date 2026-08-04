function sys=rhotrafo(sy,rho);
%rho weighting for state-space systems
%transforms [A,B,C,D] into [1/rho*A,1/rho*B,C,D]

if rho<=0;
    error('rho must be positive');
else
    rhoi=1/rho;
    % if isa(sy, 'genss')
        sys=ss(rhoi*sy.A,rhoi*sy.B,sy.C,sy.D, 1);
    % else
        % [A,B,C,D]=ssdata(sy);
        % sys=ss(rhoi*A,rhoi*B,C,D,sy.Ts);
    % end
end
