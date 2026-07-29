classdef alg_sim_out
    %ALG_SIM_OUT output of an alg_sim simulation routine
    
    properties        
      k; %time index
      z; %output of the operators
      w; %input to the operators
      wp; %performance input
      zp; %performance output 
      u; %controller output/system input
      y; %controller input/system output
      xn; %states of the network
      xc; %states of the controller
      mode; %mode of the switched system   
      param; %parameters used
      f; %function values (if applicable)
      res_w; %optimality error norm(sum(w))
      res_z; %consensus error norm(z - average(z))
      eq; %Equality constraint error norm(Ez - b)    
    end
    
    methods
        function obj = alg_sim_out()
            %ALG_SIM_OUT construct this container for output. This gets filled by alg_sim.            
        end
        

    end
end

