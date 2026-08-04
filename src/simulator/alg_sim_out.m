classdef alg_sim_out
    %ALG_SIM_OUT output of an alg_sim simulation routine
    
    properties        
      s; %number of oracles (in bind)
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


      %regulator equation tracking
      xnerr; %error in state of network
      xcerr; %error in state of network
      yerr; %error in controller input
      uerr; %error in controller output

      sq_xnerr; %square of error in state of network
      sq_xcerr; %error in state of network
      sq_yerr; %error in controller input
      sq_uerr; %error in controller output
    end
    
    methods
        function obj = alg_sim_out()
            %ALG_SIM_OUT construct this container for output. This gets filled by alg_sim.            
        end
        

    end
end

