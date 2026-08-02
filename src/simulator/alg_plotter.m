classdef alg_plotter
    %ALG_PLOTTER Plot trajectories of an algorithm execution
    
    
    properties
        sim_out;  %data in the simulation        
        FS = 16;  %fontsize for axis labels
        FST = 20; %fontsize for title  
        EQUALITY = false; %is an equality constraint used? 
        visible=true; % should the plot be shown?
    end
    
    methods
        function obj = alg_plotter(sim_out)
            %ALG_PLOTTER Construct a plotter for a trajectory
            %
            % Args:
            %   sim_out (alg_sim_out): the trajectory from alg_sim

            obj.sim_out = sim_out;
            obj.EQUALITY = ~isempty(sim_out.eq);            
        end
        
        function fig = plot(obj,traces, fignum)
            %PLOT multi-pane display 
            %
            %Args:
            %   traces: signals to plot (e.g. {'x', 'w', 'f'})
            %   fignum: figure number to display
            %Returns:
            %   fig:    figure environment
                        

            if obj.visible
                if (nargin == 3) && isnumeric(fignum) && ~isempty(fignum)                
                    fig = figure(fignum);
                else
                    fig = figure();                    
                end
            else
                fig = figure('visibility', 'off');
            end               
            clf(fig, 'reset');

            nplt = length(traces);
            if nplt ==3
                nrows = 1;
                ncols = 3;
            else
                nrows = floor(sqrt(nplt));
                ncols = ceil(sqrt(nplt));
            end
            tiledlayout(nrows, ncols)
            axlist = cell(nplt, 1);
            for r = 1:nplt
                axcurr = nexttile;
                axcurr = obj.plot_tile(axcurr, traces{r});
                axlist{r} = axcurr;
            end
        end

        function fig = plot_6f(obj, fignum)
            %PLOT_6f plot the signals ('x', 'w', 'res_w', 'f', 'z','res_z')            
            %Args:            
            %   fignum: figure number to display
            %Return:
            %   fig:    figure environment


            if nargin < 2
                fignum = [];
            end
            sigs = {'x', 'w', 'res_w', 'f', 'z', 'res_z'};
            fig = obj.plot(sigs, fignum);
        end

        function fig = plot_6(obj, fignum)
            %PLOT_6 plot the signals ('xn', 'w', 'res_w', 'xc', 'z','res_z')            
            %Args:            
            %   fignum: figure number to display
            %Return:
            %   fig:    figure environment

            
            if nargin < 2
                fignum = [];
            end
            sigs = {'xn', 'w', 'res_w', 'xc', 'z', 'res_z'};
            fig = obj.plot(sigs, fignum);
        end

        function fig = plot_4(obj, fignum)
            %PLOT_4 plot the signals ('w', 'res_w', 'z','res_z')            
            %Args:            
            %   fignum: figure number to display
            %Return:
            %   fig:    figure environment

            
            if nargin < 2
                fignum = [];
            end
            sigs = {'w', 'res_w', 'z', 'res_z'};
            fig = obj.plot(sigs, fignum);
        end


        function fig = plot_4_err(obj, fignum)
            %PLOT_4_err plot the error signals/regulated quantities ('xnerr', 'uerr', 'xcerr', 'yerr')            
            %Args:            
            %   fignum: figure number to display
            %Return:
            %   fig:    figure environment

            
            if nargin < 2
                fignum = [];
            end
            sigs = {'xnerr', 'uerr', 'xcerr', 'yerr'};
            fig = obj.plot(sigs, fignum);
        end

        function fig = plot_3_err(obj, fignum)
            %PLOT_3_err plot the error signals/regulated quantities ('xerr', 'uerr', 'yerr')            
            %Args:            
            %   fignum: figure number to display
            %Return:
            %   fig:    figure environment
            if nargin < 2
                fignum = [];
            end
            sigs = {'xerr','uerr', 'yerr'};
            fig = obj.plot(sigs, fignum);
        end

        function fig = plot_3_sq_err(obj, fignum)
            %PLOT_3_sq_err plot the squared error signals/regulated quantities ('xerr', 'uerr', 'yerr')            
            %Args:            
            %   fignum: figure number to display
            %Return:
            %   fig:    figure environment
            if nargin < 2
                fignum = [];
            end
            sigs = {'sq_xerr','sq_uerr', 'sq_yerr'};
            fig = obj.plot(sigs, fignum);
        end

        function fig = plot_4_sq_err(obj, fignum)
            %PLOT_4_sq_err plot the squared error signals/regulated quantities ('sq_xnerr', 'sq_uerr', 'sq_xcerr', 'sq_yerr')            
            %Args:            
            %   fignum: figure number to display
            %Return:
            %   fig:    figure environment
            if nargin < 2
                fignum = [];
            end
            sigs = {'sq_xnerr', 'sq_uerr', 'sq_xcerr', 'sq_yerr'};
            fig = obj.plot(sigs, fignum);
        end


        function obj = add_opt_sig(obj, reg_cl, dstar)
            %ADD_OPT_SIG Use the optimal trajectory to define the error signals            
            %Args:            
            %   reg_cl: closed-loop regulator equation, output from regulator.check_regulator()            
            %   dstar:  properties of optimal solution :math:`(-\beta^*, \hat{w}^*)`
            %Return:
            %   obj: the plotter
            
            obj.sim_out.xnerr = obj.sim_out.xn -  reg_cl.Pi * dstar;
            obj.sim_out.xcerr = obj.sim_out.xc  - reg_cl.Th * dstar;
            obj.sim_out.yerr = obj.sim_out.y - reg_cl.Phi * dstar;
            obj.sim_out.uerr = obj.sim_out.u- reg_cl.Gam * dstar;            

            fs = @(sig) squeeze(sum(sig.^2, [1, 2]));

            obj.sim_out.sq_xnerr = fs(obj.sim_out.xnerr);
            obj.sim_out.sq_xcerr = fs(obj.sim_out.xcerr);
            obj.sim_out.sq_yerr = fs(obj.sim_out.yerr);
            obj.sim_out.sq_uerr = fs(obj.sim_out.uerr);            


        end

        function ax = plot_tile(obj, ax, sig)
            %PLOT_TILE plot the signal 'sig' v.s. time
            %Args:
            %   ax: axis object in plot
            %   sig: the signal to use
            %Returns
            %   ax: axis with the signal

            k = obj.sim_out.k;
            T = length(k);
            if ismember(sig, fieldnames(obj.sim_out))
                sig_curr = getfield(obj.sim_out, sig);
            elseif strcmp(sig, 'x')
                sig_curr = [obj.sim_out.xn; obj.sim_out.xc];
            elseif strcmp(sig, 'payoff')
                sig_curr = [obj.sim_out.f];
            elseif strcmp(sig, 'xerr')
                sig_curr = [obj.sim_out.xnerr; obj.sim_out.xcerr];
            elseif strcmp(sig, 'sq_xerr')
                sig_curr = obj.sim_out.sq_xnerr +  obj.sim_out.sq_xcerr;
            elseif strcmp(sig, 'delay')
                sig_curr = obj.sim_out.mode - 1;
            elseif strcmp(sig, 'coord')
                sig_curr = obj.sim_out.mode;
            end
               

                sz_curr = size(sig_curr);
                sig_flat = reshape(permute(sig_curr, [length(sz_curr), 1:(length(sz_curr)-1)]), T,  []);
    
                
                szflat = size(sig_flat, 2);
                hold on
                for i =1:szflat
                    sig_plot = sig_flat(:, i);
                    

                    plot(obj.sim_out.k, sig_plot)
                end
    
                xlabel('$k$', 'interpreter', 'latex', 'fontsize', obj.FS)
                ylabel(obj.get_name(sig), 'interpreter', 'latex', 'fontsize', obj.FS)
                title(obj.get_title(sig), 'interpreter', 'latex', 'fontsize', obj.FST)
                if ismember(sig, {'res_w', 'res_z', 'sq_xcerr', 'sq_xerr', ...
                        'sq_xnerr', 'sq_uerr', 'sq_yerr'})
                    set(ax, 'YScale', 'log');
                end
                xlim([k(1), k(end)])
                
                if strcmp(sig, 'delay') || strcmp(sig, 'mode') || strcmp(sig, 'coord')
                    yticks(1:max(sig_plot));                    
                end
            
        end

        function name = get_title(obj, sig)
            %GET_TITLE get the title for the plot
            %Args:
            %   sig: the signal that is plotted
            %Return:
            %   name: the title to use
            switch sig
                case 'f' 
                    name = 'Function Value';
                case 'x'
                    name = 'State';
                case 'z'
                    name = 'Iterate';
                case 'w'
                    name = 'Oracle Output';
                case 'wp'
                    name = 'Performance Input';
                case 'zp'
                    name = 'Performance Output';
                case 'mode'
                    name = 'Switching Mode';
                case 'delay'
                    name = 'Time Delay';
                case 'coord'
                    name = 'Coordinate Block';
                case 'res_w'
                    name = 'Optimality Error';
                case 'res_z'
                    name = 'Consensus Error';
                case 'xc'
                    name = 'State (Controller)';
                case 'xn'
                    name = 'State (Network)';
                case 'eq'
                    name = 'Primal Feasibility';
                case 'payoff'
                    name = 'Payoff';
                case 'xerr'
                    name = 'State Error';
                case 'xcerr'
                    name = 'State (Controller) Error';
                case 'xnerr'
                    name = 'State (Network) Error';
                case 'uerr'
                    name = 'Input Error';
                case 'yerr'
                    name = 'Output Error';

                %squared errors
                case 'sq_xerr'
                    name = 'State Error';
                case 'sq_xcerr'
                    name = 'State (Controller) Error';
                case 'sq_xnerr'
                    name = 'State (Network) Error';
                case 'sq_uerr'
                    name = 'Input Error';
                case 'sq_yerr'
                    name = 'Output Error';

            end
        end

        function name = get_name(obj, sig)
            %GET_NAME the name of the signal in latex-formatted strings
            %for use in y-axis labels
            %
            %Args:
            %   sig: the signal that is plotted
            %Return:
            %   name: the title to use

            %specific channels
            if length(sig)==1
                name_mid = ['$', sig, '$'];
            elseif strcmp(sig(2), '_')
                name_mid = ['$', sig(1), '_{', sig(3:end), '}$'];
            elseif strcmp(sig, 'wp')
                name_mid = '$w_p$';
            elseif strcmp(sig, 'zp')
                name_mid = '$z_p$';
            elseif strcmp(sig, 'xc')
                name_mid = '$x_c$';
            elseif strcmp(sig, 'xn')
                name_mid = '$x_{N}$';
            elseif strcmp(sig, 'eq')
                name_mid = '$|| E z - b||_2$';

            %switching
            elseif strcmp(sig, 'mode')
                name_mid = 'mode';
            elseif strcmp(sig, 'delay')
                name_mid = 'delay';
            elseif strcmp(sig, 'coord')
                name_mid = 'coord.';
            elseif strcmp(sig, 'payoff')
                name_mid = '$f$';

            %tracking    
            elseif strcmp(sig, 'xcerr')
                name_mid = '$x_c - x_c^*$';
            elseif strcmp(sig, 'xerr')
                name_mid = '$x - x^*$';
            elseif strcmp(sig, 'xnerr')
                name_mid = '$x_N - x_N^*$';
            elseif strcmp(sig, 'uerr')
                name_mid = '$u - u^*$';
            elseif strcmp(sig, 'yerr')
                name_mid = '$y - y^*$';

            %tracking residuals
            elseif strcmp(sig, 'sq_xcerr')
                name_mid = '$||x_c - x_c^*||^2_2$';
            elseif strcmp(sig, 'sq_xerr')
                name_mid = '$||x - x^*||^2_2$';
            elseif strcmp(sig, 'sq_xnerr')
                name_mid = '$||x_N - x_N^*||^2_2$';
            elseif strcmp(sig, 'sq_uerr')
                name_mid = '$||u - u^*||^2_2$';
            elseif strcmp(sig, 'sq_yerr')
                name_mid = '$||y - y^*||^2_2$';


            %residuals
            elseif strcmp(sig, 'res_w')
                if obj.EQUALITY
                    % name_mid = '$||\text{Proj}_{\text{null} \ E} (1^{\top} w)||_2$';
                    name_mid = '$||$ Proj $_{null(E)} [1^{\top} w] ||_2$';
                else
                    name_mid = '$||1^{\top} w||_2$';
                end
            elseif strcmp(sig, 'res_z')
                name_mid = '$||z - z_{avg}||_2$';
            end

            name = name_mid;            
        end
    end
end

