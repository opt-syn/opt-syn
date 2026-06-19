classdef alg_plotter
    %ALG_PLOTTER Plot trajectories of an algorithm execution
    
    
    properties
        sim;      %data in the simulation
        FS = 16;  %fontsize for axes
        FST = 20; %fontsize for title  
        EQUALITY;
        opt_sig; 
    end
    
    methods
        function obj = alg_plotter(sim)
            %ALG_PLOTTER Construct an instance of this class
            %   Detailed explanation goes here
            obj.sim = sim;
            obj.EQUALITY = ~isempty(sim.eq);            
        end
        
        function fig = plot(obj,traces, fignum)
            %PLOT: multi-pane display 
            %Input:
            %   traces: signals to plot (e.g. {'x', 'w', 'f'})
                        
            if (nargin == 3) && isnumeric(fignum) && ~isempty(fignum)
                fig = figure(fignum);
            else
                fig = figure;
            end
            clf

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
            if nargin < 2
                fignum = [];
            end
            sigs = {'x', 'w', 'res_w', 'f', 'z', 'res_z'};
            fig = obj.plot(sigs, fignum);
        end

        function fig = plot_6(obj, fignum)
            if nargin < 2
                fignum = [];
            end
            sigs = {'xn', 'w', 'res_w', 'xi', 'z', 'res_z'};
            fig = obj.plot(sigs, fignum);
        end

        function fig = plot_4(obj, fignum)
            if nargin < 2
                fignum = [];
            end
            sigs = {'w', 'res_w', 'z', 'res_z'};
            fig = obj.plot(sigs, fignum);
        end


        function fig = plot_4_err(obj, fignum)
            if nargin < 2
                fignum = [];
            end
            sigs = {'xnerr', 'uerr', 'xierr', 'yerr'};
            fig = obj.plot(sigs, fignum);
        end

        function fig = plot_3_err(obj, fignum)
            if nargin < 2
                fignum = [];
            end
            sigs = {'xerr','uerr', 'yerr'};
            fig = obj.plot(sigs, fignum);
        end

        function fig = plot_4_sq_err(obj, fignum)
            if nargin < 2
                fignum = [];
            end
            sigs = {'sq_xnerr', 'sq_uerr', 'sq_xierr', 'sq_yerr'};
            fig = obj.plot(sigs, fignum);
        end

        function fig = plot_3_sq_err(obj, fignum)
            if nargin < 2
                fignum = [];
            end
            sigs = {'sq_xerr','sq_uerr', 'sq_yerr'};
            fig = obj.plot(sigs, fignum);
        end

        function obj = add_opt_sig(obj, reg_cl, dstar)
            %add the optimal trajectory
            
            obj.sim.xnerr = obj.sim.xn -  reg_cl.Pi * dstar;
            obj.sim.xierr = obj.sim.xi  - reg_cl.Th * dstar;
            obj.sim.yerr = obj.sim.y - reg_cl.Phi * dstar;
            obj.sim.uerr = obj.sim.u- reg_cl.Gam * dstar;            

            fs = @(sig) squeeze(sum(sig.^2, [1, 2]));

            obj.sim.sq_xnerr = fs(obj.sim.xnerr);
            obj.sim.sq_xierr = fs(obj.sim.xierr);
            obj.sim.sq_yerr = fs(obj.sim.yerr);
            obj.sim.sq_uerr = fs(obj.sim.uerr);            


        end

        function ax = plot_tile(obj, ax, sig)
            %PLOT_TILE plot the signal 'sig' v.s. time


            k = obj.sim.k;
            T = length(k);
            if ismember(sig, fieldnames(obj.sim))
                sig_curr = getfield(obj.sim, sig);
            elseif strcmp(sig, 'x')
                sig_curr = [obj.sim.xn; obj.sim.xi];
            elseif strcmp(sig, 'xerr')
                sig_curr = [obj.sim.xnerr; obj.sim.xierr];
            elseif strcmp(sig, 'sq_xerr')
                sig_curr = obj.sim.sq_xnerr +  obj.sim.sq_xierr;
            elseif strcmp(sig, 'delay')
                sig_curr = obj.sim.mode - 1;
            end
               
                %TODO: plot the regulation signals

                sz_curr = size(sig_curr);
                sig_flat = reshape(permute(sig_curr, [length(sz_curr), 1:(length(sz_curr)-1)]), T,  []);
    
                
                szflat = size(sig_flat, 2);
                hold on
                for i =1:szflat
                    sig_plot = sig_flat(:, i);
                    

                    plot(obj.sim.k, sig_plot)
                end
    
                xlabel('$k$', 'interpreter', 'latex', 'fontsize', obj.FS)
                ylabel(obj.get_name(sig), 'interpreter', 'latex', 'fontsize', obj.FS)
                title(obj.get_title(sig), 'interpreter', 'latex', 'fontsize', obj.FST)
                if ismember(sig, {'res_w', 'res_z', 'sq_xierr', 'sq_xerr', ...
                        'sq_xnerr', 'sq_uerr', 'sq_yerr'})
                    set(ax, 'YScale', 'log');
                end
                xlim([k(1), k(end)])
                
                if strcmp(sig, 'delay')
                    yticks(1:max(sig_plot));                    
                end
            
        end

        function name = get_title(obj, sig)
            %GET_TITLE get the title for the plot
            switch sig
                case 'f' 
                    name = 'Function Value';
                case 'x'
                    name = 'State';
                case 'z'
                    name = 'Iterate';
                case 'w'
                    name = 'Oracle Output';
                case 'w_p'
                    name = 'Performance Input';
                case 'z_p'
                    name = 'Performance Output';
                case 'mode'
                    name = 'Switching Mode';
                case 'delay'
                    name = 'Time Delay';
                case 'res_w'
                    name = 'Optimality Error';
                case 'res_z'
                    name = 'Consensus Error';
                case 'xi'
                    name = 'State (Controller)';
                case 'xn'
                    name = 'State (Network)';
                case 'eq'
                    name = 'Primal Feasibility';
                case 'xerr'
                    name = 'State Error';
                case 'xierr'
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
                case 'sq_xierr'
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

            %specific channels
            if length(sig)==1
                name_mid = ['$', sig, '$'];
            elseif strcmp(sig(2), '_')
                name_mid = ['$', sig(1), '_{', sig(3:end), '}$'];
            elseif strcmp(sig, 'wp')
                name_mid = '$w_p$';
            elseif strcmp(sig, 'zp')
                name_mid = '$z_p$';
            elseif strcmp(sig, 'xi')
                name_mid = '$\xi$';
            elseif strcmp(sig, 'xn')
                name_mid = '$x_{N}$';
            elseif strcmp(sig, 'eq')
                name_mid = '$|| E z - b||_2$';

            %switching
            elseif strcmp(sig, 'mode')
                name_mid = 'mode';
            elseif strcmp(sig, 'delay')
                name_mid = 'delay';

            %tracking    
            elseif strcmp(sig, 'xierr')
                name_mid = '$\xi - \xi^*$';
            elseif strcmp(sig, 'xerr')
                name_mid = '$x - x^*$';
            elseif strcmp(sig, 'xnerr')
                name_mid = '$x_N - x_N^*$';
            elseif strcmp(sig, 'uerr')
                name_mid = '$u - u^*$';
            elseif strcmp(sig, 'yerr')
                name_mid = '$y - y^*$';

            %tracking residuals
            elseif strcmp(sig, 'sq_xierr')
                name_mid = '$||\xi - \xi^*||^2_2$';
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
            % name = ['$', name_mid, '$'];
        end
    end
end

