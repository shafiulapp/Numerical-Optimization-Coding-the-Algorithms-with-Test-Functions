classdef plotsum < handle
    properties
        fvals = [];          % Stores objective function values
        gnorm = [];          % Stores gradient norms
        step_lengths = [];   % Stores step lengths
        timings = [];        % Stores timing for each step
    end

    methods
        % Add data for a single step
        function add_data(obj, fval, grad_norm, step_length, timing)
            obj.fvals(end + 1) = fval;
            obj.gnorm(end + 1) = grad_norm;
            obj.step_lengths(end + 1) = step_length;
            obj.timings(end + 1) = timing;
        end

        % Plot performance metrics
        function plot_metrics(obj)
            if isempty(obj.fvals)
                disp('No data to plot.');
                return;
            end

            obj.plot_objective_function();

            obj.plot_gradient_norm();
            obj.plot_step_lengths();


            
        end

        function plot_step_lengths(obj)
            if isempty(obj.step_lengths)
                disp('No step length data to plot.');
                return;
            end

            % Create a new figure
            figure;
            yyaxis left;
            plot(obj.step_lengths, 'LineWidth', 2); % Linear scale plot
            ylabel('Linear Scale');
            xlabel('Iterations');

            yyaxis right;
            semilogy(obj.step_lengths, 'LineWidth', 2); % Logarithmic scale plot
            ylabel('Log Scale');

            title('Step Lengths');
            grid on;
        end

        function plot_objective_function(obj)
            if isempty(obj.fvals)
                disp('No objective function data to plot.');
                return;
            end

            figure;
            yyaxis left;
            plot(obj.fvals, 'LineWidth', 2); 
            ylabel('Linear Scale');
            xlabel('Iterations');

            yyaxis right;
            semilogy(obj.fvals, 'LineWidth', 2);
            ylabel('Log Scale');

            title('Objective Function Value');
            grid on;
        end

        function plot_gradient_norm(obj)
            if isempty(obj.gnorm)
                disp('No gradient norm data to plot.');
                return;
            end

            figure;
            yyaxis left;
            plot(obj.gnorm, 'LineWidth', 2); 
            ylabel('Linear Scale');
            xlabel('Iterations');

            yyaxis right;
            semilogy(obj.gnorm, 'LineWidth', 2); 
            ylabel('Log Scale');

            title('Gradient Norm');
            grid on;
        end

        function reset(obj)
            obj.fvals = [];
            obj.gnorm = [];
            obj.step_lengths = [];
            obj.timings = [];
        end
    end
end