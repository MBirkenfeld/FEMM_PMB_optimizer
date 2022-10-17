classdef Optimization_FEMM < handle
    properties
        % parameters for optimization
        sim_handle 
        m_rotor
        num_bearings 
        f_goal
        air_gap
        R % radius at the middle of the air gap
        num_points_per_eval
        lower_boundaries
        upper_boundaries
        starting_point
        mat_Stator 
        mat_Rotor
        
        % Cache (as far as I know this is the only way to do this):
        sim
        Area
        f_z_max
        f_max_values
        sf_f_max    
    end
    
    
    methods
        function obj = Optimization_FEMM(sim_handle, lower_boundaries, upper_boundaries, starting_point, varargin)
            obj.sim_handle = sim_handle;
            obj.lower_boundaries = lower_boundaries;
            obj.upper_boundaries = upper_boundaries;
            obj.starting_point = starting_point;
            
            check_if_num  = @(x) isnumeric(x) && isscalar(x);
            
            p = inputParser;
            addOptional(p, 'm_rotor', 1000, check_if_num);
            addOptional(p, 'num_bearings', 1, check_if_num);
            addOptional(p, 'air_gap', 5, check_if_num);
            addOptional(p, 'R', 500, check_if_num);
            addOptional(p, 'num_points_per_eval', 12, check_if_num);
            addOptional(p, 'mat_Stator', 'N38')
            addOptional(p, 'mat_Rotor', 'N38')
            
            parse(p, varargin{:});
            
            obj.m_rotor = p.Results.m_rotor;
            obj.num_bearings = p.Results.num_bearings; 
            obj.air_gap = p.Results.air_gap;
            obj.R = p.Results.R;
            obj.num_points_per_eval = p.Results.num_points_per_eval;
            obj.mat_Stator = p.Results.mat_Stator;
            obj.mat_Rotor = p.Results.mat_Rotor;
            
            obj.f_goal = 1.5*obj.m_rotor*9.81/obj.num_bearings; 
        end
        
        
        function E = scale2fit_f_max(obj, x)
            sf = x;
            obj.sim.scale_magnets(sf, sf);
            
           % calculate f_max
            lower_bound = 0.2*obj.sim.mag_list(2).h;
            upper_bound = 1.4*obj.sim.mag_list(2).h;

            obj.f_max_values = obj.sim.search_fmax(obj.air_gap, lower_bound, upper_bound, obj.num_points_per_eval, false, obj.m_rotor, obj.num_bearings);
            obj.f_z_max = obj.f_max_values(2);    

            obj.Area = obj.sim.get_Area;
            E = abs(1-(abs(obj.f_z_max)/obj.f_goal));
            obj.sim.scale_magnets(1/sf, 1/sf);
        end
        
        
        function Area =  calc_variant(obj, x)    
            handle = obj.sim_handle;            
            
            geometry = [10, x*10];

            obj.sim = handle(obj.R, geometry, obj.mat_Stator, obj.mat_Rotor);

            options_scale2fit = optimset('MaxIter', 15, 'TolFun', 5e-3);    
            [sf_f_max, Error] = fminsearch(@obj.scale2fit_f_max, 1, options_scale2fit); % minimizes deviation of f_goal/f_max
            obj.sf_f_max = sf_f_max;
            Area = obj.Area;
        end
        
        
        function results = optimize(obj, varargin)
            
            p = inputParser;
            addOptional(p, 'save_results', false);
            addOptional(p, 'file_name', 'temp');
            parse(p, varargin{:});
            
            save_results = p.Results.save_results;
            
            A = [];
            b = [];
            Aeq = [];
            beq = [];       
            lb = obj.lower_boundaries;
            ub = obj.upper_boundaries;
            nonlcon = [];
            options = optimset('PlotFcns',{'psplotbestf', 'psplotbestx'}, 'TolFun', 8e-3, 'MaxIter', 40, 'Display', 'iter');
            
            optimization_fun = @(x) obj.calc_variant(x);
            [x, fval, exitflag, output] = patternsearch(optimization_fun, obj.starting_point ,A,b,Aeq,beq,lb,ub,nonlcon, options); % ,A,b,Aeq,beq,lb,ub,nonlcon,

            results.x = x;
            results.Area = fval;
            results.exitflag = exitflag;
            results.ouput = output;
            results.optimization = struct(obj);
            results.sim = struct(obj.sim);
            results.sf_f_max = obj.sf_f_max;
            
            if save_results
                save({'results\', file_name, '.mat'}, 'result')
            end
        end 
    end
end
