classdef Simulation_FEMM < handle
    properties
        type % planar/axi
        depth = 0;
        rad = 0;
        mag_list
    end
    
 
    methods
        function obj = Simulation_FEMM(type, param, mag_list)  
            if not(strcmp(type, 'planar') || strcmp(type, 'axi')) 
                error('Type must be planar or axi');
            end
            
            % Constructor
            obj.type = type; 
            
            if strcmp(type, 'axi') 
                obj.rad = param;
            else
                obj.depth = param;
            end
            
            obj.mag_list = mag_list;
        end
        
        
        function start_femm(obj, num, varargin)
           
            check_hidden = @(x) (x==0) || (x==1);
            p = inputParser;
            addOptional(p, 'hidden', 0, check_hidden);
            parse(p, varargin{:});
            hidden = p.Results.hidden;
            
            addpath('C:\femm42\mfiles')
            
            % open document and save as 'temp.fem'
            openfemm(hidden)
            newdocument(0)
            mi_saveas(strjoin({'C:\Temp\temp' num2str(num) '.fem'},''));

            mi_probdef(0, 'millimeters', obj.type, 1e-8, obj.depth, 18);
            mi_smartmesh(1); %
        end
        
        
        function get_materials(~) 
            % Loading materials directly from the *.txt-files did not work
            % properly

            % read materials from catalogue, for this you have to use the
            % matlib.dat-file from the Materials-folder or add the *.txt from
            % the Materials manually. Of course these materials are just
            % examples. You can use your own data. Try
            % https://apps.automeris.io/wpd/ to extract data from manufacturer
            % bh-curves
            
            mi_getmaterial('Air');
            mi_getmaterial('Pure Iron');
            mi_getmaterial('N38');
            mi_getmaterial('N38H_060C');
            mi_getmaterial('N38H_080C');
            mi_getmaterial('N38H_120C');
            mi_getmaterial('N48H_020C');
            mi_getmaterial('N48H_060C');
            mi_getmaterial('N48H_080C');
            mi_getmaterial('N48H_100C');
            mi_getmaterial('N48H_120C');
            mi_getmaterial('N35SH_080C');
            mi_getmaterial('N35SH_120C');
            mi_getmaterial('N38SH_080C');
            mi_getmaterial('N38SH_120C');
            mi_getmaterial('N48SH_080C');
            mi_getmaterial('N48SH_120C');
            mi_getmaterial('N40SH_080C');
            mi_getmaterial('N40SH_120C');
            mi_getmaterial('N45SH_020C');
            mi_getmaterial('N45SH_060C');
            mi_getmaterial('N45SH_080C');
            mi_getmaterial('N45SH_100C');
            mi_getmaterial('N45SH_120C');
        end
        
        
        function [size_magnets, acc] = set_environment(obj)
            % automatically scale the environment to the total size of the
            % magnets
            
            if obj.mag_list(1).h/2+abs(obj.mag_list(1).y0) > obj.mag_list(1).b
                size_magnets = obj.mag_list(1).h/2+abs(obj.mag_list(1).y0);
            else
                size_magnets = obj.mag_list(1).b;
            end
            
            if obj.mag_list(2).b > obj.mag_list(2).h
                acc_size = obj.mag_list(2).b;
            else
                acc_size = obj.mag_list(2).h;
            end
            
            
            size = size_magnets*3;
            acc = acc_size*0.05;
            
            % a round Boundary is drawn. There could not be noticed any
            % difference in an farfield (A=0) and an open BC as long as the
            % distance to the bc is about 5* the size of the model. 
            
            mi_drawarc(obj.rad, -size, obj.rad, size, 180, 1);
            mi_drawarc(obj.rad, size, obj.rad, -size, 180, 1);
            R = size;
            boundary = 0; % 2/(4*pi*10e-7*(R/1000));
            mi_addboundprop('A0', 0, 0, 0, 0, boundary, 0, 2, 0, 0); % R in m angeben!
            
            mi_selectcircle(obj.rad, 0, R, 4); 
            mi_setsegmentprop('A0', 0, 0, 0, 0);
            mi_clearselected

            mi_addblocklabel(obj.rad+1.5*size_magnets, 0);      
            mi_selectlabel(obj.rad+1.5*size_magnets, 0);       
            mi_setblockprop('Air', 0, acc, 0, 0, 0, 0);    

        end
        
        
        function set_refinement_air_gap(obj, r, size_magnets, acc) % optional
            mi_drawrectangle(obj.rad-r/2, -1.5*size_magnets, obj.rad+r/2, 1.5*size_magnets);

            mi_addboundprop('A0',0,0,0,0,0,0,0,0,0);
            mi_selectrectangle(0,0,0,0,2);
            mi_setsegmentprop('A0',0,0,0,0);
            mi_clearselected

            mi_addblocklabel(obj.rad,0);
            mi_selectlabel(obj.rad,0);
            mi_setblockprop('Air', 0, 0.4*acc,0,0,0,0);
        end
        
        
        function draw_magnet(obj, x, y, mag_num)
            
            % drawing magnets and grouping
            magnet = obj.mag_list(mag_num);
            x1 = obj.rad-magnet.b/2 + magnet.x0 + x;
            x2 = obj.rad+magnet.b/2 + magnet.x0 + x;
            y1 = -magnet.h/2 + magnet.y0 + y;
            y2 = magnet.h/2 + magnet.y0 + y;
            
            % automatic refinement:
            if (obj.mag_list(2).h > obj.mag_list(2).b)
                size_magnets = obj.mag_list(2).h;
            else
                size_magnets = obj.mag_list(2).b;
            end
            
            % for Halbach arrays the simulation needed to be much finer. If you
            % have a large diameter of the bearing you can also just use a
            % planar simulation with the diameter as the depth to reach better
            % convergence. 
            
            if strcmp(magnet.type, 'halbach') && strcmp(obj.type, 'axi')
                acc = 0.1;
            else
                acc = size_magnets*0.05;
            end
            
            % draw geometry
            mi_drawrectangle(x1, y1, x2, y2);

            % assign materials
            mi_addblocklabel(x1+(x2-x1)/2, y1+(y2-y1)/2);      
            mi_selectlabel(x1+(x2-x1)/2, y1+(y2-y1)/2);    
            mi_setblockprop(magnet.mat_name, 0, acc, 0, magnet.mag_dir, magnet.group, 0);            % mag_dir = Magnetisierungsrichtung

            % assign group
            mi_selectrectangle(x1, y1, x2, y2, 4);
            mi_setgroup(magnet.group);
        end
          
        
        function force = start_sim(~)
            % tell FEMM to start simulation
            mi_analyze();
            mi_loadsolution();

            mo_groupselectblock(1)
            % mo_groupselectblock(2)
            f_x = mo_blockintegral(18);
            f_y = mo_blockintegral(19);
            force = [f_x, f_y];
        end
        
        
        % calculate total area
        function Area = get_Area(obj)       
            Area = 0;
            for a = 1:length(obj.mag_list)
                    magnet = obj.mag_list(a);
                    if not(strcmp(magnet.mat_name, 'Pure Iron'))
                        Area = Area + magnet.b*magnet.h;
                    end
            end
        end
        
        
        % calculate total volume
        function Volume = get_Volume(obj)
            Volume = 0;
            if strcmp(obj.type, 'axi')
                for i = 1:length(obj.mag_list)
                    magnet = obj.mag_list(i);
                    if not(strcmp(magnet.mat_name, 'Pure Iron'))
                        if magnet.group == 1
                            r_i = obj.rad-magnet.b-4.5/2;
                            r_a = obj.rad-4.5/2; % gilt nur bei Luftspalt r=4.5mm
                        elseif magnet.group == 2
                            r_i = obj.rad+4.5/2;
                            r_a = obj.rad+magnet.b+4.5/2; % gilt nur bei Luftspalt r=4.5mm
                        end
                        Volume = Volume + pi*magnet.h*(r_a^2-r_i^2);
                    end
                end
            else
                for i = 1:length(obj.mag_list)
                    magnet = obj.mag_list(i);
                    if not(strcmp(magnet.mat_name, 'Pure Iron'))
                        Volume = Volume + magnet.h*magnet.b*obj.depth;
                    end
                end
            end
        end
        
        
        function b = get_b_on_line(obj, x1, y1, x2, y2, num_points)
            
            % defining line
            x = linspace(obj.rad+x1, obj.rad+x2, num_points);
            y = linspace(y1, y2, num_points);
            pt = [x', y'];
            b_list = ones(size(pt,1),2);
            
            % extract data on every point on line
            for i = 1:size(pt,1)
               b_list(i,:) = mo_getb(pt(i,1), pt(i,2));
            end
            b = [pt b_list];
        end
        
        
        function result = do_calculations(obj, z_ind, z, r, verbose, hidden, refinement)
            % magnets should all be the same size so refinment can be based on
            % the first one:
            
            b = obj.mag_list(1).b;
            
            obj.start_femm(z_ind, hidden); 
            

            obj.get_materials
            [size_magnets, acc] = obj.set_environment;       

            if refinement
                set_refinement_air_gap(obj, r, size_magnets, acc)
            end

            h = 0;
            % drawing magnets
            for i = 1:length(obj.mag_list)
                mag = obj.mag_list(i);
                h = h+mag.h/2;
                if mag.group == 1
                    obj.draw_magnet(-mag.b/2-r/2, z/2, i) 
                elseif mag.group == 2
                    obj.draw_magnet(mag.b/2+r/2, -z/2, i)
                end
            end
            force = obj.start_sim();                    
            if verbose
                % calculate flux density (beta)
                % b_1 under the bearing
                % b_2 right of the bearing
                % b_3 above the bearing

                b_1 = obj.get_b_on_line(-2*b, -h, 2*b, -h, 20);
                b_2 = obj.get_b_on_line(-2*b, h, -2*b, -h, 20);
                b_3 = obj.get_b_on_line(-2*b, h, 2*b, h, 20);

                result = [z, r, force, {b_1, b_2, b_3}];
            else
                result = [force(:,2)];
            end
            closefemm
        end
        
        
        % calculate forces as a function of sag (z direction) and air_gap (r direction)
        function result = make_study(obj, sag_list, air_gap_list, varargin)
            
            % input parsing
            validationFcn = @(x) islogical(x) && isscalar(x);
            check_hidden = @(x) (x == 0) || (x == 1);
            
            p = inputParser;
            addParameter(p, 'verbose', true, validationFcn);
            addParameter(p, 'save_results', false, validationFcn);
            addParameter(p, 'hidden', 0, check_hidden);
            addParameter(p, 'refinement', false);
            parse(p, varargin{:});
                        
            save_results = p.Results.save_results;   
            verbose = p.Results.verbose;
            hidden = p.Results.hidden;
            refinement = p.Results.refinement;
            
            result = [];
            % do calculatins for every combination of air_gap and sag
            res_index = 1;
            for air_gap = air_gap_list  
                if ~(length(sag_list) == 1)
                    parfor z_ind = 1:length(sag_list)
                        z = sag_list(z_ind);
                        res(z_ind,:) = obj.do_calculations(z_ind, z, air_gap, verbose, hidden, refinement)
                    end
                else 
                    z = sag_list;
                    res = obj.do_calculations(1, z, air_gap, verbose, hidden, refinement);
                end
                result = [result; res];
                res_index = res_index+1;
            end
            
            if save_results
                obj.save_as_mat(result)
            end
        end
         
        
        % calculating radial stiffness
        function radial_stiffness = get_k_r(obj, z, r, num_r)
            % because in axialsymmetric simulations the resulting radial force
            % always equals 0, the radial stiffness is calculated by two planar
            % simulations. One with a air_gap that gets bigger, one with an
            % air_gap getting smaller. The stiffness is calulated with a line
            % fit of the calculated curves. 
            %
            % Because in a round bearing, the force in x-direction does not mach
            % the force as in two straight rows of magnets, ther has to be a
            % correction factor. The shown correction factor only corrects for
            % the difference in direction. However the distance between
            % the two magnet rows is different from axial to planar what 
            % also has an effect on the forces. This effect was not accounted for. 
            
            if ~(isnumeric(r) && isscalar(r) && isnumeric(z) && isscalar(z))
                error('z und r muessen EIN skalarer Wert sein')
            end
            
            original_type = obj.type;
            original_depth = obj.depth;
            obj.type = 'planar';
            
            if strcmp(original_type, 'axi')
                obj.depth = pi*obj.rad; % calculations for half of the diameter
            else 
                obj.depth = original_depth/2; % again
            end    
                        
            r_list_left = linspace(1*r, 0.8*r, num_r);
            r_list_right = linspace(1*r, 1.2*r, num_r);
            res_left = obj.make_study(z, r_list_left, 'hidden', 1, 'refinement', true);
            res_right = obj.make_study(z, r_list_right, 'hidden', 1, 'refinement', true);

            f_r_planar_res = res_left(:,3)-res_right(:,3);
            f_r_planar_fit = fit(r_list_left', f_r_planar_res, 'poly1');
            k_r_planar = differentiate(f_r_planar_fit, r); 
            corr = 2/pi;   % because F_x = sin(phi)*F: corr = 1/(pi-0)*integral(sin(phi) d*phi) = 2/pi, intervall: (0 ... pi)
                                % F_r_res has a linear relationship with r
                                % so the stiffness can be multiplied by: corr = (2/pi)^2
  
            radial_stiffness = corr*k_r_planar;
            obj.type = original_type;
            obj.depth = original_depth;
        end
           
             
        function result = search_fmax(obj, r, z_min, z_max, num_pts, get_k, m_rotor, num_bearings)
            % worked faster than fminbnd() because the actual optimization is
            % done in matlab
            
            if isempty(get_k)
                get_k = false;
            end
            
            if isempty(m_rotor)
                m_rotor = 515;
            end
            
            if isempty(num_bearings)
                num_bearings = 2;
            end
            
            z_temp = linspace(z_min, z_max, num_pts);
            force_temp = abs(obj.make_study(z_temp, r, 'verbose', false, 'save_results', false, 'hidden', 1, 'refinement', true));
            f_fit = fit(z_temp', force_temp, 'poly5');    
            z_q = linspace(z_min, z_max, 1e5);
            f_max = max(f_fit(z_q));  
            z = z_q(f_fit(z_q) == f_max);
            
            
            if ~get_k
                result = [z, -f_max];                
            elseif get_k
                % searching the working point of the bearing (F=G_Rotor/number_of_bearings)
                g_rotor = 9.81*m_rotor/num_bearings;
                step = z_q(2)-z_q(1); 
                [~, idx] = min(abs(abs(f_fit(z_lower_max))-g_rotor));
                z_working_point = z_q(idx);
                k_r = obj.get_k_r(z_working_point, r, 2);
                k_z = differentiate(f_fit, z_working_point);
                result = [z, -f_max, z_working_point, k_z, k_r];
            end
        end
        
        
        function obj = scale_magnets(obj, scaling_factor_r, scaling_factor_s, varargin)
            % input parsing
            check_scaling = @(x) isnumeric(x);

            p = inputParser;
                        
            addParameter(p, 'scaling_b_r', scaling_factor_r, check_scaling);
            addParameter(p, 'scaling_h_r', scaling_factor_r, check_scaling);
            addParameter(p, 'scaling_x0_r', scaling_factor_r, check_scaling);
            addParameter(p, 'scaling_y0_r', scaling_factor_r, check_scaling);
            addParameter(p, 'scaling_b_s', scaling_factor_s, check_scaling);
            addParameter(p, 'scaling_h_s', scaling_factor_s, check_scaling);
            addParameter(p, 'scaling_x0_s', scaling_factor_s, check_scaling);
            addParameter(p, 'scaling_y0_s', scaling_factor_s, check_scaling);

            parse(p, varargin{:});
            
            scaling_b_r = p.Results.scaling_b_r;
            scaling_h_r = p.Results.scaling_h_r;
            scaling_x0_r = p.Results.scaling_x0_r;
            scaling_y0_r = p.Results.scaling_y0_r;
            scaling_b_s = p.Results.scaling_b_s;
            scaling_h_s = p.Results.scaling_h_s;
            scaling_x0_s = p.Results.scaling_x0_s;
            scaling_y0_s = p.Results.scaling_y0_s;
                    
            for mag_num = 1:length(obj.mag_list)
                magnet = obj.mag_list(mag_num);
                
                if magnet.group == 1
                    obj.mag_list(mag_num).b = magnet.b*scaling_b_s; 
                    obj.mag_list(mag_num).h = magnet.h*scaling_h_s;
                    obj.mag_list(mag_num).x0 = magnet.x0*scaling_x0_s;
                    obj.mag_list(mag_num).y0 = magnet.y0*scaling_y0_s;
                    
                elseif magnet.group == 2
                    obj.mag_list(mag_num).b = magnet.b*scaling_b_r;
                    obj.mag_list(mag_num).h = magnet.h*scaling_h_r;
                    obj.mag_list(mag_num).x0 = magnet.x0*scaling_x0_r;
                    obj.mag_list(mag_num).y0 = magnet.y0*scaling_y0_r;
                end
            end
        end
    end
end

