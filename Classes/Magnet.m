classdef Magnet < handle
    %% 
    % Es wird eine Klasse definiert, da so sichergestellt werden kann, dass
    % alle Magneten gleiche Eigenschaften besitzen. Das mag zwar erst
    % umstaendlich sein, hilft aber bei der weiteren Programmierung
    properties
        b
        h
        x0
        y0
        mat_name
        mag_dir
        group % use group = 1 for left magnets, group = 2 for right ones
        type = 'normal';
    end
   
    methods
       
        function obj = Magnet(b, h, x0, y0, mat_name, mag_dir, group, type)
           obj.b = b;
           obj.h = h;
           obj.x0 = x0;
           obj.y0 = y0;
           obj.mat_name = mat_name;
           obj.mag_dir = mag_dir;
           obj.group = group;
           
           if exist('type', 'var')
               obj.type = type;
           end
        end   
        
        function bh_magnets = make_halbach(obj, start_dir, end_dir, number) 
            num_magnet = 1;
            dir = linspace(start_dir, end_dir, number);
            len = obj.h*number;
            y_list = linspace(-len/2+obj.h/2, len/2-obj.h/2, number);
            bh_magnets = [];
            
            for angle = dir
                x = obj.x0;
                y = obj.y0 + y_list(num_magnet);
                bh_magnets =   [bh_magnets
                                Magnet(obj.b, obj.h, x, y, obj.mat_name, angle, obj.group, 'halbach')];
                num_magnet = num_magnet+1;
            end            
        end

    end
   
end
