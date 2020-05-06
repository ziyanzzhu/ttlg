classdef Layer
    % Basic datastructure to store geometric info for a 2D layer.
    
    properties
        layer_index     % index in the heterostructure
        A               % global lattice
        G               % global recip lattice
        A_local         % local lattice (theta = 0)
        G_local         % local recip lattice (theta = 0)
        theta           % twist/orientation angle
    end
    
    methods
        function obj = Layer(index,theta)
            %Layer Construct an instance of this class
            %   Detailed explanation goes here
            obj.layer_index = index;
            obj.theta = theta;
            
            alpha = 1.42*sqrt(3);
            
            obj.A_local = alpha*[   1       1/2
                                    0 sqrt(3)/2];
            obj.G_local = getRecip(obj.A_local);
            
            R_base = [cos(theta) -sin(theta);
                       sin(theta) cos(theta)];
            
            obj.A = R_base*obj.A_local;
            obj.G = getRecip(obj.A);
            
        end
    end
end

