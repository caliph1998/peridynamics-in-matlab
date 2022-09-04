classdef EllipseClass
    properties
        x_center
        y_center
        radius_major
        radius_minor
    end
    methods
        function obj = ellipseClass(x, y, R, r)
            obj.x_center = x;
            obj.y_center = y;
            obj.radius_major = R;
            obj.radius_minor = r;
        end
        
        function bool = inEllipse(obj, x, y)
            bool = (((x - obj.x_center) / obj.radius_major) ^ 2 + ((y - obj.y_center) / obj.radius_minor) ^ 2) <= 1;
        end
    end
end