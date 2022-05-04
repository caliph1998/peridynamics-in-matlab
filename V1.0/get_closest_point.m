function [closestPoint] = get_closest_point(given_x, given_y, coord)

    S = size(coord); %THE TOTAL NUMBER OF MATERIAL POINTS HAS TO BE UPDATED
    TotalNumMatPoint = S(1);
    %if (left_tip.inEllipse(given_x, given_y) || right_tip.inEllipse(given_x, given_y))
    minDistance = inf;   
    for i = 1: TotalNumMatPoint            
            distance = sqrt((given_x - coord(i,1))^2 + (given_y - coord(i,2))^2);
            if (distance <= minDistance)
                closestPoint = i;
                minDistance = distance;
            end
        end
    %end
    
    