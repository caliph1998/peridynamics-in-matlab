%%%The function below calculates the stresses in the x and y direction of
%%%all the given material points and returns a stress array of size = [totalNumMatPoint, 2].
function [stress] = CalculateStressforPoint(coord,TotalNumMatPoint,numfam,nodefam, pointfam, thick)
    stress = zeros(TotalNumMatPoint, 2); %Predifining the size of the stress array to speed up the calculations.

    for i = 1: (TotalNumMatPoint)
        f_x = 0; %Sum of all forces in X direction
        f_y = 0; %Sum of all forces in Y direction
        satisfy_x = 0;
        satisfy_y = 0;
        point_x = coord(i,1); %% Extract coordinates of point i
        point_y = coord(i,2);
        for j = 1: TotalNumMatPoint %% Iterating every other material point for conditions
            if (coord(j,2) == point_y && coord(j,1) < point_x + thick / 2) %% coord(j,2) == point_y &&
                    satisfy_x = satisfy_x + 1;
                    for k = 1:numfam(j,1) %%Iterate for all the neighbors of point j that satisfy the conditions.
                        cnode = nodefam(pointfam(j,1)+k-1,1);
                        if (coord(cnode,1) > point_x + thick / 2) % Sum bondforce of the neighbors that satisfy this condition
                            f_x = f_x + nodefam(pointfam(j,1)+k-1, 2); %%% Requires Saving BondForce_x in nodefam(:,2)
                        end
                    end
            end
            %Same below for stress in y direction
            if (coord(j,1) == point_x && coord(j,2) < point_y + thick / 2) %%%%coord(j,1) == point_x &&
                    satisfy_y = satisfy_y + 1;
                    for k = 1:numfam(j,1) %%Iterate for all the neighbors of point j that satisfy the conditions.
                        cnode = nodefam(pointfam(j,1)+k-1,1);
                        if (coord(cnode,2) > point_y + thick / 2) % Sum bondforce of the neighbors that satisfy this condition
                            f_y = f_y + nodefam(pointfam(j,1)+k-1, 3); %%% Requires Saving BondForce_x in nodefam(:,2)
                        end
                    end
            end
            stress_x = f_x * thick;
            stress_y = f_y * thick;
        end
        %Saving the stress for point i before moving on to the next point
        %on the plate
        stress(i,1) = stress_x;
        stress(i,2) = stress_y;
    end
    satisfy_x
    satisfy_y
end
            