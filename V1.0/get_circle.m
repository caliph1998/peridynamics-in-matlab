function [coord_excess] = get_circle(X, Y, R, dx, midCircle)
    
    x = X - R;
    y = Y + (R ^ 2 - (x - X) ^ 2) ^ 0.5;
    total = 0;
    while (x < X + R)
        segment_line = 2 * y;
        num_div_y = segment_line / dx;
        total = total + floor(num_div_y);
        x = x + dx;
        y = Y + (R ^ 2 - (x - X) ^ 2) ^ 0.5;
    end
    %%
    coord_excess = zeros(total, 2);
    x = X - R;
    y = Y + (R ^ 2 - (x - X) ^ 2) ^ 0.5;
    nnum = 0;

    while (x < X + R)
        segment_line = 2 * y;
        num_div_y = segment_line / dx;
        for j = 1: floor(num_div_y)
            coordx = x;
            coordy = (-1) * y + (j - 1) * dx;
             if (midCircle.inEllipse(coordx, coordy))
                continue;
             end
            
            nnum = nnum + 1;
            coord_excess(nnum, 1) = coordx;
            coord_excess(nnum, 2) = coordy;

        end
        x = x + dx;
        y = Y + (R ^ 2 - (x - X) ^ 2) ^ 0.5; 

    end
    coord_excess = coord_excess(1:nnum, :);
end

