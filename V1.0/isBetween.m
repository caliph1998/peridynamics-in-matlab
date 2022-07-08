function [isBetween] = isBetween(coord, s, f, q)
    qs = coord(q,:) - coord(s,:);
    fs = coord(f,:) - coord(s,:);
    qs(1,3) = 0;
    fs(1,3) = 0;
    
    crossProduct = cross(qs, fs);
    if (norm(crossProduct) ~= 0)
        isBetween = false;
        return;
    end
    dotProduct = dot(qs, fs);
    if (dotProduct < 0)
        isBetween = false;
        return;
    end
    squareLength = norm(fs) ^ 2;
    if (dotProduct > squareLength)
        isBetween = false;
        return;
    end
    isBetween = true;
end