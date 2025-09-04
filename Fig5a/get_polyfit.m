function [miny, minx, polycoeff] = get_polyfit(x_near_min, y_near_min)
    polycoeff = polyfit(x_near_min, y_near_min, 2);
    minx = -polycoeff(2)/(2*polycoeff(1));
    miny = polyval(polycoeff, minx);
end