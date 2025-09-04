function [x,y] = sort_by_x(x,y)
    [x,sortID] = sort(x,'ascend');
    y = y(sortID);
end