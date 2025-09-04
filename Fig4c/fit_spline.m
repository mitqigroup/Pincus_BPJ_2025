function yy = fit_spline(x,y,xx, remove_points, varargin)

% parse inputs and set defaults
args = varargin;
nargs = numel(args);
weights_range = 3;
weights_val = 1e3;
use_weights = false;
weights = ones(size(y));
k = 1;
while k<=nargs
    % if we get a figure handle, use it
    if strcmpi(args{k},'WeightsRange')||strcmpi(args{k},'Weights_Range')
        weights_range = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'UseWeights')||strcmpi(args{k},'Use_Weights')
        use_weights = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'WeightsEndValues')||strcmpi(args{k},'Weights_End_Values')
        weights_val = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'Weights')
        weights = args{k+1};
        k = k+1;
    else
        error('intput not recognised!');
    end
    k = k+1;
end

assert(~isempty(x));
assert(~isempty(y));
assert(~all(isnan(x)));
assert(~all(isnan(y)));

y(remove_points) = nan;
if use_weights
    weights(1:1+weights_range) = weights_val;
    weights(end-weights_range:end) = weights_val;
end
[pp,~] = csaps(x, y, 1-1e-5, [], weights);
yy = fnval(pp,xx);