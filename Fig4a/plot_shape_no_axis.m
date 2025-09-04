function f = plot_shape_no_axis(R, d, phi, kappa, Sigma, varargin)


% parse inputs and set defaults
args = varargin;
nargs = numel(args);
debug_on = false;
plot_colour = 'r';
k = 1;
while k<=nargs
    % if we get a figure handle, use it
    if strcmpi(args{k},'Debug_on')||strcmpi(args{k},'DebugOn')
        debug_on = args{k+1};
        k = k+1;
    elseif strcmpi(args{k},'colour')||strcmpi(args{k},'color')
        plot_colour = args{k+1};
        k = k+1;
    else
        error('intput not recognised!');
    end
    k = k+1;
end

% just call nonlinear
out = free_shape_nonlinear_free_h(R, d, phi, kappa, Sigma, debug_on);

r = out.y(2,:);
h = out.y(3,:);

figure('Position',[400,100,800,600]);
hold on
axis equal
xlabel('$r$')
ylabel('$h$')
plot(r, h-h(1), '-', 'displayname', 'free surface', 'Color',plot_colour, 'linewidth', 5);
plot(-r, h-h(1), '-', 'displayname', 'free surface', 'Color',plot_colour, 'linewidth', 5);
t = linspace(-pi,pi,1000);
x = cos(t)*R;
% y = sin(t)*R+(R*cos(phi)+h(1));
y = sin(t)*R+R*cos(phi);
plot(x,y,':', 'displayname', 'microbead', 'Color',plot_colour, 'linewidth', 5)
t = linspace(-pi/2,-pi/2+phi,1000);
x = cos(t)*R;
% y = sin(t)*R+(R*cos(phi)+h(1));
y = sin(t)*R+R*cos(phi);
plot(x,y,'-', 'displayname', 'microbead', 'Color',plot_colour, 'linewidth', 5)
plot(-x,y,'-', 'displayname', 'microbead', 'Color',plot_colour, 'linewidth', 5)
xlim([-3*R,3*R])
axis off


end