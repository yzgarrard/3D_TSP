function lh = update3DPlot(lh,xopt,idxs,stops_x,stops_y,stops_z)
% Plotting function for tsp_intlinprog example

%   Copyright 2014-2016 The MathWorks, Inc. 

% Modified by YiZhuang Garrard for 3D plotting

if ( lh ~= zeros(size(lh)) ) % First time through lh is all zeros
    set(lh,'Visible','off'); % Remove previous lines from plot
end

segments = find(round(xopt)); % Indices to trips in solution

% Loop through the trips then draw them
x = zeros(3*length(segments),1);
y = zeros(3*length(segments),1);
z = zeros(3*length(segments),1);
for ii = 1:length(segments)
    start = idxs(segments(ii),1);
    stop = idxs(segments(ii),2);
    
    % Separate data points with NaN's to plot separate line segments
    x(3*ii-2:3*ii) = [stops_x(start); stops_x(stop); NaN];
    y(3*ii-2:3*ii) = [stops_y(start); stops_y(stop); NaN]; 
    z(3*ii-2:3*ii) = [stops_z(start); stops_z(stop); NaN];
end

lh = plot3(x,y,z,'r-','LineWidth',2);

set(lh,'Visible','on'); drawnow; % Add new lines to plot