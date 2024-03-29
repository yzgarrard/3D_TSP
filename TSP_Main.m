clear;
close;

startTime = cputime;

%% Create a map with stops
figureHandle = figure;
hold;
[X,Y] = meshgrid(-5*pi:1:105*pi,-5*pi:1:105*pi);
Z = 400*sin(X./100) + 400*(cos(Y./100-pi)+1);
surface(X,Y,Z,'EdgeColor','none','FaceAlpha',0.87);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Increase spacing to decrease number of nodes and solve problem faster
spacing = 25*pi;
[U,V] = meshgrid(0:spacing:100*pi,0:spacing:100*pi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xypoints = [U(:), V(:)];
stops_x = xypoints(:,1);
stops_y = xypoints(:,2);
stops_z = 400*sin(stops_x./100) + 400*(cos(stops_y./100-pi)+1)+50;
surfPoints = [stops_x,stops_y,stops_z];
scatter3(stops_x, stops_y, stops_z);
k=1:length(stops_x); 
text(stops_x,stops_y,stops_z,num2str(k'))
view(3)


%% Create all possible edges between two stops
edges = nchoosek(1:length(surfPoints), 2);

%% Calculate distances between each point
distances = sqrt(  (stops_x(edges(:,1)) - stops_x(edges(:,2))).^2 ...
                  +(stops_y(edges(:,1)) - stops_y(edges(:,2))).^2 ...
                  +(stops_z(edges(:,1)) - stops_z(edges(:,2))).^2);
lengthDistance = length(distances);  %How many elements in distances (number of edges)
nStops = length(surfPoints);

%% Equality constraints

%These equality constraints makes it so two edges are
%attached to each node. The edges represents the drone going in and out
% of the node
%
%If node ii is part of an edge, make it part of the constraint
%e.g. (1,2), (1,3), (1,4), (2,3), (2,4), (3,4) are all edges
%row 1 of Aeq would be 1, 1, 1, 0, 0, 0, for edges that contain node 1
%row 2 of Aeq would be 1, 0, 0, 1, 1, 0, for edges that contain node 2
%row 3 of Aeq would be 0, 1, 0, 1, 0, 1, for edges that contain node 3
%row 4 of Aeq would be 0, 0, 1, 0, 1, 1, for edges that contain node 4
%rows 1-4 of beq would be 2, to show that only two edges are allowed per
%node
Aeq = [];
%[Aeq; spalloc(nStops, length(edges), nStops * (nStops - 1))]; %Extends Aeq
for ii = 1:nStops
    whichEdges = (edges == ii); %Find the trips that include stop ii
    whichEdges = sparse(sum(whichEdges, 2));    %include trips where ii is at either end
    Aeq(ii, :) = whichEdges';%Append to constraint matrix
end
beq = [2 * ones(nStops, 1)];

%% Binary bounds

%Make all the variables binary so that for each edge between two nodes, it
% either traversed or not traversed
intcon = 1:lengthDistance;  %Set all edge variables to integers
lb = zeros(lengthDistance, 1);  %Lower bound of 0
ub = ones(lengthDistance, 1);   %Upper bound of 1

%Force some edges to never be traveled because they would collide with 
% the ground
syms x y z
surfEq = 400*sin(x/100) + 400*(cos(y/100-pi)+1) - z;
for i = 1:size(edges,1)
    edge = edges(i,:);
    doesCollide = checkCollision(stops_x(edge(1)),stops_x(edge(2)),stops_y(edge(1)),stops_y(edge(2)),stops_z(edge(1)),stops_z(edge(2)),surfEq);
    if doesCollide == true
        ub(i) = 0;
        fprintf("Forbidding edge %d %d\n", edge(1), edge(2))
    end        
end

fprintf("Elapsed time: %fs\n",  cputime - startTime);

%% Optimize using intlinprog

opts = optimoptions('intlinprog', 'Display', 'off', ...
    'CutMaxIterations', 50, 'CutGeneration', 'advanced');
[x_tsp, costopt, exitflag, output] = intlinprog(distances, intcon, [], [], Aeq, beq, lb, ub, opts);

%% Visualize the solution

segments = find(x_tsp); % Get indices of lines on optimal path
lh = zeros(nStops,1); % Use to store handles to lines on plot
lh = update3DPlot(lh,x_tsp,edges,stops_x,stops_y,stops_z);
title('Solution with Subtours');

%% Detect subtours

subTours = detectSubtours(x_tsp, edges);
numTours = length(subTours); %Number of subtours
fprintf('# of subtours: %d\n', numTours);

%% Subtour constraints

A = spalloc(0, lengthDistance, 0);
b = [];
while numTours > 1
    %Add the subtour constraints
    A = [A; spalloc(numTours, lengthDistance, nStops)];
    b = [b; zeros(numTours, 1)];
    for ii = 1:numTours
        rowEdge = size(A, 1) + 1;
        subTourEdge = subTours{ii};
        subTourEdgeVariations = nchoosek(1:length(subTourEdge), 2);
        for jj = 1:length(subTourEdgeVariations)
            whichVariation = (sum(edges==subTourEdge(subTourEdgeVariations(jj,1)),2)) & ...
                (sum(edges==subTourEdge(subTourEdgeVariations(jj,2)),2));
            A(rowEdge, whichVariation) = 1;
        end
        b(rowEdge) = length(subTourEdge) - 1;
    end
    
    %Optimize again
    [x_tsp, costopt, exitflag, output] = intlinprog(distances, intcon, A, b, Aeq, beq, lb, ub, opts);
    
    %Update plot
    lh = update3DPlot(lh, x_tsp, edges,stops_x,stops_y,stops_z);
    
    %Recount subtours
    subTours = detectSubtours(x_tsp, edges);
    numTours = length(subTours);
    fprintf('# of subtours: %d\n', numTours);
end

D = sum(distances.*x_tsp);
fprintf("Shortest path between %d nodes. Elapsed time: %fs\n", nStops, cputime - startTime);
title('Solution with Subtours Eliminated');
fprintf("Total distance: %8.4f m\n", D)
hold off

%% Second optimization for time vs energy consumption
A = 1;  % Weight to prioritize time
B = 1;  % Weight to prioritize energy consumption
fun = @(V) (A*D/V)+(B*(3.1214*sqrt(96.236+0.0296*V^2)^(3/2)+0.0296*V^3));
x0 = 1;
A = -1;
b = 0;
lb = realmin;   % Don't want zero, that's infinite time.
ub = 20;    % 20m/s is a reasonable max speed. Not good for photography though

options = optimoptions('fmincon','Display','off');
V = fmincon(fun, x0, A, b, [], [], lb, ub, [], options);
fprintf("Velocity: %12.4f m/s\nEnergy Consumption: %12.4f Watts\nTime: %12.4f s\n",...
    V, (3.1214*sqrt(96.236+0.0296*V^2)^(3/2)+0.0296*V^3), D/V)
