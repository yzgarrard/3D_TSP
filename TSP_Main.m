clc;
clear;
close;

startTime = cputime;

%% Create a map with stops
figure;
rng('shuffle'); %Random number generator using current time as seed
%x = [-400, 400, 400, -400];    %x-coordinates of polygon
%y = [400, 400, -400, -400];    %y-coordinates of polygon
nStops = 150;    %Number of stops
stopsLon = zeros(nStops, 1);    %Initialize array for lons
stopsLat = stopsLon;            %Initialize array for lats
randomPoints = true;
n = 1;  %Counter for stop generation
while (n <= nStops)
    xp = 0;
    yp = 0;
    %if inpolygon(xp, yp, x, y)
    %if (n+2 > nStops)
    %break;
    %end
    if (~randomPoints)
        stopsLon(n) = 0;
        stopsLat(n) = (n) / nStops;
        stopsLon(n+1) = 1;
        stopsLat(n+1) = (n) / nStops;
        n = n + 2;
    end
    if (randomPoints)
        stopsLon(n) = rand;
        stopsLat(n) = rand;
        n = n + 1;
    end
    
    %end
end
%stopsLon(nStops) = -1;
%stopsLat(nStops) = -1;
%plot(x, y, 'Color', 'red');
hold on
plot(stopsLon, stopsLat, '*b')
hold off

%pause(1)

%% Create all possible edges between two stops
edges = nchoosek(1:nStops, 2);


%% Calculate distances between each point
distances = hypot(stopsLat(edges(:,1)) - stopsLat(edges(:,2)), ...
    stopsLon(edges(:,1)) - stopsLon(edges(:,2)));
lengthDistance = length(distances);  %How many elements in distances (number of edges)

%% Equality constraints

%This makes a sparse array. No idea why its done this way
%First type of equality constraint. States that the number of
%edges must equal the number of nodes
Aeq = spones(1:length(edges));  %Creates array from 1 to length(edges), then sets all nonzeros to 1
beq = nStops;

%Second type of equality constraint makes it so two edges are
%attached to each node.
%Essentially, if node ii is part of an edge, make it part of the constraint
%e.g. (1,2), (1,3), (1,4), (2,3), (2,4), (3,4) are all edges
%row 2 of Aeq would be 1, 1, 1, 0, 0, 0, for edges that contain node 1
%row 3 of Aeq would be 1, 0, 0, 1, 1, 0, for edges that contain node 2
%row 4 of Aeq would be 0, 1, 0, 1, 0, 1, for edges that contain node 3
%row 5 of Aeq would be 0, 0, 1, 0, 1, 1, for edges that contain node 4
%rows 2-5 of beq would be 2, to show that only two edges are allowed per
%node
Aeq = [Aeq; spalloc(nStops, length(edges), nStops * (nStops - 1))]; %Extends Aeq
for ii = 1:nStops
    whichEdges = (edges == ii); %Find the trips that include stop ii
    whichEdges = sparse(sum(whichEdges, 2));    %include trips where ii is at either end
    Aeq(ii + 1, :) = whichEdges';%Append to constraint matrix
end
beq = [beq; 2 * ones(nStops, 1)];

%% Binary bounds

%Make all the variables binary. Finally, something I don't have to analyze
%for half an hour!
intcon = 1:lengthDistance;  %Set all edge variables to integers
lb = zeros(lengthDistance, 1);  %Lower bound of 0
ub = ones(lengthDistance, 1);   %Upper bound of 1

%% Optimize using intlinprog

opts = optimoptions('intlinprog', 'Display', 'off', ...
    'CutMaxIterations', 50, 'CutGeneration', 'advanced');
[x_tsp, costopt, exitflag, output] = intlinprog(distances, intcon, [], [], Aeq, beq, lb, ub, opts);

%% Visualize the solution

segments = find(x_tsp); % Get indices of lines on optimal path
lh = zeros(nStops,1); % Use to store handles to lines on plot
lh = updateSalesmanPlot(lh,x_tsp,edges,stopsLon,stopsLat);
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
    lh = updateSalesmanPlot(lh, x_tsp, edges, stopsLon, stopsLat);
    
    %Recount subtours
    subTours = detectSubtours(x_tsp, edges);
    numTours = length(subTours);
    fprintf('# of subtours: %d\n', numTours);
end


fprintf("Shortest path between %d nodes. Elapsed time: %fs\n", nStops, cputime - startTime);
nStops = nStops + 2;
title('Solution with Subtours Eliminated');
hold off