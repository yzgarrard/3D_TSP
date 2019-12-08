clc;
clear;
close;

startTime = cputime;

%% Create a map with stops
figureHandle = figure;
hold;
[X,Y] = meshgrid(0:1:100*pi,0:1:100*pi);
Z = 400*sin(X./100) + 400*(cos(Y./100-pi)+1);
surface(X,Y,Z,'EdgeColor','none','FaceAlpha',0.5);
[A,B] = meshgrid(0:20*pi:100*pi,0:20*pi:100*pi);
xypoints = [A(:), B(:)];
A = xypoints(:,1);
B = xypoints(:,2);
C = 400*sin(A./100) + 400*(cos(B./100-pi)+1)+50;
surfPoints = [A,B,C];
scatter3(A, B, C);
k=1:length(A); 
text(A,B,C,num2str(k'))
view(3)

stops_x = A;
stops_y = B;
stops_z = C;

%% Create all possible edges between two stops
edges = nchoosek(1:length(surfPoints), 2);

%% Calculate distances between each point
distances = sqrt(  (A(edges(:,1)) - A(edges(:,2))).^2 ...
                  +(B(edges(:,1)) - B(edges(:,2))).^2 ...
                  +(C(edges(:,1)) - C(edges(:,2))).^2);
lengthDistance = length(distances);  %How many elements in distances (number of edges)

% In case I need an example of how to plot a line
% testX = [surfPoints(1,1), surfPoints(8,1)];
% testY = [surfPoints(1,2), surfPoints(8,2)];
% testZ = [surfPoints(1,3), surfPoints(8,3)];
% plot3(testX,testY,testZ);

nStops = length(surfPoints);

%% Equality constraints

%This makes a sparse array. No idea why its done this way
%First type of equality constraint. States that the number of
%edges must equal the number of nodes
%Aeq = spones(1:length(edges));  %Creates array from 1 to length(edges), then sets all nonzeros to 1
%beq = nStops;

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


fprintf("Shortest path between %d nodes. Elapsed time: %fs\n", nStops, cputime - startTime);
nStops = nStops + 2;
title('Solution with Subtours Eliminated');
hold off