function subTours = detectSubtours(prevSolution,edges)
% Returns a cell array of subtours. The first subtour is the first row of x, etc.

%   Copyright 2014 The MathWorks, Inc. 

%This algorithm essentially starts at which starts at node 1, and
%follows the connected edges in one direction until it arrives back at the
%starting node, then add 1 to the number of subtours. Repeat as along as
%there are unvisited nodes. If the answer > 1, there are subtours.

% Commented by YiZhuang Garrard

prevSolution = round(prevSolution); % correct for not-exactly integers
prevSolutionEdgeIndices = find(prevSolution); % indices of the trips that exist in the solution
prevSolutionEdges = edges(prevSolutionEdgeIndices,:); % the collection of edges in the solution
unvisited = ones(length(prevSolutionEdgeIndices),1); % keep track of edges(?) not yet visited
curr = 1; % subtour we are evaluating
startour = find(unvisited,1); % first unvisited edge
    while ~isempty(startour)%While there are unvisted edges
        home = prevSolutionEdges(startour,1); % starting point of subtour
        nextpt = prevSolutionEdges(startour,2); % next point of tour
        visited = nextpt; unvisited(startour) = 0; % update unvisited points
        while nextpt ~= home
            % Find the other trips that starts at nextpt
            [srow,scol] = find(prevSolutionEdges == nextpt);
            % Find just the new trip
            trow = srow(srow ~= startour);
            scol = 3-scol(trow == srow); % turn 1 into 2 and 2 into 1
            startour = trow; % the new place on the subtour
            nextpt = prevSolutionEdges(startour,scol); % the point not where we came from
            visited = [visited,nextpt]; % update nodes on the subtour
            unvisited(startour) = 0; % update unvisited
        end
        subTours{curr} = visited; % store in cell array
        curr = curr + 1; % next subtour
        startour = find(unvisited,1); % first unvisited trip
    end
end