
[pg_centered, xy_centre, xy_scale] = zscore( [xcoords, ycoords] );
[ptp_centered, ptp_centre, ptp_scale] = zscore( ptp, 0, 'all' );


tic
% Create optimization variables
theta_hat = optimvar("theta_hat",1,4,"LowerBound",-10,"UpperBound",10);

% Set initial starting point for the solver
initialPoint.theta_hat = theta;

% Create problem
problem = optimproblem;

% Define problem objective
problem.Objective = fcn2optimexpr(@objectiveFcn,theta_hat,ptp_centered,...
    pg_centered);

% Display problem information
show(problem);

% Solve problem
[solution,objectiveValue,reasonSolverStopped] = solve(problem,initialPoint);

% Display results
disp(solution)
disp(reasonSolverStopped)
disp(objectiveValue)
toc

function objective = objectiveFcn(theta_hat, ptp_centered, pg_centered)
% This function should return a scalar representing an optimization objective.

% Example: Concession stand profit
% revenue = 3*soda + 5*popcorn + 2*candy;
% cost = 1*soda + 2*popcorn + 0.75*candy;
% objective = revenue - cost; % profit

% Edit the lines below with your calculations.
for c = 1:size( ptp_centered, 2 )
    objective = sum( ( ptp_centered(:,c) - ...
        ( theta_hat(1)./ sqrt( sum( (pg_centered - theta_hat([2,3])).^2, 2 ) + ...
        theta_hat(4).^2 ) ) ).^2, 'all' );
end
end