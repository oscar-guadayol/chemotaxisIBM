function [bacterialDistribution, rBins, nCells, nutrientExposure, simTime] =...
    chemotaxisIBM_phycosphere(...
    nCells, timeStep, intTimes, rDomain, gradientFunction,...
    meanTumbleTime, meanRunLength, rotDiffCoeff,...
    tumbleAnglesDistribution, meanRunSpeed)


% chemotaxisIBM_phycosphere simulates a population of Escherichia coli
% chemotactic cells swimming in a 2D circular domain with a gradient of the
% chemoattractant aspartate around a spherical phytoplankton cell of given
% diameter.
% 
% The chemotactic behaviour is modelled using the chemosensory circuit of
% E. coli for aspartate, following Tu et al 2008 and Kalinin et al 2009.
%
% The motility pattern is defined by the rotational diffusivity
% coefficient of the cells and by their tumble time and angle, running
% speed and run length. The model allows for bimodal motility patterns,
% where the cells alternate between different reorientation behaviours.
%
% Both boundaries (the internal representing the phytoplankton cell wall)
% boundary are reflective: cells bounce back from them at the same angle of
% incidence with respect to the tangent of the boundary at the point of
% intersection.
%
% Input:
%   nCells is the number of E. coli cells in the simulation.
%
%   timeStep is the length of the simulation time steps in seconds.
%
%   intTimes are the interrogation times for the simulation in seconds. At
%       each interrogation time, the bacterial distribution along the
%       gradient is assessed.
%
%   rDomain is a 2 elements vector with the minimum and mximum radii of the
%       simulation domain in microns. The minimum radius represents the
%       radius of the spherical phytoplankton cell.
%
%   gradientFunction is a function handle to compute the concentration of
%       aspartate in microM as a function of radial position.
%       For example:
%       gradientFunction = @(r) Csurface*CellRadius./r;
%       where Csurface is the concentration of aspartate at the surface of
%       the phytoplankton cell, Cell radius is the radius of the cell and r
%       is the distance to the center of the cell in microns,
%
%   meanTumbleTime is the mean tumble time in s.
%
%   meanRunLength is the mean run length in s.
%
%   rotDiffCoeff is the rotational diffusion coefficient along the long
%       semiaxis in radians^2/s
%
%   tumbleAnglesDistribution is a cell array in which each cell contains a
%       vector with a population of tumble angles from which values are
%       drawn randomly in each tumble event. The number of cells marks the
%       number of reorientation modes the cells display. Each cells changes
%       the distribution of angles the on every successive reorientation.
%       This allows to model, for example, run and reverse or
%       run-flick-reverse motility patterns.
%
%   meanRunSpeed is the mean speed during runs in microns/s.
%
%
% Output:
%   bacterialDistribution is a NxM matrix with the concentration of
%       bacteria per distance to source relative to background
%       concentrations. N is the number of 1 micron bins in the rDomain, M
%       is the number of intTimes.
%
%   rBins are the bin edges, in microns, used to partition bacterial
%       distributions.
%
%   simTimes is a vector of simulation times, in seconds.
%
%   nCells is a vector with the number of cells within the simulation
%       domain at each time in simTimes.
%
%   nutrientExposure is a vector with the average nutrient exposure in
%       microM/cell/s at each time in simTimes.
%
% Dependencies: Statistics toolbox, Statistics and Machine Learning Toolbox
%
% Bibliography:
%
%     Kalinin YV, Jiang L, Tu Y, Wu M. Logarithmic Sensing in Escherichia
%         coli Bacterial Chemotaxis. 2009;96(6):2439–48. 
%     Shimizu TS, Tu Y, Berg HC. A modular gradient‐sensing network for
%         chemotaxis in Escherichia coli revealed by responses to
%         time‐varying stimuli. Molecular Systems Biology. 2010 Jan
%         1;6(1):382.
%     Tu Y, Shimizu TS, Berg HC. Modeling the chemotactic response of 
%         Escherichia coli to time-varying stimuli.
%         PNAS. 2008 Sep 30;105(39):14855–60. 
%     Vladimirov N, Løvdok L, Lebiedz D, Sourjik V. Dependence of Bacterial 
%         Chemotaxis on Gradient Shape and Adaptation Rate.
%         PLOS Computational Biology. 2008 Dec 19;4(12):e1000242. 
% 
% Copyright (C) 2023,  Oscar Guadayol oscar_at_guadayol.cat
%
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License, version 3.0, as
% published by the Free Software Foundation.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.

%% Spatial domain
rMax = max(rDomain); % Radius of outer domain.
rMin = min(rDomain); % Radius of phytoplankton cell.
dr = 1; % width of the bins.
rBins = rMin:dr:rMax; % bins used to evaluate bacterial distributions.


%% Initial bacterial population
% Initial x and y positions and orientations of bacteria assuming uniform
% distributions. Suffix "In" indicates cells that are inside the simulation
% domain. Orientations are the orientations in radians of the cells, rho is
% the radial distance in microns to the center of the phytoplanktonc cell,
% X and Y are the coordinates of the cells in the bidimensional space.
[Orientations, rho, X, Y] =...
    populateDomain(rMax, rMin, nCells(1));

% Abundance of cells:
BackgroundAbundance2D = nCells(1)/(pi*rMax^2-pi*rMin^2); % Cells m^-2


%% Chemotactic pathway constants
% From Kalinin et al 2009, Tu et al 2008, Shimizu et al 2010.

%%% Phosporilation constants %%%

% Number of responding receptor dimers in the MCP (Methyl-accepting
% Chemotactic Proteins):
nReceptors = 6;

% Dissociation constant of the ligand to the inactive receptor in microM:
K_I = 18;

% K_A = K_I/C dissociation constant of the ligand to the active receptor:
K_A = K_I/0.0062;

% Free-energy change per added methyl group (in units of kT). In Shimizu et
% al (2010) alpha = 2kT at ~22C.
ALPHA = 1.7;

% Methylation level at which the methylation level dependent
% freeEnergyDifference is 0.
M_0 = 1;


%%% Methylation constants %%%

% Methylation rate for the inactive receptors in s^-1:
k_R = 0.005;

% De-methylation rate for the active receptors in s^-1:
k_B = 0.005;


%%% Tumble probability constants %%%

% Hill coefficient:
H = 10.3;

% Steady-state activity in the absence of stimuli.
% In Shimizu et al (2010) a_0 = 1/3 at 32C.
A_0 = 1/2;

% Clockwise bias in the absence of a chemoattactant gradient:
cwBiasHomogeneous = meanTumbleTime/(meanTumbleTime+meanRunLength);

% Kinase activity that induces a clockwise bias of 0.5 in the absence of a
% chemoattactant gradient:
aHalf = A_0*(1/cwBiasHomogeneous-1)^(1/H);

meanTumbleTime = single(meanTumbleTime);

%%% Initial methylation levels %%%

% Calculates methylation levels M assuming all cells are at steady state
% for each location along the gradient in L and assigns the correct M level
% at each particle according to its position in rho.

methylationLevel =...
    equilibriumMethylationLevel(gradientFunction(rho));

%% Change random number generator to 'simdTwister',
% SIMD-oriented Fast Mersenne Twister:
rng(1,'simdTwister')
% rng('default') %reverts changes

%% Declare variables for IBM simulation

% Simulation time in seconds:
simTime = 0:timeStep:max(intTimes);

% Last iteration in the main IBM loop:
maxIteration = round(max(intTimes)/timeStep)+1;

% nCells-long vector with the iterations in which the last run started for
% each cell. If for a particular cell it is lower or equal than the current
% iteration value, then the cell is actively running. Otherwise the cell is
% tumbling. Initially all cells are running.
runStarts = ones(1, nCells(1));

% Pooled statistics
% In order to allow for a reduction of the number of cells in each
% simulation without loss of accuracy of the model, some parameters are
% calculated by pooling all the data points of the last n timesteps, where
% n is determined from the timeWindow defined below.

% Time window (s) for the accumulated distributions:
timeWindow = 10;

% Time pooled x positions for the last n iterations to compute
% bacterial distributions.
pooledRho = nan(timeWindow/timeStep, nCells);

% Cell displacements vars:
xDisplacement = nan(nCells,1); 
yDisplacement = nan(nCells,1);

% Bacterial distributions matrix
bacterialDistribution = nan(numel(rBins)-1, numel(intTimes));
bacterialDistribution(:,1) = ...
    histcounts(rho, rBins)'./...
    BackgroundAbundance2D./...
    diff(pi*rBins.^2)';

% Average nutrient exposure in microM/cell/s.
nutrientExposure = nan(maxIteration,1);
nutrientExposure(1) = mean(gradientFunction(rho))*timeStep;

% Reorientation modes of the motility pattern.
% Number of reorientation modes
nTumbleModes = numel(tumbleAnglesDistribution); 
tumbleMode = randi([1 nTumbleModes],1,nCells(1));

%% Create the exponential distribution function for tumbling time
% The exponential distribution is truncated to a minimum tumbling time
% equal to timestep:
mu = meanTumbleTime-timeStep;
if mu<=0
    tumbleDistributions = makedist('exponential', meanTumbleTime);
else
    tumbleDistributions = makedist('exponential', meanTumbleTime);
    tumbleDistributions = truncate(tumbleDistributions, timeStep, inf);
end

%% Main timestepping loop.

%    At any given time step each bacteria in the population is categorized
%    by the following logical variables:
%
%       runnnersIdx, 1 for cells presently running, 0 for cells tumbling.
%
%       newTumblersIdx, 1 for cells that have started a tumble in the
%           current timeStep.
%
%       newRunnersIdx, 1 for cells that have started a run in the current
%           timeStep.
%
%       bumpersIdx, 1 for cells currently running that have encountered a
%       boundary anytime during the run.
%
%     The position and orientation of each cell is described by variables
%     x, y, rho (radial distance), and Orientations.

for iteration = 2:maxIteration

    % Determine methylation state and tumbling probability of each cell:
    [methylationLevel, tumblingProbability] =...
        chemotacticPathwayModel(gradientFunction(rho), methylationLevel);

    % Randomly assign cells that could be starting tumbling in current
    % time step according to tumbling probability:
    newTumblersIdx = poissrnd(tumblingProbability)>0;

    % Remove from the newTumblers pool cells that have not been running for
    % longer than timeStep:
    newTumblersIdx = newTumblersIdx...
        & (iteration-runStarts)>=1;

    % Cells that are not tumbling are running:
    runnersIdx = ~newTumblersIdx & runStarts<=iteration;

    %% Tumbles

    % Assign tumble durations (in number of iterations) to cells that
    % started to tumble this time step. Drawn randomly from a truncated
    % exponential described above. Tumbles shorter than the minimum tumble
    % time are reset to minimum tumble time:
    tumbleDurations =...
        round(random(tumbleDistributions,1, sum(newTumblersIdx))/timeStep);

    % Assign the next tumbling mode to the new tumblers.
    tumbleMode(newTumblersIdx) =...
        mod(tumbleMode(newTumblersIdx)+1, nTumbleModes);
    tumbleMode(tumbleMode==0) = nTumbleModes;

    % Assign random tumble orientations to cells that started to tumble
    % this time step, according to their current tumble mode:
    for iMode = 1: nTumbleModes
        tumbleAngles = datasample(...
            tumbleAnglesDistribution{iMode},...
            sum(tumbleMode(newTumblersIdx)==iMode));
        Orientations(tumbleMode==iMode & newTumblersIdx) =...
            Orientations(tumbleMode==iMode & newTumblersIdx)...
            + tumbleAngles';
    end

    %% Runs

    % Apply rotational brownian motion component to cells actively running:
    Orientations(runnersIdx) =...
        Orientations(runnersIdx) +...
        sqrt(2*rotDiffCoeff .* timeStep) .* randn(1,sum(runnersIdx));

    % Calculate cell displacements accounting for the running speed of each
    % cell and the change in direction caused by Brownian motion:
    xDisplacement(runnersIdx) =...
        meanRunSpeed* timeStep*cos(Orientations(runnersIdx));
    yDisplacement(runnersIdx) =...
        meanRunSpeed* timeStep*sin(Orientations(runnersIdx));

    % New position of the running cells at the end of the timestep:
    X(runnersIdx) = X(runnersIdx) + xDisplacement(runnersIdx);
    Y(runnersIdx) = Y(runnersIdx) + yDisplacement(runnersIdx);
    rho = sqrt(Y.^2 + X.^2);

    % Reassign runStarts for cells that have started tumbling:
    runStarts(newTumblersIdx) = iteration+tumbleDurations;

    %% Boundary behaviour
    % Bacteria encountering either of the boundaries bounce away
    % conserving the angle of incidence.
    
    %%% Phytoplankton cell wall %%%
    ParticlesInBoundaries = rho < rMin;

    if sum(ParticlesInBoundaries)>0
        [X(ParticlesInBoundaries),...
            Y(ParticlesInBoundaries),...
            Orientations(ParticlesInBoundaries),...
            rho(ParticlesInBoundaries)] = bounce(...
            X(ParticlesInBoundaries),...
            Y(ParticlesInBoundaries),...
            xDisplacement(ParticlesInBoundaries),...
            yDisplacement(ParticlesInBoundaries),...
            rMin);
    end

    %%% External domain %%%
    ParticlesInBoundaries = rho > rMax;

    if sum(ParticlesInBoundaries)>0

        [X(ParticlesInBoundaries),...
            Y(ParticlesInBoundaries),...
            Orientations(ParticlesInBoundaries),...
            rho(ParticlesInBoundaries)] = bounce(...
            X(ParticlesInBoundaries),...
            Y(ParticlesInBoundaries),...
            xDisplacement(ParticlesInBoundaries),...
            yDisplacement(ParticlesInBoundaries),...
            rMax);
    end

    %% Bacterial distributions
    % If current simulation time is within the last timeWindow before the
    % next user-defined intTimes, add current cell positions to pooledRho.
    if any(ismember(iteration+(0:timeWindow/timeStep-1), intTimes/timeStep))
        idx = uint64(mod(iteration,timeWindow/timeStep)+1);
        pooledRho(idx,:) = rho;

        % Resets the container variables after 1 intTime has been
        % reached:
    elseif rem(iteration, intTimes)==1
        pooledRho = pooledRho*nan;
    end


    %% Assess bacterial distribution along rho at user-defined intTimes
    if ismember(iteration*timeStep,intTimes)
        [~,idx] = intersect(intTimes,iteration*timeStep);
        bacterialDistribution(:,idx+1) = ...
            histcounts(pooledRho(:), rBins)./...
            (BackgroundAbundance2D)./...
            diff(pi*rBins.^2)/...
            (timeWindow/timeStep);
    end

    % Cellwise average exposure to nutrients during this timestep:
    nutrientExposure(iteration) = mean(gradientFunction(rho))*timeStep;
end

%% Model of the Escherichia coli chemosensory circuit
    function [methylationLevel, tumblingProbability] =...
            chemotacticPathwayModel(...
            aspartateConcentration, methylationLevel)

        %% Phosphorylation

        % Methylation level dependent free energy difference:
        freeEnergyDifference = ALPHA*(M_0 - methylationLevel);
        freeEnergy = nReceptors*(freeEnergyDifference +...
            log(1+aspartateConcentration/K_I) -...
            log(1+aspartateConcentration/K_A));

        % Average kinase activity.
        activity = (1+exp(freeEnergy)).^-1;

        %% Methylation, using linear model (Kalinin et al 2009)

        % Function of dependence of methylation level:
        F = k_R*(1-activity)-k_B*activity;

        % Average methylation level:
        methylationLevel = methylationLevel + F*timeStep;

        %% Tumble probability (Cluzel et al 2000)

        % Probability state for a cell of going into tumbling:
        cwBias = 1./((aHalf./activity).^H+1);

        % Expected run length for each cell at each time step in seconds:
        runLength = meanTumbleTime./cwBias-meanTumbleTime;

        % Frequency of the Poisson process of run-tumble (Vladimirov et al
        % 2008):
        frequency = 1./runLength;

        % Probability of tumbling:
        tumblingProbability = timeStep*frequency;

    end


    function [Orientations, rho, X, Y] =...
            populateDomain(rMax, rMin, nCells)

        % Generates uniform distributions of nCells randomly oriented
        % within an annulus of small radius rMin and large radius rMax
        % Outputs are the Orientations of the cells, the radial distance
        % rho from the center of the simulation, the radial angle theta
        % and the x and y coordinates in the bidimensional space.

        Orientations = 2*pi*rand(1, nCells);
        rho = sqrt(rand(1, nCells)*(rMax^2-rMin^2)+rMin^2);
        theta = 2*pi*rand(1, nCells);
        X = rho .* cos(theta);
        Y = rho .* sin(theta);
    end


    function methylationLevel = equilibriumMethylationLevel(nutrients)

        % Calculates methylation levels at steady state assuming nutrient
        % concentration is constant for each location along
        % the nutrient gradient.

        methylationLevel =...
            ((ALPHA*M_0+log(1+nutrients/K_I)-log(1+nutrients/K_A))...
            *nReceptors-log(k_B/k_R))/(ALPHA*nReceptors);
    end

    function [X, Y, Orientations, rho] = bounce(X, Y, xDisplacement, yDisplacement, radius)
        X = X(:);
        Y = Y(:);
        xDisplacement = xDisplacement(:);
        yDisplacement = yDisplacement(:);

        % Calculate intersection of trajectory with the boundary. It is a
        % simplification of function linecirc
        slopes = yDisplacement./xDisplacement; 
        intercepts = Y-slopes.*X;
        a = 1+slopes.^2;
    	c = intercepts.^2-radius.^2;
        b = 2*slopes.*intercepts;

        xintersect = nan(2, size(a,2));
        yintersect = nan(2, size(a,2));
        for ii = 1:length(a)
            xintersect(:,ii) = roots([a(ii),b(ii),c(ii)]);
            yintersect(:,ii) = intercepts(ii)+slopes(ii).*xintersect(:,ii);
        end
        distances = sqrt((X-xintersect').^2+(Y-yintersect').^2);
        [~,minimum] = min(distances,[],2);
        xintersect = xintersect(sub2ind(size(xintersect), minimum, (1:size(xintersect,2))'));
        yintersect = yintersect(sub2ind(size(yintersect), minimum, (1:size(yintersect,2))'));

        % Get tangent at this intersection (slope/intercept form)
        m = - 1 ./ (yintersect ./ xintersect);
        b = yintersect - m .* xintersect;

        % "Reflect" outside points across the tangent line to "bounce" them
        % Equations from: https://stackoverflow.com/a/3307181/670206
        d = (X + (Y - b) .* m) ./ (1 + m.^2);

        X = (2 * d - X)';
        Y = (2 .* d .* m - Y + 2 .* b)';

        % In some cases, when the trajectory of the cell outside the
        % boundary is very close to the tangent at the intersect point with
        % the boundary, the end position after bouncing is still outside
        % the boundary, though very close to it. This can later on throw an
        % error, because if the cell remains outside in the following
        % iteration there is no intersect. Thus, in this cases the celll is
        % repositioned at the boundary.
        stillOut = vecnorm([X; Y])>rMax;
        slopes = Y(stillOut)./X(stillOut);
        X(stillOut) = sqrt(rMax^2./(1+slopes.^2));
        Y(stillOut) = slopes.*sqrt(rMax^2./(1+slopes.^2));

        % Recompute the direction of the particles that "bounced"
        Orientations = atan2(Y - yintersect', X - xintersect');
        rho = vecnorm([X; Y]);
    end
end
