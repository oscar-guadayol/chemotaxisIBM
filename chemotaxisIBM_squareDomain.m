function [beta, CI, R2, regressionTimes, betaSS, timeSS, nutrientExposure,...
    PopulationParameters,...
    bacterialDistribution, xBins] = ...
    chemotaxisIBM_squareDomain(...
    nCells, timeStep, domainWidth, gradientFunction,...
    meanTumbleTime, meanRunLength, rotDiffCoeff, tumbleAnglesDistribution,...
    runSpeeds, Options, boundaryBehaviour)

% chemotaxisIBM_squareDomain simulates chemotaxis in a population of bacteria in
% a square bidimensional space with a gradient of chemoattractant aspartate. It
% is designed to replicate an agarose-based microfluidic device that creates
% linear, stable and quiescent chemical gradients.
%
% The chemotactic behaviour is modelled using the chemosensory circuit of E.
% coli for aspartate, following Tu et al 2008 and Kalinin et al 2009.
%
% The motility pattern is defined by the rotational diffusivity coefficient of
% the cells and by their tumble time and angle, running speed and run length.
% The model allows for bimodal motility patterns, where the cells alternate
% between different reorientation behaviours.
%
% It regularly fits an exponential distribution and then uses spline smoothing
% of its parameters over time to assess whether the simulation has reached an
% operational steady state.
%
% Cells encountering one of the boundary perpendicular to the chemoattractant
% behave as described below in boundaryBehaviour parameter. Longitudinal
% boundaries are periodic, that is cells encountering them appear at the oposite
% side of the domain without change in speed or direction.
%
% Input variables:
%
%       nCells is the number of cells in the simulation.
%
%       timeStep is the length of the simulation time steps in seconds.
%
%       domainWidth is the width in microns of the square domain.
%
%       gradientFunction is a function handle to compute the ligand
%           concentration as a function of x position. The gradient can have any
%           shape but must be monotonically changing. For example:
%              gradientFunction = @(x) L0 + gradient*x;
%           where L0 is the concentration of aspartate at the source wall, and x
%           is the distance to it.
%
%       meanTumbleTime is the mean tumble time in s.
%
%       meanRunLength is the mean run length in s.
%
%       rotDiffCoeff is the rotational diffusion coefficient along the
%          long semiaxis in radians^2/s
%
%       tumbleAnglesDistribution is a cell in which each element contains a
%          vector with a population of tumble angles from which values are drawn
%          randomly in each tumble event. The number of cells marks the number
%          of reorientation modes the cells display. Each cells changes the
%          distribution of angles the on every successive reorientation. This
%          allows to model, for example, run-flick-reverse motility patterns.
%
%       runSpeeds is an nCellsX1 vector with the run speeds in microns/s for
%           each cell in the simulation.
%
%       Options is a structure with options for the steady state
%           analysis. Fields are:
%
%           regressionTimeStep is the time in seconds between computations
%               of beta. Default is 100 s.
%
%           spline_start is the initial time for the computation of steady
%               state. Skipping over initial rapidly varying region may save
%               computational time.
%
%           interrogation_iterations is a vector with the times in seconds
%               at which it checks for steady state attainment. Default is
%               2.^(1:1:30)./dt; The last value marks the end of the simulation
%               if steady state is not reached.
%
%           nknots0 is a 2 elements vector with the initial guesses for how
%               many spline knots to use when fitting splines to fitted
%               beta(1:2) vs time, probably never need to change. Default is [3
%               3];
%
%           steady_state_tol is the acceptable relative tolerance, in terms
%               of relative change per hour of the two exponential distribution
%               parameters. Default is 0.1 per hour.
%
%           maxTime is the maximum time in seconds for the simulation. If
%               empty or nonexistent, the simulations runs until steady state is
%               reached. If existing, simulation will run until maxTime, even if
%               steady state has been reached.
%
%           minTime is the minimum time in seconds for the simulation. Even
%               if steady state is reached before, the simulation will continue
%               until minTime. Default is 3600 seconds.
%
%           margin is the length in microns of the area to exclude from the
%               exponential fit to bacterial distribution in microns. This is to
%               account for the fact that in the experimental microfluidic
%               system cells near the walls can not be visualized. Default is 50
%               microns.
%
%           populationParametersBins is a vector of x positions that are used to
%               evaluate spatuale the population-level parameters.  The default
%               consists of three bins, one for the left margin, one for the
%               center of the channel, and one for the right margin.
%
%       boundaryBehaviour is a categorical coding for the behaviour of
%           cells when they encounter a wall. Values are:
%
%           0 Reverses running cells that bump into the x boundaries with
%                the same angle. Default.
%
%           1 As in Kalinin et al 2009, when a cell hits a boundary, it
%               sticks to the boundary for 1 s, and then leaves the boundary
%               with and arbitrary angle.
%
%           2 As case 1, but cells sticks for the average tumble time.
%
%           3 When a cell hits a boundary, it starts a tumble.
%
%           4 cell tries to keep running in the same direction without moving
%               until a tumble occurs.
%
% Output variables:
%
%       beta is an NX2 matrix, where N is the number of times
%           the regression was evaluated, of the estimated coefficients from
%           nonlinear regression of the exponential distribution of cells along
%           the ligand gradient. Second column corresponds to the chemotactic
%           precision length in microns.
%
%       CI is a 2X2XN matrix of the 95% confidence intervals of beta.
%
%       regressionTimes is a 1XN vector of the times in seconds at which
%           the regression was performed.
%
%       betaSS is 1X2 vector of the estimates of beta at the operational
%           steady states.
%
%       timeSS is a 1X2 vector of the operational steady state times in
%           seconds. The operational steady state is defined as relative change
%           per hour defined by options.steady_state_tol.
%
%       nutrientExposure is a vector with the average nutrient exposure in
%           microM/cell/s at each time in the simulation.
%
%       PopulationParameters is a table that contains spatially resolved
%           estimates and standard errors for the following population level
%           parameters: chemotactic velocities, random motility coefficients,
%           and chemotactic sensitivities.
%
%       bacterialDistribution is a vector of frequencies of bacteria (relative
%         to nCells) along the spatial gradient at the last iteration point.
%
%       xBins is a vector with the bin edges used to compute
%           bacterialDistribution.
%
%
% Dependencies: Statistics toolbox
%               Statistics and Machine Learning Toolbox SLMtools (D'Errico 2022)
%
% Bibliography:
%
%     Ahmed T, Stocker R. Experimental Verification of the Behavioral
%         Foundation of Bacterial Transport Parameters Using Microfluidics.
%         2008;95(9):4481–93.
%     Berg HC. E. coli in motion. New York, NY: Springer New York;
%         2004. 152 p. (Biological and Medical Physics, Biomedical Engineering).
%     D’Errico J. SLM - Shape Language Modeling [Internet].
%         [cited 2022 Sep 16]. Available from:
%         https://www.mathworks.com/matlabcentral/fileexchange/24443-slm-shape-language-modeling
%     Kalinin YV, Jiang L, Tu Y, Wu M. Logarithmic Sensing in Escherichia
%         coli Bacterial Chemotaxis. 2009;96(6):2439–48.
%     Shimizu TS, Tu Y, Berg HC. A modular gradient‐sensing network for
%         chemotaxis in Escherichia coli revealed by responses to time‐varying
%         stimuli. Molecular Systems Biology. 2010 Jan 1;6(1):382.
%     Tu Y, Shimizu TS, Berg HC. Modeling the chemotactic response of
%         Escherichia coli to time-varying stimuli. PNAS. 2008 Sep
%         30;105(39):14855–60.
%     Vladimirov N, Løvdok L, Lebiedz D, Sourjik V. Dependence of Bacterial
%         Chemotaxis on Gradient Shape and Adaptation Rate. PLOS Computational
%         Biology. 2008 Dec 19;4(12):e1000242.
%
%  Copyright (C) 2023, Oscar Guadayol oscar_at_guadayol.cat
%
%  This program is free software; you can redistribute it and/or modify it under
%  the terms of the GNU General Public License, version 3.0, as published by the
%  Free Software Foundation.
%
%  This program is distributed in the hope that it will be useful, but WITHOUT
%  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
%  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
%  details.
%
%  You should have received a copy of the GNU General Public License along with
%  this program.  If not, see <http://www.gnu.org/licenses/>.

%% Steady state analyses options
addpath SLMtools

if ~isfield(Options,'regressionTimeStep')...
        || isempty(Options.regressionTimeStep)
    Options.regressionTimeStep = 100/timeStep;
else
    Options.regressionTimeStep = Options.regressionTimeStep/timeStep;
end

if ~isfield(Options,'spline_start')...
        || isempty(Options.spline_start)
    Options.spline_start = 150;
end
if  ~isfield(Options,'interrogation_iterations') ...
        ||isempty(Options.interrogation_iterations)
    Options.interrogation_iterations = 2.^(1:1:20)./timeStep;
end

if isfield(Options,'maxTime') && ~isempty(Options.maxTime)
    maxIteration = Options.maxTime/timeStep;
    Options.interrogation_iterations = ...
        [Options.interrogation_iterations(...
        Options.interrogation_iterations<maxIteration)...
        maxIteration];
    checkSteadyState = false;
else
    maxIteration = Options.interrogation_iterations(end);
    checkSteadyState = true;
end

if ~isfield(Options,'options') || isempty(Options.options)
    Options.nknots0 = [3 3];
end

if ~isfield(Options,'steady_state_tol')...
        || isempty(Options.steady_state_tol)
    Options.steady_state_tol = 0.1;
end

if ~isfield(Options,'minTime') || isempty(Options.minTime)
    Options.minTime = 60*60;
end

if ~isfield(Options,'margin') || isempty(Options.margin)
    margin = 50; % margins width in micrometers
else
    margin = Options.margin;
end

if ~isfield(Options,'populationParametersBins') ||...
        isempty(Options.populationParametersBins)
    populationParametersBins = [0, margin, domainWidth-margin, domainWidth];
end

warning('off','optimlib:lsqlin:WillBeRemoved');  % for slm fitting, the
% default interior-point algorithm seems to ignore the spline shape options for
% some reason, but active-set will eventually be removed.

% Spline options for beta(1):
splineOptions(1) = slmset('concavedown','on',...
    'increasing','on','extrapolation','NaN',...
    'interiorknots','fixed','knots',Options.nknots0(1),...
    'LsqlinAlgorithm','interior-point');

% Spline options for beta(2):
splineOptions(2) = slmset('concaveup','on',...
    'decreasing','on','extrapolation','NaN',...
    'interiorknots','fixed','knots',Options.nknots0(2),...
    'LsqlinAlgorithm','interior-point');

% Declare variable containing spline fits
splineFit = nan(2,1);


%% Boundary behaviour
if ~exist('boundaryBehaviour','var') || isempty('boundaryBehaviour')
    boundaryBehaviour = 0;
end

%% Spatial domain
% Removes the margins of the gradient from the distribution in the nlinfit
% regression later on:

% Width of the bins for assess bacterial distribution:
dx = 1;

% Spatial bins used to evaluate bacterial distributions:
xBins = margin:dx:domainWidth-margin; 
midXbins = xBins(1:end-1)+diff(xBins)/2;

%% Model for nonlinear fitting of the bacterial distribution
modelfun = @(b,x)(b(1)*exp(-x/b(2)));

%% x and y positions and orientations of initial population

% Populates the spatial domain with NCells of random position and orientation:
x = 0 + domainWidth*rand(1, nCells);
y = 0 + domainWidth*rand(1, nCells);
Orientations = 2*pi*rand(1, nCells);

%% Chemotactic pathway constants
% From Kalinin et al 2009, Tu et al 2008, Shimizu et al 2010.


%%% Phosporilation constants %%%

% Number of responding receptor dimers in the MCP (Methyl-accepting Chemotactic
% Proteins):
nRECEPTORS = 6;

% Dissociation constant of the ligand to the inactive receptor in microM:
K_I = 18;

% K_A = K_I/C dissociation constant of the ligand to the active receptor:
K_A = K_I/0.0062;

% Free-energy change per added methyl group (in units of kT). In Shimizu et al
% (2010) alpha = 2kT at ~22C.
ALPHA = 1.7;

% Methylation level at which the methylation level dependent
% freeEnergyDifference is 0.
M_0 = 1;


%%% Methylation constants %%%

% Methylation rate for the inactive receptors in s^-1:
K_R = 0.005;

% De-methylation rate for the active receptors in s^-1:
K_B = 0.005;


%%% Tumble probability constants %%%

% Hill coefficient:
H = 10.3;

% Steady-state activity in the absence of stimuli. In Shimizu et al (2010) a_0 =
% 1/3 at 32C.
a_0 = 1/2;

% Clockwise bias in the absence of a chemoattactant gradient:
cwBiasHomogeneous = meanTumbleTime/(meanTumbleTime+meanRunLength);

% Kinase activity that induces a clockwise bias of 0.5 in the absence of a
% chemoattactant gradient:
aHalf = a_0*(1/cwBiasHomogeneous-1)^(1/H);

% Mean tumble time in seconds
meanTumbleTime = single(meanTumbleTime);

%%% Initial methylation levels %%%

% Calculates methylation levels M assuming all cells are at steady state for
% each location along the gradient in L and assigns the correct M level at each
% particle according to its position in x
ligandConc = gradientFunction(x);
methylationLevel =...
    ((ALPHA*M_0+log(1+ligandConc/K_I)-log(1+ligandConc/K_A))...
    *nRECEPTORS-log(K_B/K_R))/(ALPHA*nRECEPTORS);

%% Preallocation of output parameters

% Time to steady state, defined as the time at which the relative change in the
% beta parameters of the nlinfit is 
timeSS = nan(1,2);
betaSS = nan(1,2);
regressionTimes = double(...
    (Options.regressionTimeStep:Options.regressionTimeStep:maxIteration)...
    *timeStep)';
if isempty(regressionTimes)
    regressionTimes = maxIteration*timeStep;
end

nRegressions = numel(regressionTimes);
beta = nan(nRegressions,2);
CI =  nan(2,2,nRegressions);
R2 =  nan(nRegressions,1);

% Population level parameters 
% These are estimated with some degree of spatial discrimination, to minimize
% biases near the walls, and to be able to explore dependences of these
% parameters on local gradients. x bins within which the population parameters
% are estimated.
chemotacticVelocity = nan(nRegressions,numel(populationParametersBins));
chemotacticSensitivity = nan(nRegressions,numel(populationParametersBins));
motilityCoefficient = nan(nRegressions,numel(populationParametersBins));
chemotacticVelocitySE = nan(nRegressions,numel(populationParametersBins));
chemotacticSensitivitySE = nan(nRegressions,numel(populationParametersBins));
motilityCoefficientSE = nan(nRegressions,numel(populationParametersBins));
nutrientExposure = mean(gradientFunction(x))*timeStep;

%% Set nlinfit warnings to temporarily issue errors (exceptions)
s = warning('error', 'stats:nlinfit:ModelConstantWRTParam');
warning('error', 'MATLAB:nearlySingularMatrix')
warning('error', 'MATLAB:rankDeficientMatrix')
warning('error', 'stats:nlinfit:UnableToDecreaseSSE')

%% Change random number generator to 'simdTwister':
% SIMD-oriented Fast Mersenne Twister
rng(1,'simdTwister')
% rng('default') %reverts changes

%% Declare variables for IBM simulation

% nCells long vector with the iterations in which the last run started for each
% cell. If for a particular cell it is lower or equal than the current iteration
% value, then the cell is actively running. Otherwise the cell is tumbling. To
% start with, all cells are running.
runStarts = ones(1,nCells);

% Matrix of size nCellsX2 with the pooledulated sines and cosines for the
% current or last run.
averageRunOrientations = [sin(Orientations); cos(Orientations)];

% Vector of length nCells with the pooled x positions for the current or
% last run of each cell .
averageRunPositions = nan(1, nCells);

hasBumped = false(1, nCells);

% In order to reduce the number of cells in each simulation without reducing the
% accuracy of the model, some parameters are calculated by pooling all the data
% points of the las N timesteps, where n is determined from the timeWindow
% defined below.

% Time window (s) for the pooledulated distributions:
timeWindow = 60;

% Time pooled x positions for the last n iterations to compute bacterial
% distributions.
pooledX = nan(timeWindow/timeStep, nCells);

% Time-pooled run lengths and orientations used to calculate population-level
% chemotactic parameters positions.
pooledRunLengths = cell(timeWindow/timeStep,1);
pooledRunOrientations = cell(timeWindow/timeStep,1);
pooledxRuns = cell(timeWindow/timeStep,1);

% Randomly assign initial tumble modes for each cells:
nTumbleModes = numel(tumbleAnglesDistribution);
tumbleMode = randi(nTumbleModes,1,nCells(1));

% Make sure runSpeeds is a row vector;
if ~isrow(runSpeeds), runSpeeds = runSpeeds';end

%% Create the exponential distribution function for tumbling probability.
% The exponential distribution is truncated to a minimum tumbling time equal to
% timestep:
mu = meanTumbleTime-single(timeStep);
if mu<=0
    mu=meanTumbleTime;
end
tumbleDistributions = makedist('exponential', mu);
tumbleDistributions = truncate(tumbleDistributions, timeStep, inf);


%% Symbolic functions to calculate chemotactic parameter with uncertainties

% Equations to calculate chemotactic parameters from individual tracks are taken
% from Ahmed and Stocker 2008.

%%% Equations constants %%%

% KD is the receptor/ligand dissociation constant in microMol For E. coli
% exposed to alpha-methylaspartate, KD has been estimated as 0.125 mM or 0.160
% mM:
KD = 0.125*1000;

% Mean and standard error run speeds (microns/s)
swimmingSpeed = mean(runSpeeds(:),'omitnan');
swimmingSpeedSE = std(runSpeeds(:),'omitnan')/...
    sqrt(sum(~isnan(runSpeeds(:))));

% Mean and standard error directional persistence
directionalPersistence =...
    mean(cos(cat(1,tumbleAnglesDistribution{:})),'omitnan');
directionalPersistenceSE =...
    std(cos(cat(1,tumbleAnglesDistribution{:})),'omitnan')/...
    sqrt(sum(~isnan(cos(cat(1,tumbleAnglesDistribution{:})))));


%%% Equations variables %%%

% Create symbolic variables into the static workspace for analytical solution of
% the equations
sChemotacticSensitivity = sym('sChemotacticSensitivity');
eChemotacticSensitivity = sym('eChemotacticSensitivity');
sMotilityCoefficient = sym('sMotilityCoefficient');
eMotilityCoefficient = sym('eMotilityCoefficient');
sDirectionalPersistence = sym('sDirectionalPersistence');
eDirectionalPersistence = sym('eDirectionalPersistence');
sChemotacticVelocity = sym('sChemotacticVelocity');
eChemotacticVelocity = sym('eChemotacticVelocity');
sSwimmingSpeed = sym('sSwimmingSpeed');
sTPlus = sym('sTPlus');
sTMinus = sym('sTMinus');
eSwimmingSpeed = sym('eSwimmingSpeed');
eTPlus = sym('eTPlus');
eTMinus = sym('eTMinus');
sMeanC = sym('sMeanC');
sGradient = sym('sGradient');


%%% Equations %%%

sChemotacticVelocity =...
    8*sSwimmingSpeed/3/pi*(sTPlus-sTMinus)/(sTPlus+sTMinus);
eChemotacticVelocity = sqrt(...
    eSwimmingSpeed^2*diff(sChemotacticVelocity,sSwimmingSpeed)^2 +...
    eTPlus^2*diff(sChemotacticVelocity,sTPlus)^2 +...
    eTMinus^2*diff(sChemotacticVelocity,sTMinus)^2);

sChemotacticSensitivity =...
    (atanh(3*pi*sChemotacticVelocity/8/sSwimmingSpeed))/...
    (pi/8/sSwimmingSpeed*KD/(KD+sMeanC)^2*(-sGradient))*1e-8;
eChemotacticSensitivity = sqrt(...
    eSwimmingSpeed^2*diff(sChemotacticSensitivity,sSwimmingSpeed)^2 +...
    eTPlus^2*diff(sChemotacticVelocity,sTPlus)^2 +...
    eTMinus^2*diff(sChemotacticVelocity,sTMinus)^2);

sMotilityCoefficient =...
    16*sSwimmingSpeed.^2*sTMinus/3/pi^2./(1-sDirectionalPersistence);
eMotilityCoefficient = sqrt(...
    eSwimmingSpeed^2*diff(sMotilityCoefficient,sSwimmingSpeed)^2 +...
    eTMinus^2*diff(sMotilityCoefficient,sTMinus)^2 +...
    eDirectionalPersistence^2*...
    diff(sMotilityCoefficient,sDirectionalPersistence)^2);

%% IBM SIMULATION
%    At any given time step each bacteria in the population is categorized by
%    the following logical variables:
%
%      runnnersIdx, 1 for cells presently running, 0 for cells tumbling.
%
%      newTumblersIdx, 1 for cells that have started a tumble in the
%       current timeStep.
%
%      newRunnersIdx, 1 for cells that have started a run in the current
%       timeStep.
%
%      bumpersIdx, 1 for cells currently running that have bumped
%       anytime during the run.
%
%   The position and orientation of each cell is described by variables x,
%      y and Orientations.

for iteration = 2:maxIteration

    % Determine methylation state and tumbling probability of each cell:
    [methylationLevel, tumblingProbability] =...
        chemotacticPathwayModel(gradientFunction(x), methylationLevel);

    % Randomly assign cells that could be starting tumbling in current time step
    % according to tumbling probability:
    newTumblersIdx = poissrnd(tumblingProbability)>0;

    % Remove from the newTumblers pool cells that have not been running for
    % longer than timeStep:
    newTumblersIdx = newTumblersIdx...
        & (iteration-runStarts)>=1;

    % Cells that are not tumbling are running:
    runnersIdx = ~newTumblersIdx & runStarts<=iteration;

   
    %% Tumbles
    % Assign random tumble durations (in number of iterations) to cells that
    % started to tumble this time step. Tumbles shorter than the minimum tumble
    % time are reset to minimum tumble time:
    tumbleDurations =...
        round(random(tumbleDistributions,1, sum(newTumblersIdx))/timeStep);

    % Assign the next tumbling mode to the new tumblers.
    tumbleMode(newTumblersIdx) =...
        mod(tumbleMode(newTumblersIdx)+1, nTumbleModes);
    tumbleMode(tumbleMode==0) = nTumbleModes;

    % Assign random tumble orientations to cells that started to tumble this
    % time step, according to their current tumble mode:
    for iMode = 1: nTumbleModes
        iModeTumblesIdx = tumbleMode==iMode & newTumblersIdx;
        tumbleAngles = datasample(...
            tumbleAnglesDistribution{iMode},...
            sum(iModeTumblesIdx));
        Orientations(iModeTumblesIdx) =...
            Orientations(iModeTumblesIdx) + tumbleAngles;
    end


    %% Runs
    % Apply rotational brownian motion component only to cells actively running:
    Orientations(runnersIdx) =...
        Orientations(runnersIdx) +...
        sqrt(2*rotDiffCoeff .* timeStep) .* randn(1,sum(runnersIdx));

    % Calculate cell displacements accounting for the running speed of each cell
    % and the change in direction caused by Brownian motion:
    xDisplacement =...
        runSpeeds(runnersIdx) .* timeStep .*cos(Orientations(runnersIdx));
    yDisplacement =...
        runSpeeds(runnersIdx) .* timeStep .*sin(Orientations(runnersIdx));

    % Average nutrient exposure per cell during this timestep. It assumes that
    % both the trajectory of the cell and the gradient is lineal, which it may
    % only be true for very short timesteps.
    nutrientExposure(iteration) = mean(gradientFunction(x))*timeStep;

    % New position of the running cells at the end of the timestep:
    x(runnersIdx) = x(runnersIdx) + xDisplacement;
    y(runnersIdx) = y(runnersIdx) + yDisplacement;

    % Cells crossing upper or lower boundaries reappear in the opposite side
    % (periodic boundary condition)
    y(y>domainWidth) = y(y>domainWidth)-domainWidth;
    y(y<0) = y(y<0)+domainWidth;

    % Cells that started a tumble and have not bumped during the previous run,
    % excluding the runs that were started before the simulation (i.e.,
    % runStarts>1):
    completedRunsIdx = newTumblersIdx & ~hasBumped & runStarts>1;

    % Compute the length in number of iterations of the previous runs for the
    % cells that started a tumble and have not bumped during the previous run
    % excluding the runs that were started before the simulation:
    runLengths = iteration-runStarts(completedRunsIdx);

    % Reassign runStarts to cells that have started tumbling:
    runStarts(newTumblersIdx) = iteration+tumbleDurations;

    %% Boundary behaviour
    % For cells that have bumped agaionst the wall during this timestep, execute
    % the preset behaviour

    switch boundaryBehaviour
        case 0

            % Reverses running cells that bump into the x boundaries with the
            % same angle.

            bumpersIdx =...
                all([runnersIdx; any([x<0; x>domainWidth],1)], 1);
            Orientations(bumpersIdx) = Orientations(bumpersIdx)+pi;
            x(x<0) = 0-x(x<0);
            x(x>domainWidth) = domainWidth*2-x(x>domainWidth);

        case 1

            % As in Kalinin et al 2009, when a cell hits a boundary, it sticks
            % to the boundary for 1 s (bumpDurations = 1/timeStep), and then
            % leaves the boundary with and arbitrary angle.

            bumpersIdx =...
                all([runnersIdx; any([x < 0; (x > domainWidth)],1)], 1);
            Orientations(runnersIdx & x < 0) =...
                pi*rand(sum(runnersIdx & x < 0),1)-pi/2;
            Orientations(runnersIdx & x > domainWidth) =...
                pi*rand(sum(runnersIdx & x > domainWidth),1)+pi/2;
            x(x<0) = 0-x(x<0);
            x(x>domainWidth) = domainWidth*2-x(x>domainWidth);
            bumpDurations = 1/timeStep;

            % Reassign run_starts to cells that have just bumped into a
            % boundary:
            runStarts(bumpersIdx) = iteration+1/bumpDurations+1;

        case 2

            % As case 1, but cells sticks for the average tumble time.

            bumpersIdx =...
                all([runnersIdx; any([x < 0; (x > domainWidth)],1)], 1);
            Orientations(runnersIdx & x < 0) =...
                pi*rand(sum(runnersIdx & x < 0),1)-pi/2;
            Orientations(runnersIdx & x > domainWidth) =...
                pi*rand(sum(runnersIdx & x > domainWidth),1)+pi/2;
            x(x<0) = 0-x(x<0);
            x(x>domainWidth) = domainWidth*2-x(x>domainWidth);
            bumpDurations = meanTumbleTime/timeStep;

            % Reassign run_starts to cells that have just bumped into a
            % boundary:
            runStarts(bumpersIdx) = iteration+1/bumpDurations+1;

        case 3

            % When a cell hits a boundary, it starts a tumble

            bumpersIdx =...
                all([runnersIdx; any([x < 0; (x > domainWidth)],1)], 1);
            Orientations(runnersIdx & x < 0) =...
                pi*rand(sum(runnersIdx & x < 0),1)-pi/2;
            Orientations(runnersIdx & x > domainWidth) =...
                pi*rand(sum(runnersIdx & x > domainWidth),1)+pi/2;
            x(x<0) = 0-x(x<0);
            x(x>domainWidth) = domainWidth*2-x(x>domainWidth);
            bumpDurations =...
                round(exprnd(meanTumbleTime/timeStep,sum(bumpersIdx),1));
            bumpDurations(bumpDurations<0.1/timeStep) =...
                timeStep/timeStep;

            % Reassign runStarts to cells that have just bumped into a boundary:
            runStarts(bumpersIdx) = iteration+1/bumpDurations+1;

        case 4

            % Particle tries to keep running in the same direction without
            % moving, until a tumble occurs.

            bumpersIdx = [];
            x(x<0) = 0;
            x(x>domainWidth) = domainWidth;
    end


    %% Resets state of bumped cells

    % Cells that have bumped during this iteration are removed from the run
    % statistics:
    hasBumped(bumpersIdx) = true;

    % Cells that start a run are reinstated in the run statistics
    hasBumped(runStarts==iteration) = false;


    %% Model statistics
    % This is to collect the data to compute chemotactic parameters from the
    % theoretical model.

    % Records the increase in sine and cosine of the present time step for cells
    % actively running. If this is the first time step of the run, it resets the
    % runOrientations.
    newRunners = runnersIdx &...
        runStarts==iteration &...
        runStarts>1;
    oldRunners = runnersIdx &...
        runStarts<iteration &...
        ~hasBumped &...
        runStarts>1;
    averageRunOrientations(:,newRunners) =...
        [sin(Orientations(newRunners)); cos(Orientations(newRunners))];
    averageRunOrientations(:,oldRunners) =...
        averageRunOrientations(:,oldRunners)...
        + [sin(Orientations(oldRunners)); cos(Orientations(oldRunners))];

    averageRunPositions(newRunners) = x(newRunners);
    averageRunPositions(oldRunners) = ...
        averageRunPositions(oldRunners) + x(oldRunners);

    % Registers x positions for the computation of bacterial distribution only
    % for the last timeWindow/timeStep times before the next beta estimation.
    if any(~rem(iteration+(1:timeWindow/timeStep),...
            Options.regressionTimeStep)) &&...
            ~rem(iteration,0.1/timeStep)
        idx = uint64(mod(floor(iteration*timeStep/0.1),timeWindow/0.1)+1);
        pooledX(idx,:) = x;
        pooledRunLengths{idx} = uint16(runLengths);

        % Computes average orientation of the previous run only for those cells
        % that have tumbled this timestep:
        pooledRunOrientations{idx} = atan2(...
            averageRunOrientations(1,completedRunsIdx),...
            averageRunOrientations(2,completedRunsIdx));
        pooledxRuns{idx} = ...
            averageRunPositions(completedRunsIdx)./runLengths;

        % Resets the container variables after 1 regressionTimeStep has been
        % reached:
    elseif rem(iteration, Options.regressionTimeStep)==1
        pooledX = pooledX*nan;
        pooledRunLengths = cell(timeWindow/0.1,1);
        pooledRunOrientations = cell(timeWindow/0.1,1);
        pooledxRuns = cell(timeWindow/0.1,1);
    end


    %% Computes parameters of exponential distribution:

    if ~rem(iteration, Options.regressionTimeStep)

        %%% nonlinear fitting %%%
        bacterialDistribution = histcounts(pooledX(:),xBins)';
        bacterialDistribution =...
            bacterialDistribution./sum(bacterialDistribution)/dx;
        idx = iteration/Options.regressionTimeStep;
        try
            % Initial estimates of the exponential distribution parameters using
            % linear regression:
            beta0 = regress(...
                log(bacterialDistribution(bacterialDistribution>0)),...
                [ones(size(midXbins(bacterialDistribution>0),1),1),...
                midXbins(bacterialDistribution>0)]);

            % Nonlinear fit of the cell distribution along the gradient.
            [beta(idx,:),R,~,CovB] =...
                nlinfit(midXbins,...
                bacterialDistribution,...
                modelfun,...
                [exp(beta0(1)) -1/beta0(2)]);

            % Confidence intervals of the parameters
            CI(:,:,idx) = nlparci(beta(idx,:),R,'covar',CovB);

            % Calculates the R2 as a crude good of fitness parameter.
            yHat = beta(idx,1)*exp(-midXbins/beta(idx,2));
            ESS =... % Explained sum of squares.
                sum((yHat-mean(bacterialDistribution)).^2);
            TSS =... %Total sum of squares.
                sum((bacterialDistribution-mean(bacterialDistribution)).^2);
            R2(idx) = ESS/TSS;

            % If any of the CI is lower than 0, that is if the parameters are
            % statistically non different from 0, then disregard the estimates.
            if any(CI(:,1,idx)<0) || any(any(isnan(CI(:,:,idx))))
                beta(idx,:) = nan;
                CI(:,:,idx) = nan;
                R2(idx) = nan;
            end

        catch
            beta(idx,:) = nan;
            CI(:,:,idx) = nan;
            R2(idx) = nan;
        end


        %%% Chemotactic parameters %%%

        % Computes population chemotactic parameters from the individual
        % motility parameters (run speed and length, tumble time and angle)
        % using equations in Ahmed and Stocker 2008. This only makes sense when
        % the gradient is linear.
        pooledRunOrientations = double(horzcat(pooledRunOrientations{:}))';
        pooledRunLengths = double(horzcat(pooledRunLengths{:}))'*timeStep;
        pooledxRuns = double(horzcat(pooledxRuns{:}))';


        [~,~,runsInBin] = histcounts(pooledxRuns, populationParametersBins);
        for iBin = 1:max(runsInBin)
            meanLigandGrad = diff(...
                gradientFunction(populationParametersBins(iBin:iBin+1)));
            meanLigandConc = integral(...
                gradientFunction,...
                populationParametersBins(iBin),...
                populationParametersBins(iBin+1))/...
                diff(populationParametersBins(iBin:iBin+1));
            [chemotacticVelocity(idx,iBin),...
                chemotacticVelocitySE(idx,iBin),...
                chemotacticSensitivity(idx, iBin),...
                chemotacticSensitivitySE(idx, iBin),...
                motilityCoefficient(idx, iBin),...
                motilityCoefficientSE(idx, iBin)] =...
                chemotactic_parameters(...
                pooledRunOrientations(runsInBin==iBin),...
                pooledRunLengths(runsInBin==iBin),...
                swimmingSpeed, swimmingSpeedSE,...
                directionalPersistence, directionalPersistenceSE,...
                meanLigandGrad,...
                meanLigandConc);
        end

    end
    if checkSteadyState

        %% Assess if steady state has been reached

        % Check if current time is an interrogation time:
        ind = find(Options.interrogation_iterations==iteration,1);

        if ~isempty(ind) % && ~exist('maxTime','var')

            % Remove the early transient part of the beta(time) curve from the
            % spline fit:
            filter = regressionTimes >= Options.spline_start;

            % Remove nan tails from the vectors for the spline fit:
            filter = filter(:) & ~isnan(beta(:,2));

            % Check there are at least 4 beta(time) data points. If not, skip
            % the steady state assessment this time.
            if sum(filter) < 4
                continue
            end

            % Initialize best estimate of steady state beta:
            betaSS = [NaN NaN];
            for bi = 1:2 % for beta(1:2)
                splineOptions(bi).Knots = linspace(...
                    log10(regressionTimes(find(filter,1))),...
                    log10(regressionTimes(find(filter,1,'last'))),...
                    max([5 floor(sum(filter)/20)]));
                splineFit(bi) = slmengine(log10(regressionTimes(filter)),...
                    log10(beta(filter,bi)), splineOptions(bi)); % actually does the spline fit

                try
                    splineFitTime =...
                        10^splineFit(bi).knots(1):1:10^splineFit(bi).knots(end);
                    splineFitbeta =...
                        10.^slmeval(log10(splineFitTime), splineFit(bi));
                    relativeChange =...
                        abs(gradient(splineFitbeta,splineFitTime))./splineFitbeta.*3600; % relative change in hours^-1
                    timeSS(bi) = splineFitTime(...
                        find(relativeChange>Options.steady_state_tol,1,'last'));
                    betaSS(bi) = 10.^slmeval(log10(timeSS(bi)), splineFit(bi));
                catch
                    continue
                end
            end

            % Ensure steady state was not reached within the last 10 minutes, so
            % there is no false positives, and that the simulation has surpassed
            % the minimum time before interrupting the timestepping loop:
            if all(timeSS < Options.interrogation_iterations(ind)*timeStep-600)...
                    && iteration*timeStep>Options.minTime
                % we reached steady state or maximum simulation time, stop
                % simulation!
                break
            end
        end
    end
end  % Main timestepping loop

if checkSteadyState
    if any(timeSS > iteration*timeStep-600)||...
            any(isnan(timeSS))
        disp('Warning - steady state not reached.')
    end
end

%% Restore the warnings back to their previous (non-error) state
warning(s)


%% Population parameters
% To simplify the output of the function, creates a single table with all the
% population chemotactic parameters.
PopulationParameters = table(...
    regressionTimes, chemotacticVelocity, chemotacticVelocitySE,...
    chemotacticSensitivity, chemotacticSensitivitySE,...
    motilityCoefficient, motilityCoefficientSE);


%% Remove empty rows and reorganize output data
lastRow = iteration/Options.regressionTimeStep;
beta(lastRow+1:end,:) = [];
nutrientExposure(iteration+1:end) = [];
CI = permute(CI,[3,1,2]);
CI(lastRow+1:end,:,:) = [];
CI = permute(CI,[2,3,1]);
R2(lastRow+1:end) = [];
PopulationParameters(lastRow+1:end,:) = [];
regressionTimes(lastRow+1:end) = [];


%% SUBFUNCTIONS %%

%% Model of Escherichia coli chemosensory circuit
    function [methylationLevel, tumblingProbability] =...
            chemotacticPathwayModel(...
            aspartateConcentration, methylationLevel)

        %% Phosphorylation

        % Methylation level dependent free energy difference:
        freeEnergyDifference = ALPHA*(M_0 - methylationLevel);
        freeEnergy = nRECEPTORS*(freeEnergyDifference +...
            log(1+aspartateConcentration/K_I) -...
            log(1+aspartateConcentration/K_A));

        % Average kinase activity.
        activity = (1+exp(freeEnergy)).^-1;

        %% Methylation, using linear model (Kalinin et al 2009)

        % Function of dependence of methylation level:
        F = K_R*(1-activity)-K_B*activity;

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

%% Compute population level chemotactic population parameters with uncertainties
    function [chemotacticVelocity, chemotacticVelocitySE,...
            chemotacticSensitivity,chemotacticSensitivitySE,...
            motilityCoefficient, motilityCoefficientSE] =...
            chemotactic_parameters(...
            runOrientations, runLengths,...
            swimmingSpeed, swimmingSpeedSE,...
            directionalPersistence, directionalPersistenceSE,...
            meanLigandGrad, meanLigandConc)

        % T is a logical variable that discriminates cells runnning upgradient
        % from those running downgradient. It assumes the chemoattractant source
        % is at the left.

        % If chemoattractant source is on the left:
        if meanLigandGrad<0
            T = cos(runOrientations)<0;
            % If chemoattractant source is on the right:
        else
            T = cos(runOrientations)>0;
        end

        % Mean and standard error run times (s) when bacteria travel up the
        % chemical gradient:
        TPlus = mean(runLengths(T));
        TPlusSE = std(runLengths(T))/sqrt(numel(runLengths(T)));

        % Mean and standard error run times (s) when bacteria travel down the
        % chemical gradient:
        TMinus = mean(runLengths(~T));
        TMinusSE = std(runLengths(~T))/sqrt(numel(runLengths(~T)));

        % Mean and standard error chemotactic speed upgradient (microns/s):
        chemotacticVelocity = double(subs(sChemotacticVelocity,...
            {sSwimmingSpeed,sTPlus,sTMinus},{swimmingSpeed,TPlus,TMinus}));
        chemotacticVelocitySE = double(subs(eChemotacticVelocity,...
            {sSwimmingSpeed, sTPlus, sTMinus,...
            eSwimmingSpeed,  eTPlus,   eTMinus},...
            {swimmingSpeed,   TPlus,  TMinus,...
            swimmingSpeedSE,  TPlusSE, TMinusSE}));

        % Mean and standard error chemotactic sensitivity in micrometre^2/s:
        chemotacticSensitivity = double(subs(sChemotacticSensitivity,...
            {sSwimmingSpeed ,sTPlus, sTMinus, sGradient, sMeanC},...
            {swimmingSpeed, TPlus, TMinus,...
            meanLigandGrad, meanLigandConc}));
        chemotacticSensitivitySE = double(subs(eChemotacticSensitivity,...
            {sSwimmingSpeed, sTPlus, sTMinus,...
            eSwimmingSpeed, eTPlus, eTMinus, sGradient, sMeanC},...
            {swimmingSpeed, TPlus, TMinus,...
            swimmingSpeedSE, TPlusSE, TMinusSE,...
            meanLigandGrad, meanLigandConc}));

        % Mean and standard error random motility coefficients in
        % micrometre^2/s. Assumes run_lengths don't decrease when swimming
        % downgradient (optimistic behaviour, Berg 2004):

        motilityCoefficient = double(subs(sMotilityCoefficient,...
            {sSwimmingSpeed, sTMinus, sDirectionalPersistence},...
            {swimmingSpeed, TMinus, directionalPersistence}));
        motilityCoefficientSE = double(subs(eMotilityCoefficient,...
            {sSwimmingSpeed, sTMinus, sDirectionalPersistence,...
            eSwimmingSpeed, eTMinus, eDirectionalPersistence},...
            {swimmingSpeed,TMinus, directionalPersistence,...
            swimmingSpeedSE,TMinusSE, directionalPersistenceSE}));
    end
end