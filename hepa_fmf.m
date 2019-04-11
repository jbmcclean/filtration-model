% Simple free molecular flow HEPA filter model

clear;

% Define bulk filter properties
mu = 1E-5; % dynamic viscosity (Pa s)
U_f = 5; % filtration velocity (cm/s)
Z = 0.400; % thickness of filter (mm)
A = 0.23*0.08; % area of filter (m^2)

% Define filter properties
muFibreDia = 1; % mean fibre diameter
sigmaFibreDia = 0.1; % standard deviation of fibre diameter
muFibreSpacing = 10; % mean fibre spacing (centre to centre)
sigmaFibreSpacing = 0.1; % standard deviation of fibre spacing (centre to centre)
nFibres = 100; % number of fibres

% Define dust particle properties
muDustParticleDia = 1; % mean dust particle diameter
sigmaDustParticleDia = 0.1; % standard deviation of dust particle diameter
nDustParticlesToCollect = 1000; % how many particles to collect?

% Create an empty array ready to record the diameter of collected particles
capturedDustParticleDia = zeros(nDustParticlesToCollect,1);
nDustParticleCollected = 0;

% Create a seed array to capture the volume of collected particles with
% timestep as well as the pressure drop
volCollectedParticles = zeros(1,2);
Delta_P = zeros(1,1);

%% Generate a range of fibre centre points

% Generate a normal distribution of fibre centre points
fibreSpacingDist = normrnd(muFibreSpacing,sigmaFibreSpacing,nFibres,1);

% Carry out a cumulative sum to generate the fibre centre point coordinates
fibreSpacingCumDist = cumsum(fibreSpacingDist);

%% Generate the fibres

% Generate a normal distribution of fibre diametres
fibreDiaDist = normrnd(muFibreDia,sigmaFibreDia,nFibres,1);

% Calculate the total volume of the fibres
V_f = sum(0.08*pi*(fibreDiaDist/2).^2);

% Calcuate the constant in the pressure drop equation (Davies 1973)
U_f = U_f*1E-2; % convert from cm/s to m/s
Z = Z*1E-3; % convert from mm to m
V = A*Z; % calcuate filter face area
d_f = muFibreDia*1E-6; % convert from um to m
Delta_P_const = (64*mu*U_f*Z*(V^(-3/2))*V_f)/d_f;

%% Generate the starting position of the incident particle

% Get the centre point of the final fibre
filterEndPos = max(fibreSpacingCumDist);

timeStep = 0;


axis([0,1000,1E7,5E7]);
hold on;

while nDustParticleCollected < nDustParticlesToCollect

timeStep = timeStep + 1;    
    
% Generate a random x-coordinate betweeen 0 and the centre point of the
% final fibre
dustParticleCentrePos = filterEndPos*rand(1,1);

%% Generate the size of the particle

% Generate the diameter of the dust particle from a normal distribution
dustParticleDia = normrnd(muDustParticleDia,sigmaDustParticleDia);

%% Determine which two fibres the particle is between

result = dustParticleCentrePos - fibreSpacingCumDist;

% Any comparisons less than zero mean that the dust particle is to the
% right of that fibre
result(result<0) = NaN;
[dustParticleGapPos,fibreLeftNum] = min(result);

% Derive the positions of the fibres
fibreLeftPos  = fibreSpacingCumDist(fibreLeftNum);
fibreRightNum = fibreLeftNum + 1;
fibreRightPos = fibreSpacingCumDist(fibreRightNum);

% Derive the diameters of the fibres
fibreLeftDia = fibreDiaDist(fibreLeftNum);
fibreRightDia = fibreDiaDist(fibreRightNum);

%% Work out the geometry of the encounter to determine whether the particle sticks or not

fibreRightPos = fibreRightPos - fibreLeftPos;
fibreLeftPos = 0;

% clf;
% hold on;
% axis([0,20,-10,10]);
% plot([0,fibreLeftDia/2],[0,0],'-+k');
% plot([fibreRightPos-fibreRightDia/2,fibreRightPos],[0,0],'-+k');
% plot(dustParticleGapPos,0,'or');
% plot([dustParticleGapPos-dustParticleDia/2,dustParticleGapPos+dustParticleDia/2],[0,0],'-+r');
% 
% % Try plotting circles
% viscircles([0,0],fibreLeftDia/2);
% viscircles([fibreRightPos,0],fibreRightDia/2);
% viscircles([dustParticleGapPos,0],dustParticleDia/2);
% hold off;

% Check for capture on leftmost fibre
if((dustParticleGapPos-dustParticleDia/2)<=fibreLeftDia/2);
    %plot(dustParticleGapPos,0.5,'*r');
    nDustParticleCollected = nDustParticleCollected + 1;
    capturedDustParticleDia(nDustParticleCollected,1) = dustParticleDia;
    fprintf('Collected particle number %d.\n',nDustParticleCollected);
end

% Check for capture on rightmost fibre
if((dustParticleGapPos+dustParticleDia/2)>=(fibreRightPos-fibreRightDia/2));
    %plot(dustParticleGapPos,0.5,'*g');
    nDustParticleCollected = nDustParticleCollected + 1;
    capturedDustParticleDia(nDustParticleCollected,1) = dustParticleDia;
    fprintf('Collected particle number %d.\n',nDustParticleCollected);
end

volCollectedParticles(end+1,1) = timeStep;
volCollectedParticles(end,2) = sum((4/3)*pi*(capturedDustParticleDia/2).^3);

V_p = volCollectedParticles(end,2);
Delta_P(end+1,1) = Delta_P_const*(V_f + V_p)^(1/2);

plot(timeStep,Delta_P(end,1),'+k');
pause(0.01)

end
