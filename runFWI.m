%% This script is used for running the Full Waveform Inversion method.
% We're plotting the model & gradient in  this script.
%% Setting up dimensions
dims.dy =     10; % [m]
dims.dx =     10; % [m]
dims.dt = 1.0e-3; % [s]

dims.ny = 201; % Cells in y-direction
dims.nx = 301; % Cells in x-direction
dims.nt = 801; % Amount of time steps

%% Model dimensions
dims.modely = 100:150;
dims.modelx = 100:200;
dims.my = length(dims.modely);
dims.mx = length(dims.modelx);

%% Source locations
sx = min(dims.modelx):max(dims.modelx);
sy = min(dims.modely)*ones(1,length(sx));
dims.srcPos = sy + dims.ny*sx;

%% Receiver locations
rx = min(dims.modelx):max(dims.modelx);
ry = min(dims.modely)*ones(1,length(rx));
dims.recPos = ry+dims.ny*rx;

%% Creating background model
bg = zeros(dims.ny,dims.nx,'single');
bg(:) = 2.0e3;         % [m/s] - Background
bg(115:end,:) = 2.3e3; % [m/s] - Layer

%% Begin iteration
model = bg;     % Starting model
dims.ds = 5;   % Grid point distance between sources
maxIter = 10;   % Maximum number of iterations per frequency
freqs = [4,6,8,10,12];  % Frequencies to use in inversion

%% Initial value of u(x,t)
unext = zeros(dims.ny,dims.nx,'single');
u = zeros(dims.ny,dims.nx,'single');
uprev = zeros(dims.ny,dims.nx,'single');

errVec = zeros(1,maxIter*length(freqs));
it = 1; 
start=tic;
fprintf('Iteration \t Error \t\t Runtime(s) \n');
%% FWI loop
for f = freqs
    %% Generating ricker source signature wavelet 
    source = rickerWave(f,dims);
    %% Load true recording
    load (['trueRec_',num2str(f),'Hz.mat']);   
    for i = 1:maxIter
        %% Calculate gradient
        [gradient,err] = calculateGradient(dims,source,trueRec,model,unext,u,uprev);
        %% Taper gradient
        taperg = taperGradient(gradient);
        %% Calculate error & step length
        [stepLength,err] = calculateStepLength(dims,model,taperg,err,source,trueRec,unext,u,uprev);
        %% Update model
        model = model+stepLength*taperg;
        %% Gradient & Model plotting
        figure(2)
        subplot(2,1,1);
        imagesc(taperg(dims.modely,dims.modelx));
        colormap(jet);
        title(['Gradient ',num2str(f),' Hz Plot']);
        axis('image');
        colorbar();
        drawnow();
        subplot(2,1,2);
        imagesc(model(dims.modely,dims.modelx));
        colormap(jet);
        title(['Model ',num2str(f),' Hz Plot']);
        axis('image');
        colorbar();
        drawnow();
        
        errVec(it) = err;
        it = it + 1;
        runtime=toc(start);
        fprintf('%i \t\t\t %6.4f \t %6.4f \n',i,err,runtime);
    end
end