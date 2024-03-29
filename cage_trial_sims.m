% Script for running simulations of the (updated) D. suzukii cage trial
% code in cage_trial_split.m and/or cage_trial_full.m. These simulations 
% use the most recent set of parameters (01/22/23). 
%
% Also generates figures for showing average gen. of extinction. 

%% Test sims

% conversion efficiencies
%            DRIVE 1 (dominant)     DRIVE 2 (recessive)
%  MALES          95.03%                  97.36%
% FEMALES          nan                    93.88%
multiRelease = false;
rho = 1;
MALE_CONV_RATE = 0.9736; % 0.9503; 
FEMALE_CONV_RATE = 0.9388; % 
fitnessCostVec = [0.05, 0];
RELATIVE_FECUNDITY = 0.6274;
graphBool = true; 

% cage_trial_split for split drive simulations, cage_trial_full for
% autonomous drive simulations
data = cage_trial_full(multiRelease,rho,MALE_CONV_RATE,FEMALE_CONV_RATE,...
    fitnessCostVec,RELATIVE_FECUNDITY,graphBool);
disp(data.extinctGens)
disp(data.femaleMat)
disp(data.maleMat)

%% Plots of generation of extinction for drive type 1 (dominant female sterile)

% release ratio vector
rhoVec = 0.05:0.05:1;
% no. of replicates per experiment
numReps = 100;
% store generation of extinction by driveType...
extinctGenMat = nan(numReps,length(rhoVec),2);
% ...and generation of extinction for no drive
extinctGenMat_noDrive = nan(numReps,length(rhoVec),2);

% parameters for simulation
multiRelease = true;
MALE_CONV_RATE = 0.9503;
FEMALE_CONV_RATE = NaN;
fitnessCostVec = [0.05, 0];
RELATIVE_FECUNDITY = 0.6274;
graphBool = false; 

for i = 1:length(rhoVec)
    sprintf("Running sim %.0f of %.0f",i,length(rhoVec))
    % setup experiment
    rho = rhoVec(i); 
    for j = 1:numReps
        data = cage_trial_split(multiRelease,rho,MALE_CONV_RATE,FEMALE_CONV_RATE,...
            fitnessCostVec,RELATIVE_FECUNDITY,graphBool);
        extinctGenMat(j,i,1) = data.extinctGens;
    end
end

% plot driveType 1 results
plot(rhoVec, mean(extinctGenMat(:,:,1)),'-o','LineWidth',2,'color','black',...
    'MarkerFaceColor','black');
% title({'avg. generation of extinction for', 'dominant female sterile drive and',...
%     'multiple releases'},...
%     'fontweight','bold')
xlim([0,max(rhoVec)+0.05]);
ylim([0,max(extinctGenMat(:,:,1),[],'all')]);
ylabel('avg. generation of extinction','interpreter','latex');
xlabel('release ratio','interpreter','latex');
set(gca,'FontSize',16);
set(gca,'Layer','top')
f = gcf;
exportgraphics(f,'dominant_batch_runs.pdf','Resolution',600);
% hold on

%% Plots of generation of extinction for drive type 2 (recessive female sterile)

% parameters for simulation
MALE_CONV_RATE = 0.9736;
FEMALE_CONV_RATE = 0.9388;

for i = 1:length(rhoVec)
    sprintf("Running sim %.0f of %.0f (drive Type 2)",i,length(rhoVec))
    % setup experiment
    rho = rhoVec(i); 
    for j = 1:numReps
        data = cage_trial_split(multiRelease,rho,MALE_CONV_RATE,FEMALE_CONV_RATE,...
            fitnessCostVec,RELATIVE_FECUNDITY,graphBool);
        extinctGenMat(j,i,2) = data.extinctGens;
    end
end

% plot driveType 2 results
figure
plot(rhoVec, mean(extinctGenMat(:,:,2)),'-o','LineWidth',2,'color','black',...
    'MarkerFaceColor','black');
% title({'avg. generation of extinction for', 'recessive female sterile drive and',...
%     'multiple releases'},...
%     'fontweight','bold')
xlim([0,max(rhoVec)+0.05]);
ylim([0,max(extinctGenMat(:,:,2),[],'all')]);
ylabel('avg. generation of extinction','interpreter','latex');
xlabel('release ratio','interpreter','latex');
set(gca,'FontSize',16);
set(gca,'Layer','top')
f = gcf;
exportgraphics(f,'recessive_batch_runs.pdf','Resolution',600);
% hold on

%% Batch of runs for drive 1 using a single set of conditions
multiRelease = true;
MALE_CONV_RATE = 0.9503;
FEMALE_CONV_RATE = NaN;
fitnessCostVec = [0.05, 0];
RELATIVE_FECUNDITY = 0.6274;
graphBool = false; 
rho = 0.2; % 1:4
alpha = 0.75;

% no. of runs to plot (add 1 for the initial plot)
numRuns = 30;
numRunsMat = zeros(numRuns,20);

data = cage_trial_split(multiRelease,rho,MALE_CONV_RATE,FEMALE_CONV_RATE,...
    fitnessCostVec,RELATIVE_FECUNDITY,graphBool);
numRunsMat(1,1:length(data.popVec)) = data.popVec;

figure
plot(0:(length(data.popVec)-1), data.popVec,'LineWidth',2,'color',[0,0,0]+alpha);
hold on

for i = 1:(numRuns-1)
    data = cage_trial_split(multiRelease,rho,MALE_CONV_RATE,FEMALE_CONV_RATE,...
        fitnessCostVec,RELATIVE_FECUNDITY,graphBool);
    
    numRunsMat(i+1,1:length(data.popVec)) = data.popVec;
    plot(0:(length(data.popVec)-1), data.popVec,'LineWidth',2,'color',[0,0,0]+alpha);
end

% now plot the mean
meanPopVec = mean(numRunsMat);

% plot driveType 1 results
plot(0:(length(meanPopVec)-1), meanPopVec,'-k','LineWidth',2);
% title({'batch of runs for dominant female','sterile drive (1:4 release)'},...
%     'fontweight','bold')
xlim([0,12]);
ylim([0,600]);
ylabel('total population','interpreter','latex');
xlabel('generation','interpreter','latex');
set(gca,'FontSize',16);
set(gca,'Layer','top')
f = gcf;
exportgraphics(f,'dominant_single_run.pdf','Resolution',600);

%% Batch of runs for drive 2 using a single set of conditions
multiRelease = true;
MALE_CONV_RATE = 0.9736;
FEMALE_CONV_RATE = 0.9388;
fitnessCostVec = [0.05, 0];
RELATIVE_FECUNDITY = 0.6274;
graphBool = false; 
rho = 0.2; % 1:4
alpha = 0.75;

% no. of runs to plot (add 1 for the initial plot)
numRuns = 30;
numRunsMat = zeros(numRuns,15);

data = cage_trial_split(multiRelease,rho,MALE_CONV_RATE,FEMALE_CONV_RATE,...
    fitnessCostVec,RELATIVE_FECUNDITY,graphBool);
numRunsMat(1,1:length(data.popVec)) = data.popVec;

figure
plot(0:(length(data.popVec)-1), data.popVec,'LineWidth',2,'color',[0,0,0]+alpha);
hold on

for i = 1:(numRuns-1)
    data = cage_trial_split(multiRelease,rho,MALE_CONV_RATE,FEMALE_CONV_RATE,...
        fitnessCostVec,RELATIVE_FECUNDITY,graphBool);
    
    numRunsMat(i+1,1:length(data.popVec)) = data.popVec;
    plot(0:(length(data.popVec)-1), data.popVec,'LineWidth',2,'color',[0,0,0]+alpha);
end

% now plot the mean
meanPopVec = mean(numRunsMat);

% plot driveType 1 results
plot(0:(length(meanPopVec)-1), meanPopVec,'-k','LineWidth',2);
% title({'batch of runs for recessive female','sterile drive (1:4 release)'},...
%     'fontweight','bold')
xlim([0,12]);
ylim([0,600]);
ylabel('total population','interpreter','latex');
xlabel('generation','interpreter','latex');
set(gca,'FontSize',16);
set(gca,'Layer','top')
f = gcf;
exportgraphics(f,'recessive_single_run.pdf','Resolution',600);

%% Bar plot for frequency of extinctions resulting from single release for
%  both drive types

rhoVec = 1:10;
% no. of replicates per set release
numReps = 100;
% frequency of extinction within 20 generations
extinctFreqMat = nan(numReps,length(rhoVec),2);

% parameters for simulation
multiRelease = false;
MALE_CONV_RATE = 0.9503;
FEMALE_CONV_RATE = NaN;
fitnessCostVec = [0.05, 0];
RELATIVE_FECUNDITY = 0.6274;
graphBool = false; 

for i = rhoVec
    sprintf("Running sim %.0f of %.0f",i,length(rhoVec))
    
    rho = i; 
    for j = 1:numReps
        data = cage_trial_split(multiRelease,rho,MALE_CONV_RATE,FEMALE_CONV_RATE,...
            fitnessCostVec,RELATIVE_FECUNDITY,graphBool); 
        extinctFreqMat(j,i,1) = ~isnan(data.extinctGens);
    end
end

% store frequencies for drive type 1
extinctFreqVec_1 = mean(extinctFreqMat(:,:,1));

% parameters for simulation
MALE_CONV_RATE = 0.9736;
FEMALE_CONV_RATE = 0.9388;

for i = rhoVec
    sprintf("Running sim %.0f of %.0f (drive type 2)",i,length(rhoVec))
    
    rho = i; 
    for j = 1:numReps
        data = cage_trial_split(multiRelease,rho,MALE_CONV_RATE,FEMALE_CONV_RATE,...
            fitnessCostVec,RELATIVE_FECUNDITY,graphBool); 
        extinctFreqMat(j,i,2) = ~isnan(data.extinctGens);
    end
end

% store frequencies for drive type 2
extinctFreqVec_2 = mean(extinctFreqMat(:,:,2));

figure
bar([extinctFreqVec_1; extinctFreqVec_2]','BarWidth',1);
ylim([0,1]);
xlabel('release ratio','interpreter','latex');
ylabel('extinction probability','interpreter','latex');
set(gca,'FontSize',16);
set(gca,'Layer','top')
f = gcf;
exportgraphics(f,'single_release_bar_plot.pdf','Resolution',600);

%% Bar plot for frequency of extinctions resulting from single release for
%  both drive types, assuming an autonomous drive

rhoVec = 0.01:0.01:0.1;
% no. of replicates per set release
numReps = 100;
% frequency of extinction within 20 generations
extinctFreqMat = nan(numReps,length(rhoVec),2);

% parameters for simulation
multiRelease = false;
MALE_CONV_RATE = 0.9503;
FEMALE_CONV_RATE = NaN;
fitnessCostVec = [0.05, 0];
RELATIVE_FECUNDITY = 0.6274;
graphBool = false; 

for i = 1:length(rhoVec)
    sprintf("Running sim %.0f of %.0f",i,length(rhoVec))
    
    rho = rhoVec(i); 
    for j = 1:numReps
        data = cage_trial_full(multiRelease,rho,MALE_CONV_RATE,FEMALE_CONV_RATE,...
            fitnessCostVec,RELATIVE_FECUNDITY,graphBool); 
        extinctFreqMat(j,i,1) = ~isnan(data.extinctGens);
    end
end

% store frequencies for drive type 1
extinctFreqVec_1 = mean(extinctFreqMat(:,:,1));

% parameters for simulation
MALE_CONV_RATE = 0.9736;
FEMALE_CONV_RATE = 0.9388;

for i = 1:length(rhoVec)
    sprintf("Running sim %.0f of %.0f (drive type 2)",i,length(rhoVec))
    
    rho = rhoVec(i); 
    for j = 1:numReps
        data = cage_trial_full(multiRelease,rho,MALE_CONV_RATE,FEMALE_CONV_RATE,...
            fitnessCostVec,RELATIVE_FECUNDITY,graphBool); 
        extinctFreqMat(j,i,2) = ~isnan(data.extinctGens);
    end
end

% store frequencies for drive type 2
extinctFreqVec_2 = mean(extinctFreqMat(:,:,2));

figure
bar([extinctFreqVec_1; extinctFreqVec_2]','BarWidth',1);
xticks(1:1:length(rhoVec));
xticklabels(0.01:0.01:0.1); 
ylim([0,1]);
xlabel('release ratio','interpreter','latex');
ylabel('extinction probability','interpreter','latex');
% legend('dominant female sterile','recessive female sterile','location','best');
set(gca,'FontSize',16);
set(gca,'Layer','top')
f = gcf;
exportgraphics(f,'single_release_bar_plot.pdf','Resolution',600);