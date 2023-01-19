% Script for running simulations of the (updated) D. suzukii cage trial
% code in cage_trial_split.m. These simulations use the most recent set of
% parameters (01/19/23). 
%
% Also generates figures for showing average gen. of extinction. 

%% Test sims

% 
% multiRelease = true;
% rho = 0.05;
% MALE_CONV_RATE = 0.9736;
% FEMALE_CONV_RATE = 0.9388;
% fitnessCostVec = [0.05, 0];
% RELATIVE_FECUNDITY = 0.6274;
% graphBool = false; 
% 
% % random_num = 581027; % randi([0,1000000])
% % rng(random_num);
% 
% data = cage_trial_split(multiRelease,rho,MALE_CONV_RATE,FEMALE_CONV_RATE,...
%     fitnessCostVec,RELATIVE_FECUNDITY,graphBool);

% subplot(1,2,1)
% plot(0:19,data.femaleMat(:,2) + data.femaleMat(:,5),'-r','linewidth',2)
% xlabel("generation",'interpreter','latex');
% ylabel("female hemizygotes in cage",'interpreter','latex');
% set(gca,'fontsize',16);
% 
% subplot(1,2,2)
% plot(1:19,N,'-r','linewidth',2)
% xlabel("generation",'interpreter','latex');
% ylabel("$N=300\times(T/C)$",'interpreter','latex');
% set(gca,'fontsize',16);

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
