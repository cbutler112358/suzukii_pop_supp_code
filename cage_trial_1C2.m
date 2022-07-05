% Cage trial simulator (last updated 06/08/2022)
% Author: Cole Butler 
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Function that performs a single small cage trial experiment. Separate 
% Excel files are prepared containing the zygote frequencies of each
% pairing, assuming females have a reduced fecundity in the recessive
% female sterile drive case. 
%
% The following code snippet prepares the "raw" input for the function 
% below, taking Drive 1 (dominant female lethal) as an example:
% 
%     filename = "drive1_fsRIDL.xlsx";
%     if contains(filename,"1")
%         MALE_CONV_RATE = 0.1382;
%         FEMALE_CONV_RATE = NaN;
%         releaseInd = 2; % Aa males released, dominant female sterile
%     elseif contains(filename,"2")
%         MALE_CONV_RATE = 0.4136;
%         FEMALE_CONV_RATE = 0.2262;
%         releaseInd = 1; % AA males released, recessive female sterile
%     else
%         error("Incorrect filename.");
%     end
%     raw = readcell(filename);
% 
% Remaining inputs and outputs are discussed below:
% NUM_GENS -- no. of generations within which we are looking for pop.
%             extinction 
% multiRelease -- Boolean value, "true" for multiple releases, "false"
%                 otherwise
% rho -- release ratio of transgenic males released to (initial) wild-type
%        males, e.g. 0.2 is 1 transgenic male released for every 5
%        wild-type males
% raw -- cell array depending on Excel files, see example above
% MALE_CONV_RATE -- conversion efficiency in males (denoted 'g' in Excel
%                   file)
% FEMALE_CONV_RATE -- conversion efficiency in females (denoted 'h' in
%                     Excel file); if left as NaN, conversion does not
%                     occur in females (dominant female sterile)
% driveInd -- index denoting which drive is being simulated: 1 is dominant 
%             female sterile, Drive 2 is recessive female sterile
% FITNESS_COST -- relative fitness cost: 1 for wild-type, 1-s for
%                 hemizygotes, (1-s)^2 for homozygotes; affects hatching
%                 probability
% graphBool -- Boolean value; do you want the results to be plotted? 
%
% Outputs:
% extinctGens -- generation of extinction, NaN if extinction was not
%                achieved
% totalMat -- matrix of total population, where genotype is by column
% femaleVec -- vector of female population
% popVec -- vector of total population, excluding releases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [extinctGens, totalMat, femaleVec, popVec] = cage_trial_1C2(NUM_GENS,multiRelease,rho,...
    raw,MALE_CONV_RATE,FEMALE_CONV_RATE,driveInd,FITNESS_COST, graphBool)

    if (driveInd==1)
        releaseInd = 2; % Aa males released, dominant female sterile
    elseif (driveInd==2)
        releaseInd = 1; % AA males released, recessive female sterile
    else
        error("Incorrect filename.");
    end

    % FITNESS_COST = 0.05; % early-acting? 

    % raw = readcell(filename);
    genoFreqs = raw(2:end,3:end); 
    % convert to string
    test = string(genoFreqs(:,:));
    % male genotype conversion rate
    genoFreqs = replace(test,'g',sprintf('%.2f',MALE_CONV_RATE));
    % female genotype conversion rate 
    if ~isnan(FEMALE_CONV_RATE)
        genoFreqs = replace(genoFreqs,'h',sprintf('%.2f',FEMALE_CONV_RATE));
    end
    genoFreqs = replace(genoFreqs,'s',sprintf('%.2f',FITNESS_COST));
    % construct new array 
    zygoteFreqs = zeros(size(genoFreqs));
    matDim = size(genoFreqs);

    for i = 1:matDim(1)
        for j = 1:matDim(2)
            zygoteFreqs(i,j) = eval(genoFreqs(i,j));
        end
    end

    % add a column for eggs that are not produced; zygoteFreqs follows male
    % first pairing specification
    deathProbVec = ones(length(zygoteFreqs),1) - sum(zygoteFreqs,2);
    zygoteFreqs = [zygoteFreqs, deathProbVec];

    %% run the sim
    % store no. of genotypes
    NUM_GENOTYPES = 3;

    maleMat = zeros(NUM_GENS+1, NUM_GENOTYPES);
    femaleMat = zeros(NUM_GENS+1, NUM_GENOTYPES);
    % matrix only for plotting (omits released individuals each gen)
    plotMat = zeros(NUM_GENS+1, NUM_GENOTYPES);
    % allele frequency vector
    alleleFreqVec = NaN(1,NUM_GENS+1); 
    fitnessVec = [(1-FITNESS_COST)^2, 1-FITNESS_COST, 1];

    LAMBDA = 50;

    % 0 generation, begin with 30 males and females 
    INIT_POP = 300;
    RELEASE_RATIO = rho; 
    maleMat(1,3) = INIT_POP;
    femaleMat(1,3) = INIT_POP;
    % excludes released males... 
    plotMat(1,:) = femaleMat(1,:) + maleMat(1,:);
    popVec = zeros(1, NUM_GENS+1); 
    popVec(1) = sum(maleMat(1,:) + femaleMat(1,:));
    % add ZYAA males to test bottle
    maleMat(1,releaseInd) = round(RELEASE_RATIO*INIT_POP);
    
    % total population 
    totalMat = maleMat + femaleMat;
    alleleFreqVec(1) = (2*totalMat(1,1) + totalMat(1,2))/(2*sum(totalMat(1,:)));
    femaleVec = zeros(1, NUM_GENS+1); 
    femaleVec(1) = sum(femaleMat(1,:)); 

    % expected no. of female eggs in control bottle (doesn't change)
    numFemales_control = (1/2)*LAMBDA*INIT_POP; 

    % prop. of surviving offspring, regardless of genotype
    % NOTE: lines 72-89 are only used for calculating expected value of
    % female offspring produced
    offspringZygoteFreqs = zygoteFreqs(:,1:3); 
    % split up by sex...
    female_offspringZygoteFreqs = (1/2)*offspringZygoteFreqs;
    male_offspringZygoteFreqs = (1/2)*offspringZygoteFreqs;
    % apply (early-acting) fitness costs, by genotype
    female_offspringZygoteFreqs = female_offspringZygoteFreqs.*fitnessVec;
    male_offspringZygoteFreqs = male_offspringZygoteFreqs.*fitnessVec;
    % back to prop. of offspring
    offspringZygoteFreqs = male_offspringZygoteFreqs + female_offspringZygoteFreqs; 
    
    fitnessMat = zeros(3,3); 
    sumZygoteFitness = sum(offspringZygoteFreqs,2);
    % male genotype by column, female genotype by row; 
    % probability offspring of certain genotype survive to adulthood in
    % next generation
    fitnessMat(:,1) = sumZygoteFitness(1:3);
    fitnessMat(:,2) = sumZygoteFitness(4:6);
    fitnessMat(:,3) = sumZygoteFitness(7:9);

    % main for loop
    for i = 2:(NUM_GENS+1)
        
        % fprintf("Simulating generation %.f of %.f\n",i-1,NUM_GENS);
        % probability of mating with male of a certain genotype
        maleProbVec = maleMat(i-1,:)/sum(maleMat(i-1,:));


        %% expected number of female offspring from test bottle
        numFemales_test = (1/2)*sum(sum((maleProbVec'*femaleMat(i-1,:)*LAMBDA).*(fitnessMat')));


        %% test bottle

        % matings of females with males; female genotype per column, male
        % genotype per row
        breedFemaleMat = mnrnd(femaleMat(i-1,:)', maleProbVec)';
        % no. of eggs produced from pairings
        numOffspringMat = poissrnd(breedFemaleMat*LAMBDA);
        % number of eggs by zygote stored in numOffspringVec
        % (organized by male genotype first)
        pairingVec = reshape(numOffspringMat',1,[]);
        offspringGenotypeVec = sum(mnrnd(pairingVec',zygoteFreqs));
        % remove eggs that are not produced
        offspringGenotypeVec = offspringGenotypeVec(1:3); 
        % determine no. of females and males produced from mating
        numMales = binornd(offspringGenotypeVec,0.5);
        numFemales = offspringGenotypeVec - numMales; 

        % determine genotypes and sexes of next generation...
        N = round(INIT_POP*(numFemales_test/numFemales_control));
        
        if sum(numMales + numFemales) < N
            % no. needed to initialize next generation is GREATER than the
            % number of eggs actually being produced, take the smaller
            % number...
            N = sum(numMales + numFemales); 
        end
        
        % genotype probability by sex
        maleGenotypeProb = numMales/sum(numMales);
        femaleGenotypeProb = numFemales/sum(numFemales);

        % next generation
        if (N == 0) || (all(numMales == 0)) || (all(numFemales == 0))
            % extinction
            maleMat(i,:) = 0; 
            femaleMat(i,:) = 0;
        else 
            % persistence
            maleMat(i,:) = mnrnd(N, maleGenotypeProb);
            femaleMat(i,:) = mnrnd(N, femaleGenotypeProb);
            % early-acting fitness cost imposed on eggs; survival to adulthood
            % assumed constant between genotypes
            maleMat(i,:) = binornd(maleMat(i,:), fitnessVec);
            femaleMat(i,:) = binornd(femaleMat(i,:), fitnessVec);
        end
        
        % update plot matrix and vector, omitting released flies 
        plotMat(i,:) = maleMat(i,:) + femaleMat(i,:);
        popVec(i) = sum(maleMat(i,:) + femaleMat(i,:));

        if (multiRelease)
            % release adults if performing multiple releases
            maleMat(i,releaseInd) = maleMat(i,releaseInd) + round(RELEASE_RATIO*INIT_POP);  
        end
        
        totalMat(i,:) = maleMat(i,:) + femaleMat(i,:);
        alleleFreqVec(i) = (2*totalMat(i,1) + totalMat(i,2))/(2*sum(totalMat(i,:)));
        femaleVec(i) = sum(femaleMat(i,:)); 
        
        % extinction occurs when either all females or all males disappear 
        femaleBool = (sum(femaleMat(i,:))  > 0);
        maleBool = (sum(maleMat(i,:)) > 0);
        extinctFlag = (femaleBool & maleBool);
        if ~extinctFlag
            % fprintf("Population extinct at generation %.f\n",i-1);
            MAX_GENS = i-1;
            extinctGens = i-1;
            break
        end
    end

    if extinctFlag
        % the population did not go extinct
        MAX_GENS=NUM_GENS;
        extinctGens = NaN;
    end


    if (graphBool)
        % plot everything
        subplot(2,1,1)
        plot(0:NUM_GENS, sum(plotMat(1:(NUM_GENS+1),:),2),'-b','LineWidth',2);
        ylim([0,max(sum(plotMat,2))+10]);
        xlim([0,NUM_GENS]);
        ylabel('total pop.','interpreter','latex');
        xlabel('generation','interpreter','latex');
        set(gca,'FontSize',16);

        subplot(2,1,2)
        plot(0:NUM_GENS, femaleVec,'-r','LineWidth',2);
        ylim([0,max(femaleVec)+10]);
        xlim([0,NUM_GENS]);
        ylabel('no. of females','interpreter','latex');
        xlabel('generation','interpreter','latex');
        set(gca,'FontSize',16);
    end
    
    return
end





