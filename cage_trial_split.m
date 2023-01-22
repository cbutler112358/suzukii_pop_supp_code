% Cage trial simulator (last updated 01/22/2023)
% Author: Cole Butler 
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Function that performs a single small cage trial experiment. Separate 
% Excel files are prepared for determination of zygote frequencies for each
% pairing. The program (and corresponding Excel sheets) assumes that the 
% Cas9 is X-linked.
%
% If a single release is used (multiRelease == false), the simulation 
% automatically terminates if the gRNA and Cas9 components go extinct.
%
% INPUTS:
%   multiRelease -- Boolean variable, true if multiple releases occur.
%   rho -- release ratio
%   MALE_CONV_RATE -- conversion efficiency in males
%   FEMALE_CONV_RATE -- conversion efficiency in females; should be left as
%       NaN in the dominant female fertile case
%   fitnessCostVec -- vector containing fitness costs of construct; (1)
%       Cas9 and (2) gRNA
%   RELATIVE_FECUNDIY -- relative fecundity of hemizygous females
%   graphBool -- Boolean variable; do you want to see plots of the results?
%
% OUTPUTS:
% A single structual array containing the following variables:
%   extinctGens -- generation of extinction of cage pop.
%   popVec -- total no. of flies per generation
%   femaleVec -- total no. of females per generation
%   femaleMat -- matrix showing females per generation by genotype
%   maleMat -- as above but for males
%   gRNA_alleleFreqVec -- frequency of the gRNA per generation in the total
%       population (including released pop.)
%   zygoteFreqsMales -- male breeding table produced by genBT algorithm
%   zygoteFreqsFemales -- as above but for females

function [dataMat] = cage_trial_split(multiRelease,rho,...
    MALE_CONV_RATE,FEMALE_CONV_RATE,fitnessCostVec,RELATIVE_FECUNDITY,graphBool)

    if isnan(FEMALE_CONV_RATE)
        % dominant female sterile
        releaseInd = 5; % ZYAa males released
        % construct zygote frequencies
        breedingTables = genBT_simple_tethered_homing(MALE_CONV_RATE, ...
            FEMALE_CONV_RATE, RELATIVE_FECUNDITY, 1);
        depoBool = false; % boolean flag for maternal Cas9 deposition
    else
        % recessive female sterile
        releaseInd = 5; % ZYAa males released
        breedingTables = genBT_simple_tethered_homing(MALE_CONV_RATE, ...
            FEMALE_CONV_RATE, RELATIVE_FECUNDITY, 2);
        depoBool = true; 
    end
    
    zygoteFreqsMales = breedingTables.male_genoArray;
    zygoteFreqsFemales = breedingTables.female_genoArray;
    
    % fitness costs of respective split drive components
    CAS9_COST = fitnessCostVec(1);
    GRNA_COST = fitnessCostVec(2); 
    
    %% run the sim
    % store no. of genotypes
    NUM_GENOTYPES_MALES = 6;
    NUM_GENOTYPES_FEMALES = 9;

    maleMat = zeros(1, NUM_GENOTYPES_MALES);
    % matrix for FERTILE females
    femaleMat = zeros(1, NUM_GENOTYPES_FEMALES);
    % matrix for females made sterile by Cas9 deposition
    sterileFemaleMat = zeros(1, NUM_GENOTYPES_FEMALES);
    % matrix only for plotting (omits released individuals each gen)
    %%% plotMat = NaN(1,NUM_GENS+1);
    % allele frequency vectors
    gRNA_alleleFreqVec = NaN(1,1); 
    CAS9_alleleFreqVec = NaN(1,1); 
    fitnessVec_males = [(1-GRNA_COST)^2, 1-GRNA_COST, 1,... 
        (1-CAS9_COST)*(1-GRNA_COST)^2, (1-CAS9_COST)*(1-GRNA_COST), 1-CAS9_COST];
    fitnessVec_females = [(1-GRNA_COST)^2, 1-GRNA_COST, 1, ...
        (1-CAS9_COST)*(1-GRNA_COST)^2, (1-CAS9_COST)*(1-GRNA_COST), 1-CAS9_COST, ...
        ((1-CAS9_COST)^2)*(1-GRNA_COST)^2, ((1-CAS9_COST)^2)*(1-GRNA_COST), (1-CAS9_COST)^2];

    LAMBDA = 50;

    % 0 generation, begin with 300 males and females 
    INIT_POP = 300;
    RELEASE_RATIO = rho; 
    maleMat(1,3) = INIT_POP;
    femaleMat(1,3) = INIT_POP;
    % excludes released males... 
    %%% plotMat(1,:) = femaleMat(1,:) + maleMat(1,:);
    popVec = zeros(1, 1); 
    popVec(1) = sum(maleMat(1,:)) + sum(femaleMat(1,:));
    % add ZYAa males to test bottle
    maleMat(1,releaseInd) = round(RELEASE_RATIO*INIT_POP);
    
    % total population 
    %%% totalMat = maleMat + femaleMat;
    gRNA_alleleFreqVec(1) = (2*maleMat(1,1) + maleMat(1,2) + 2*maleMat(1,4) + maleMat(1,5) + ...
        2*femaleMat(1,1) + femaleMat(1,2) + 2*femaleMat(1,4) + femaleMat(1,5) + 2*femaleMat(1,7) + ...
        femaleMat(1,8)) / (2*(sum(maleMat(1,:)) + sum(femaleMat(1,:)))); 
    CAS9_alleleFreqVec(1) = (maleMat(1,4) + maleMat(1,5) + maleMat(1,6) + ...
        2*femaleMat(1,7) + 2*femaleMat(1,8) + 2*femaleMat(1,9) + ...
        femaleMat(1,4) + femaleMat(1,5) + femaleMat(1,6)) / (sum(maleMat(1,:)) + 2*sum(femaleMat(1,:))); 
    femaleVec = zeros(1, 1); 
    femaleVec(1) = sum(femaleMat(1,:)); 
    
    % expected no. of female eggs in control bottle (doesn't change)
    numFemales_control = (1/2)*LAMBDA*INIT_POP; 

    % prop. of surviving offspring, regardless of genotype
    % NOTE: the following lines are only used for calculating expected value of
    % female offspring produced

    % apply early-acting fitness cost
    offspringZygoteFreqs = zygoteFreqsFemales(:,1:NUM_GENOTYPES_FEMALES).*fitnessVec_females;
    fitnessMat = zeros(9,6); 
    sumZygoteFitness = sum(offspringZygoteFreqs,2);
    % male genotype by column, female genotype by row; 
    % probability offspring of certain genotype survive to adulthood in
    % next generation
    for j = 1:(54/9)
        fitnessMat(:,j) = sumZygoteFitness(((j-1)*9+1):(j*9));
    end

    % main for loop
    extinctFlag = true; 
    i = 2; % counter
    while (extinctFlag)
        
        % fprintf("Simulating generation %.f of %.f\n",i-1,NUM_GENS);
        % probability of mating with male of a certain genotype
        maleProbVec = maleMat(i-1,:)/sum(maleMat(i-1,:));


        %% expected number of female offspring from test bottle
        % (recall that femaleMat excludes sterile females)
        numFemales_test = sum(sum((maleProbVec'*femaleMat(i-1,:)*(LAMBDA/2)).*(fitnessMat')));

        %% test bottle

        % matings of females with males; female genotype per column, male
        % genotype per row
        breedFemaleMat = mnrnd(femaleMat(i-1,:)', maleProbVec)';
        % no. of eggs produced from pairings
        numOffspringMat = poissrnd(breedFemaleMat*LAMBDA);
        % number of eggs by zygote stored in numOffspringVec
        % (organized by male genotype first)
        pairingVec = reshape(numOffspringMat',1,[]);
        % determine males and females
        numMales = binornd(pairingVec,0.5);
        numFemales = pairingVec - numMales;
        % male genotypes straightforward
        maleGenotypeVec = sum(mnrnd(numMales',zygoteFreqsMales));
        
        % disp(pairingVec);
        
        if (depoBool)
            %%% some female progeny are sterile because of maternal deposition
            % of Cas9---keep track of which
            
            % row --> motherDepoInd (crossings affected)
            % col --> gRNA          (inherit gRNA)
            femaleGenos = mnrnd(numFemales',zygoteFreqsFemales);
            % filter by crossings affected by maternal deposition... all
            % female progeny resulting from such crossings are sterile,
            % regardless of whether they actually inherited a gRNA
            sterileFemales = ((femaleGenos').*breedingTables.motherDepoInd)';
            % remainder are fertile
            fertileFemales = femaleGenos - sterileFemales; 
            sterileFemaleGenotypeVec = sum(sterileFemales);
            fertileFemaleGenotypeVec = sum(fertileFemales);
            
            femaleGenotypeVec = sterileFemaleGenotypeVec + fertileFemaleGenotypeVec;

            % remove eggs that are not produced
            maleGenotypeVec = maleGenotypeVec(1:NUM_GENOTYPES_MALES);
            femaleGenotypeVec = femaleGenotypeVec(1:NUM_GENOTYPES_FEMALES);
            sterileFemaleGenotypeVec = sterileFemaleGenotypeVec(1:NUM_GENOTYPES_FEMALES); 
            %%% unused
            %  fertileFemaleGenotypeVec = fertileFemaleGenotypeVec(1:NUM_GENOTYPES_FEMALES);
        else
            femaleGenotypeVec = sum(mnrnd(numFemales',zygoteFreqsFemales));
            maleGenotypeVec = maleGenotypeVec(1:NUM_GENOTYPES_MALES);
            femaleGenotypeVec = femaleGenotypeVec(1:NUM_GENOTYPES_FEMALES);
        end

        % determine genotypes and sexes of next generation...
        N = floor(INIT_POP*(numFemales_test/numFemales_control));
        % disp([i,INIT_POP*(numFemales_test/numFemales_control)])
        % disp([i,N]); 
        % disp(maleGenotypeVec); 
        
        if (sum(maleGenotypeVec) + sum(femaleGenotypeVec)) < N
            % no. needed to initialize next generation is GREATER than the
            % number of eggs actually being produced, take the smaller
            % number
            
            % at times, this occurs when all females mate with sterile
            % males, so that the no. in the next generation is 0 and thus
            % less than any positive N (which is controlled by the EXPECTED
            % no. of females produced)
            warning('No. of eggs selected larger than no. produced.');
            N = sum(maleGenotypeVec) + sum(femaleGenotypeVec); 
        end
        
        % genotype probability by sex
        maleGenotypeProb = maleGenotypeVec/sum(maleGenotypeVec);
        femaleGenotypeProb = femaleGenotypeVec/sum(femaleGenotypeVec);
        if (depoBool)
            % probability of selecting sterile females
            sterileFemaleGenotypeProb = sterileFemaleGenotypeVec./femaleGenotypeVec;
            sterileFemaleGenotypeProb(isnan(sterileFemaleGenotypeProb)) = 0;
        end

        % next generation
        if (N == 0) || (all(maleGenotypeVec == 0)) || (all(femaleGenotypeVec == 0))
            % extinction
            maleMat(i,:) = 0; 
            femaleMat(i,:) = 0;
            sterileFemaleMat(i,:) = 0;
        else 
            % persistence
            maleMat(i,:) = mnrnd(N, maleGenotypeProb);
            if (depoBool)
                % select eggs as usual
                femaleSelect = mnrnd(N, femaleGenotypeProb);
                % some are sterile, others fertile
                sterileFemaleMat(i,:) = binornd(femaleSelect, sterileFemaleGenotypeProb);
                femaleMat(i,:) = femaleSelect - sterileFemaleMat(i,:);
                % EA fitness cost
                sterileFemaleMat(i,:) = binornd(sterileFemaleMat(i,:), fitnessVec_females); 
            else                            
                femaleMat(i,:) = mnrnd(N, femaleGenotypeProb);
            end
            
            % early-acting fitness cost imposed on eggs; survival to adulthood
            % assumed constant between genotypes
            maleMat(i,:) = binornd(maleMat(i,:), fitnessVec_males);
            femaleMat(i,:) = binornd(femaleMat(i,:), fitnessVec_females);
        end
        
        % update plot matrix and vector, omitting released flies 
        %%% plotMat(i,:) = maleMat(i,:) + femaleMat(i,:);
        % disp(femaleMat(i,:));
        
        if (depoBool)
            popVec(i) = sum(maleMat(i,:)) + sum(femaleMat(i,:)) + sum(sterileFemaleMat(i,:));

            if (multiRelease)
                % release adults if performing multiple releases
                maleMat(i,releaseInd) = maleMat(i,releaseInd) + round(RELEASE_RATIO*INIT_POP);  
            end
            
            %%% totalMat(i,:) = maleMat(i,:) + femaleMat(i,:);
            gRNA_alleleFreqVec(i) = (2*maleMat(i,1) + maleMat(i,2) + 2*maleMat(i,4) + maleMat(i,5) + ...
                2*femaleMat(i,1) + femaleMat(i,2) + 2*femaleMat(i,4) + femaleMat(i,5) + 2*femaleMat(i,7) + ...
                femaleMat(i,8) + 2*sterileFemaleMat(i,1) + sterileFemaleMat(i,2) + 2*sterileFemaleMat(i,4) + ...
                sterileFemaleMat(i,5) + 2*sterileFemaleMat(i,7) + sterileFemaleMat(i,8)) / ...
                (2*(sum(maleMat(i,:)) + sum(femaleMat(i,:)) + sum(sterileFemaleMat(i,:)))); 
            
            CAS9_alleleFreqVec(i) = (maleMat(i,4) + maleMat(i,5) + maleMat(i,6) + ...
                2*femaleMat(i,7) + 2*femaleMat(i,8) + 2*femaleMat(i,9) + ...
                femaleMat(i,4) + femaleMat(i,5) + femaleMat(i,6) + ...
                2*sterileFemaleMat(i,7) + 2*sterileFemaleMat(i,8) + 2*sterileFemaleMat(i,9) + ...
                sterileFemaleMat(i,4) + sterileFemaleMat(i,5) + sterileFemaleMat(i,6)) / (sum(maleMat(i,:)) + 2*sum(femaleMat(i,:)) + 2*sum(sterileFemaleMat(i,:))); 
            
            femaleVec(i) = sum(femaleMat(i,:) + sterileFemaleMat(i,:)); 
        else
            popVec(i) = sum(maleMat(i,:)) + sum(femaleMat(i,:));
            
            if (multiRelease)
                % release adults if performing multiple releases
                maleMat(i,releaseInd) = maleMat(i,releaseInd) + round(RELEASE_RATIO*INIT_POP);  
            end            
            
            %%% totalMat(i,:) = maleMat(i,:) + femaleMat(i,:);
            gRNA_alleleFreqVec(i) = (2*maleMat(i,1) + maleMat(i,2) + 2*maleMat(i,4) + maleMat(i,5) + ...
                2*femaleMat(i,1) + femaleMat(i,2) + 2*femaleMat(i,4) + femaleMat(i,5) + 2*femaleMat(i,7) + ...
                femaleMat(i,8)) / (2*(sum(maleMat(i,:)) + sum(femaleMat(i,:)))); 
            
            CAS9_alleleFreqVec(i) = (maleMat(i,4) + maleMat(i,5) + maleMat(i,6) + ...
                2*femaleMat(i,7) + 2*femaleMat(i,8) + 2*femaleMat(i,9) + ...
                femaleMat(i,4) + femaleMat(i,5) + femaleMat(i,6)) / (sum(maleMat(i,:)) + 2*sum(femaleMat(i,:))); 
            
            femaleVec(i) = sum(femaleMat(i,:)); 
        end

        % extinction occurs when there are no females (not counting made 
        % sterile by maternal deposition) or males left
        femaleBool = (sum(femaleMat(i,:))  > 0);
        % only count males that are NOT released
        maleBool = (sum(maleMat(i,[1:4,6])) > 0);
        extinctFlag = (maleBool & femaleBool);
        
        if (~extinctFlag)
            % fprintf("Population extinct at generation %.f\n",i-1);
            extinctGens = i-1;
            break
        end
        
        % for single releases, drive failure occurs when (i) the population 
        % is not extinct and (ii) drive components have disappeared
        if ((multiRelease == false) && extinctFlag && (gRNA_alleleFreqVec(i) == 0) && (CAS9_alleleFreqVec(i) == 0))
            extinctGens = nan;
            graphBool = false; 
            break
        end
        i = i + 1;
    end

    % replace nans in gRNA and CAS9 vectors
    gRNA_alleleFreqVec(isnan(gRNA_alleleFreqVec)) = 0;
    CAS9_alleleFreqVec(isnan(CAS9_alleleFreqVec)) = 0;

    if (graphBool)
        % plot everything
        close all 
        figure
        subplot(3,1,1)
        plot(0:extinctGens, popVec,'-o','LineWidth',2,'color','blue',...
            'MarkerFaceColor','blue');
        ylim([0,max(popVec)+10]);
        xlim([0,extinctGens]);
        ylabel('total pop.','interpreter','latex');
        xlabel('generation','interpreter','latex');
        set(gca,'FontSize',16);

        subplot(3,1,2)
        plot(0:extinctGens, femaleVec,'-o','LineWidth',2,'color','red',...
            'MarkerFaceColor','red');
        ylim([0,max(femaleVec)+10]);
        xlim([0,extinctGens]);
        ylabel('no. of females','interpreter','latex');
        xlabel('generation','interpreter','latex');
        set(gca,'FontSize',16);
        
        subplot(3,1,3)
        plot(0:extinctGens, gRNA_alleleFreqVec,'-o','LineWidth',2,'color','black',...
            'MarkerFaceColor','black');
        ylim([0,1]);
        xlim([0,extinctGens]);
        ylabel('allele freq.','interpreter','latex');
        xlabel('generation','interpreter','latex');
        set(gca,'FontSize',16);
    end
    
    % array containing all data of interest
    dataMat = struct();
    dataMat.extinctGens = extinctGens;
    dataMat.popVec = popVec;
    dataMat.femaleVec = femaleVec;
    dataMat.femaleMat = femaleMat;
    dataMat.sterileFemaleMat = sterileFemaleMat;
    dataMat.maleMat = maleMat;
    dataMat.gRNA_alleleFreqVec = gRNA_alleleFreqVec;
    dataMat.CAS9_alleleFreqVec = CAS9_alleleFreqVec;
    dataMat.zygoteFreqsMales = zygoteFreqsMales;
    dataMat.zygoteFreqsFemales = zygoteFreqsFemales;
    
    return
end





