% Function algorithmically producing breeding tables for a simple homing 
% drive using conditions of the D. suzukii study.
% Author: Cole Butler 
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% The following function produces a breeding table given conversion
% efficiencies for either sex. Constructs (transgenes) are indicated by
% capital letters while wild-type are denoted by lower-case letters.
%
% Note that sex determination and fitness costs are incorporated externally
% (i.e. not within this function). This function includes a reduced
% fecundity for hemizygous females, and all homozygotes are sterile.
%
% INPUTS:
% -- MALE_CONVERSION_PROB: Conversion efficiency of construct in males.
% -- FEMALE_CONVERSION_PROB: Conversion efficiency of construct in females.
% -- RELATIVE_FEMALE_FECUND: Fecundity of hemizygous females relative to wt
%    females.
% -- driveType: The type of drive being considered, 1 for dominant sterile 
%   and 2 for recessive sterile
%
% OUTPUTS:
% -- breedingTable: A structure array containing the male and female
%   breeding tables
% -- motherDepoInd: Vector w/indices of pairings affected by maternal
%   deposition of Cas9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [breedingTables] = genBT_simple_homing_suzukii(MALE_CONVERSION_PROB, ...
    FEMALE_CONVERSION_PROB,RELATIVE_FEMALE_FECUND,driveType)

%% produce breeding table for zygotes
% tmpArray = readcell('SH_breeding_table.xlsx');
tmpArray = readcell('SH_breeding_table.xlsx');

% pull out the part of the breeding table we will update
genoArray = tmpArray(2:end,3:end);
arrayDim = size(genoArray);
numPairings = arrayDim(1);
numGenotypes = arrayDim(2);
% include a column for nonviable progeny
genoArray = zeros(numPairings, numGenotypes); 
% vector storing indices of pairings producing female offspring that are
% affected by maternal deposition
motherDepoInd = zeros(1,numPairings); 

%% dominant sterile--hemizygous males fertile, homozygotes and 
% hemizygous females sterile
if (driveType == 1)
    % for each possible pairing, determine zygote frequencies...
    for i = 1:numPairings
        % initialize haplotypes @ 0 (for male progeny, the father never
        % contributes a X chromosome)
        fatherHaplo = struct();    
        fatherHaplo.A = 0;
        fatherHaplo.a = 0;

        motherHaplo = fatherHaplo;    

        dadGenotype = string(tmpArray(1+i,1));
        momGenotype = string(tmpArray(1+i,2));

        % first make sure that progeny are actually produced
        nonviablePairing = (dadGenotype == "AA" || momGenotype == "AA" || ...
            momGenotype == "Aa");
        
        if (~nonviablePairing)
            % progeny are produced, check for homing in either parent
            homingDadBool = (dadGenotype == "Aa");

            if (homingDadBool)
                % with MALE_CONVERSION_PROB, then father's genotype is AA
                fatherHaplo.A = (1-MALE_CONVERSION_PROB)*0.5 + MALE_CONVERSION_PROB;
                fatherHaplo.a = (1-MALE_CONVERSION_PROB)*0.5;
            else 
                % homing is not possible
                fatherHaplo.A = fatherHaplo.A + count(dadGenotype,"A")/2;
                fatherHaplo.a = fatherHaplo.a + count(dadGenotype,"a")/2;
            end
            
            % homing is not possible in females
            motherHaplo.A = motherHaplo.A + count(momGenotype,"A")/2;
            motherHaplo.a = motherHaplo.a + count(momGenotype,"a")/2;

            % fill in each column 
            % AA
            genoArray(i,1) = (fatherHaplo.A)*(motherHaplo.A);
            % Aa
            genoArray(i,2) = (fatherHaplo.A)*(motherHaplo.a) + (fatherHaplo.a)*(motherHaplo.A);
            % aa
            genoArray(i,3) = (fatherHaplo.a)*(motherHaplo.a);
        else
            % parental sterility (no progeny are produced)
            genoArray(i,:) = 0;
        end        
        
    end % end of for 1:numPairings loop

%% recessive sterile--hemizygous females + males fertile, homozygotes sterile
elseif driveType == 2
    % recessive female sterile
    
    for i = 1:numPairings
        fatherHaplo = struct();    
        fatherHaplo.A = 0;
        fatherHaplo.a = 0;

        motherHaplo = fatherHaplo;    

        dadGenotype = string(tmpArray(1+i,1));
        momGenotype = string(tmpArray(1+i,2));

        % first make sure that progeny are actually produced
        nonviablePairing = (dadGenotype == "AA" || momGenotype == "AA");
        if (~nonviablePairing)
            % progeny are produced, check for homing in either parent
            homingDadBool = (dadGenotype == "Aa");
            homingMomBool = (momGenotype == "Aa");

            if (homingDadBool)
                % with MALE_CONVERSION_PROB, then father's genotype is AA
                fatherHaplo.A = (1-MALE_CONVERSION_PROB)*0.5 + MALE_CONVERSION_PROB;
                fatherHaplo.a = (1-MALE_CONVERSION_PROB)*0.5;
            else 
                % homing is not possible
                fatherHaplo.A = fatherHaplo.A + count(dadGenotype,"A")/2;
                fatherHaplo.a = fatherHaplo.a + count(dadGenotype,"a")/2;
            end
            
            if (homingMomBool)
                % with FEMALE_CONVERSION_PROB, then mother's genotype is CC
                motherHaplo.A = (1-FEMALE_CONVERSION_PROB)*0.5 + FEMALE_CONVERSION_PROB;
                motherHaplo.a = (1-FEMALE_CONVERSION_PROB)*0.5;
            else 
                % homing is not possible
                motherHaplo.A = motherHaplo.A + count(momGenotype,"A")/2;
                motherHaplo.a = motherHaplo.a + count(momGenotype,"a")/2;
            end             

            % fill in each column 
            % AA
            genoArray(i,1) = (fatherHaplo.A)*(motherHaplo.A);
            % Aa
            genoArray(i,2) = (fatherHaplo.A)*(motherHaplo.a) + (fatherHaplo.a)*(motherHaplo.A);
            % aa
            genoArray(i,3) = (fatherHaplo.a)*(motherHaplo.a);

            % is the mother hemizygous? 
            if (momGenotype == "Aa")
                % reduced fecundity
                genoArray(i,:) = genoArray(i,:) * RELATIVE_FEMALE_FECUND;
                % progeny are effected by maternal deposition (only occurs
                % when mom has Cas9/gRNA complex)
                motherDepoInd(i) = 1; 
            end
            
            
        else
            % parental sterility (no progeny are produced)
            genoArray(i,:) = 0;
        end         
    
    end % end of for 1:numPairings loop
    
else
    error("Incorrect driveType. Must either be 1 (dominant female sterile) or 2 (recessive sterile).");
end

% append a column for death probability
genoArray = [genoArray, 1-sum(genoArray,2)];
% correct numerical error
genoArray(abs(genoArray) < 10^(-8)) = 0; 


%% return
breedingTables = struct();
breedingTables.genoArray = genoArray;
breedingTables.motherDepoInd = motherDepoInd;

end