% Function algorithmically producing breeding tables for a simple
% tethered homing drive with dominant and recessive sterility options.
% Author: Cole Butler 
%
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% The following function produces a breeding table given conversion
% efficiencies for either sex and relative female fecundity for hemizygous
% females. Organisms homozygous for the construct, denoted "A", are
% sterile. The Cas9 is X-linked, in which case it is denoted "Z". 
%
% INPUTS:
% -- MALE_CONVERSION_PROB: Conversion efficiency of construct in males.
% -- FEMALE_CONVERSION_PROB: Conversion efficiency of construct in females.
% -- RELATIVE_FEMALE_FECUND: Relative fecundity of hemizygous females
%   (hemizygous males are assumed comparable to wild-type).
% -- driveType: The type of drive being considered, 1 for dominant sterile 
%   and 2 for recessive sterile
%
% OUTPUTS:
% -- breedingTable: A structure array containing the male and female
%   breeding tables
% -- motherDepoInd: Vector w/indices of pairings affected by maternal
%   deposition of Cas9
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [breedingTables] = genBT_simple_tethered_homing(MALE_CONVERSION_PROB, ...
    FEMALE_CONVERSION_PROB, RELATIVE_FEMALE_FECUND, driveType)

%% produce breeding table for male zygotes
tmpArray = readcell('test_male_breeding_table.xlsx');

% pull out the part of the breeding table we will update
male_genoArray = tmpArray(2:end,3:end);
arrayDim = size(male_genoArray);
numPairings = arrayDim(1);
numGenotypes = arrayDim(2);
% include a column for nonviable progeny
male_genoArray = zeros(numPairings, numGenotypes); 

%% calculate zygote frequencies producing males
if (driveType == 1)
    % dominant sterile--hemizygous males fertile, homozygous and 
    % hemizygous females sterile
    
    % for each possible pairing, determine zygote frequencies...
    for i = 1:numPairings
        % initialize haplotypes @ 0 (for male progeny, the father never
        % contributes a X chromosome)
        fatherHaplo = struct();    
        fatherHaplo.YA = 0;
        fatherHaplo.Ya = 0;

        motherHaplo = struct();
        motherHaplo.XA = 0;
        motherHaplo.Xa = 0;
        motherHaplo.ZA = 0;
        motherHaplo.Za = 0;

        dadGenotype = string(tmpArray(1+i,1));
        momGenotype = string(tmpArray(1+i,2));

        % extract locus genotypes for each parent
        dadChrGeno = extractBetween(dadGenotype,1,2);
        dadAllGeno = extractBetween(dadGenotype,3,4);
        momChrGeno = extractBetween(momGenotype,1,2);
        momAllGeno = extractBetween(momGenotype,3,4);

        % first make sure that progeny are actually produced
        nonviablePairing = (dadAllGeno == "AA" || momAllGeno == "AA" || ...
            momAllGeno == "Aa");
        if (~nonviablePairing)
            % progeny are produced, check for homing in either parent
            homingDadBool = (dadAllGeno == "Aa") && (contains(dadChrGeno,"Z"));

            if (homingDadBool)
                % with MALE_CONVERSION_PROB, then father's genotype is ZYAa;
                % note that male progeny are guaranteeed to receive a Y...
                fatherHaplo.YA = MALE_CONVERSION_PROB + fatherHaplo.YA;
                % homing fails...
                fatherHaplo.YA = (1-MALE_CONVERSION_PROB)*0.5 + fatherHaplo.YA;
                fatherHaplo.Ya = (1-MALE_CONVERSION_PROB)*0.5 + fatherHaplo.Ya;
            else 
                % homing is not possible, genotype is XYAA, XYAa, XYaa,
                % ZYaa, ZYAA
                fatherHaplo.YA = fatherHaplo.YA + count(dadAllGeno,"A")/2;
                fatherHaplo.Ya = fatherHaplo.Ya + count(dadAllGeno,"a")/2;
            end
            
            % homing does not occur in females, genotype is XXaa, ZZaa,
            % etc.
            motherHaplo.XA = motherHaplo.XA + ((count(momChrGeno,"X"))/2)*((count(momAllGeno,"A"))/2);
            motherHaplo.Xa = motherHaplo.Xa + ((count(momChrGeno,"X"))/2)*((count(momAllGeno,"a"))/2);
            motherHaplo.ZA = motherHaplo.ZA + ((count(momChrGeno,"Z"))/2)*((count(momAllGeno,"A"))/2);
            motherHaplo.Za = motherHaplo.Za + ((count(momChrGeno,"Z"))/2)*((count(momAllGeno,"a"))/2);

            % fill in each column (for males)
            % XYAA = YA + XA
            male_genoArray(i,1) =  (fatherHaplo.YA)*(motherHaplo.XA);
            % XYAa = YA + Xa, Ya + XA
            male_genoArray(i,2) =  (fatherHaplo.YA)*(motherHaplo.Xa) + (fatherHaplo.Ya)*(motherHaplo.XA);
            % XYaa = Ya + Xa
            male_genoArray(i,3) =  (fatherHaplo.Ya)*(motherHaplo.Xa);
            % ZYAA = YA + ZA
            male_genoArray(i,4) =  (fatherHaplo.YA)*(motherHaplo.ZA);
            % ZYAa = YA + Za, Ya + ZA
            male_genoArray(i,5) =  (fatherHaplo.YA)*(motherHaplo.Za) + (fatherHaplo.Ya)*(motherHaplo.ZA);
            % ZYaa = Ya + Za
            male_genoArray(i,6) =  (fatherHaplo.Ya)*(motherHaplo.Za);
        else
            % parental sterility (no progeny are produced)
            male_genoArray(i,:) = 0;
        end    


    end
elseif (driveType == 2)
    % recessive sterile--hemizygotes fertile, homozygotes sterile
    
    % for each possible pairing, determine zygote frequencies...
    for i = 1:numPairings
        % initialize haplotypes @ 0 (for male progeny, the father never
        % contributes a X chromosome)
        fatherHaplo = struct();    
        fatherHaplo.YA = 0;
        fatherHaplo.Ya = 0;

        motherHaplo = struct();
        motherHaplo.XA = 0;
        motherHaplo.Xa = 0;
        motherHaplo.ZA = 0;
        motherHaplo.Za = 0;

        dadGenotype = string(tmpArray(1+i,1));
        momGenotype = string(tmpArray(1+i,2));

        % extract locus genotypes for each parent
        dadChrGeno = extractBetween(dadGenotype,1,2);
        dadAllGeno = extractBetween(dadGenotype,3,4);
        momChrGeno = extractBetween(momGenotype,1,2);
        momAllGeno = extractBetween(momGenotype,3,4);

        % first make sure that offspring are actually produced
        nonviablePairing = (dadAllGeno == "AA" || momAllGeno == "AA");
        if (~nonviablePairing)
            % progeny are produced, check for homing in either parent
            homingDadBool = (dadAllGeno == "Aa") && (contains(dadChrGeno,"Z"));
            homingMomBool = (momAllGeno == "Aa") && (contains(momChrGeno,"Z"));

            if (homingDadBool)
                % with MALE_CONVERSION_PROB, then father's genotype is ZYAa;
                % note that male progeny are guaranteeed to receive a Y...
                fatherHaplo.YA = MALE_CONVERSION_PROB + fatherHaplo.YA;
                % homing fails...
                fatherHaplo.YA = (1-MALE_CONVERSION_PROB)*0.5 + fatherHaplo.YA;
                fatherHaplo.Ya = (1-MALE_CONVERSION_PROB)*0.5 + fatherHaplo.Ya;
            else 
                % homing is not possible, genotype is XYAA, XYAa, XYaa,
                % ZYaa, ZYAA
                fatherHaplo.YA = fatherHaplo.YA + count(dadAllGeno,"A")/2;
                fatherHaplo.Ya = fatherHaplo.Ya + count(dadAllGeno,"a")/2;
            end

            if (homingMomBool)
                % with FEMALE_CONVERSION_PROB, then mother's genotype is
                % ZZAa or ZXAa
                motherHaplo.XA = FEMALE_CONVERSION_PROB*(count(momChrGeno,"X")/2) + motherHaplo.XA;
                motherHaplo.ZA = FEMALE_CONVERSION_PROB*(count(momChrGeno,"Z")/2) + motherHaplo.ZA;
                % homing fails...
                motherHaplo.XA = (1-FEMALE_CONVERSION_PROB)*(count(momChrGeno,"X")/2)*0.5 + motherHaplo.XA;
                motherHaplo.Xa = (1-FEMALE_CONVERSION_PROB)*(count(momChrGeno,"X")/2)*0.5 + motherHaplo.Xa;
                motherHaplo.ZA = (1-FEMALE_CONVERSION_PROB)*(count(momChrGeno,"Z")/2)*0.5 + motherHaplo.ZA;
                motherHaplo.Za = (1-FEMALE_CONVERSION_PROB)*(count(momChrGeno,"Z")/2)*0.5 + motherHaplo.Za;                
            else 
                % homing is not possible, genotype is XXAA, XXAa, XXaa,
                % ZZaa, ZZAA, ZXAA, ZXaa, etc.
                motherHaplo.XA = motherHaplo.XA + ((count(momChrGeno,"X"))/2)*((count(momAllGeno,"A"))/2);
                motherHaplo.Xa = motherHaplo.Xa + ((count(momChrGeno,"X"))/2)*((count(momAllGeno,"a"))/2);
                motherHaplo.ZA = motherHaplo.ZA + ((count(momChrGeno,"Z"))/2)*((count(momAllGeno,"A"))/2);
                motherHaplo.Za = motherHaplo.Za + ((count(momChrGeno,"Z"))/2)*((count(momAllGeno,"a"))/2);
            end            



            % fill in each column (for males)
            % XYAA = YA + XA
            male_genoArray(i,1) =  (fatherHaplo.YA)*(motherHaplo.XA);
            % XYAa = YA + Xa, Ya + XA
            male_genoArray(i,2) =  (fatherHaplo.YA)*(motherHaplo.Xa) + (fatherHaplo.Ya)*(motherHaplo.XA);
            % XYaa = Ya + Xa
            male_genoArray(i,3) =  (fatherHaplo.Ya)*(motherHaplo.Xa);
            % ZYAA = YA + ZA
            male_genoArray(i,4) =  (fatherHaplo.YA)*(motherHaplo.ZA);
            % ZYAa = YA + Za, Ya + ZA
            male_genoArray(i,5) =  (fatherHaplo.YA)*(motherHaplo.Za) + (fatherHaplo.Ya)*(motherHaplo.ZA);
            % ZYaa = Ya + Za
            male_genoArray(i,6) =  (fatherHaplo.Ya)*(motherHaplo.Za);


            % is the mother hemizygous? 
            if (momAllGeno == "Aa")
                % reduced fecundity
                male_genoArray(i,:) = male_genoArray(i,:) * RELATIVE_FEMALE_FECUND;
            end

        else
            % parental sterility (no progeny are produced)
            male_genoArray(i,:) = 0;
        end    

    end % end of for numPairings loop
else
    error("Incorrect driveType. Must either be 1 (dominant female sterile) or 2 (recessive sterile).");
end

% append a column for death probability
male_genoArray = [male_genoArray, 1-sum(male_genoArray,2)];


%% now do the same for the female genotypes...
tmpArray = readcell('test_female_breeding_table.xlsx');

% pull out the part of the breeding table we will update
female_genoArray = tmpArray(2:end,3:end);
arrayDim = size(female_genoArray);
numPairings = arrayDim(1);
numGenotypes = arrayDim(2);
% include a column for nonviable progeny
female_genoArray = zeros(numPairings, numGenotypes);
% vector storing indices of pairings producing female offspring that are
% affected by maternal deposition
motherDepoInd = zeros(1,numPairings); 

%% calculate zygote frequencies producing females
if (driveType == 1)
    % dominant sterile--hemizygous males fertile, homozygotes and 
    % hemizygous females sterile

    % for each possible pairing, determine zygote frequencies...
    for i = 1:numPairings
        % initialize haplotypes @ 0 (for male progeny, the father never
        % contributes a Y chromosome)
        fatherHaplo = struct();    
        fatherHaplo.XA = 0;
        fatherHaplo.Xa = 0;
        fatherHaplo.ZA = 0;
        fatherHaplo.Za = 0;    

        motherHaplo = struct();
        motherHaplo.XA = 0;
        motherHaplo.Xa = 0;
        motherHaplo.ZA = 0;
        motherHaplo.Za = 0;

        dadGenotype = string(tmpArray(1+i,1));
        momGenotype = string(tmpArray(1+i,2));

        % extract locus genotypes for each parent
        dadChrGeno = extractBetween(dadGenotype,1,2);
        dadAllGeno = extractBetween(dadGenotype,3,4);
        momChrGeno = extractBetween(momGenotype,1,2);
        momAllGeno = extractBetween(momGenotype,3,4);

        % first make sure that offspring are actually produced
        nonviablePairing = (dadAllGeno == "AA" || momAllGeno == "AA" || ...
            momAllGeno == "Aa");
        if (~nonviablePairing)
            % progeny are produced, check for homing in either parent
            homingDadBool = (dadAllGeno == "Aa") && (contains(dadChrGeno,"Z"));

            if (homingDadBool)
                % with MALE_CONVERSION_PROB, then father's genotype is ZYAa;
                % note that male progeny are guaranteeed to receive a Y...
                fatherHaplo.ZA = MALE_CONVERSION_PROB + fatherHaplo.ZA;
                % homing fails...
                fatherHaplo.ZA = (1-MALE_CONVERSION_PROB)*0.5 + fatherHaplo.ZA;
                fatherHaplo.Za = (1-MALE_CONVERSION_PROB)*0.5 + fatherHaplo.Za;
            else 
                % homing is not possible, genotype is XYAA, XYAa, XYaa,
                % ZYaa, ZYAA
                fatherHaplo.XA = fatherHaplo.XA + count(dadChrGeno,"X")*count(dadAllGeno,"A")/2;
                fatherHaplo.Xa = fatherHaplo.Xa + count(dadChrGeno,"X")*count(dadAllGeno,"a")/2;
                fatherHaplo.ZA = fatherHaplo.ZA + count(dadChrGeno,"Z")*count(dadAllGeno,"A")/2;
                fatherHaplo.Za = fatherHaplo.Za + count(dadChrGeno,"Z")*count(dadAllGeno,"a")/2;
            end

            
            % no homing in females
            motherHaplo.XA = motherHaplo.XA + ((count(momChrGeno,"X"))/2)*((count(momAllGeno,"A"))/2);
            motherHaplo.Xa = motherHaplo.Xa + ((count(momChrGeno,"X"))/2)*((count(momAllGeno,"a"))/2);
            motherHaplo.ZA = motherHaplo.ZA + ((count(momChrGeno,"Z"))/2)*((count(momAllGeno,"A"))/2);
            motherHaplo.Za = motherHaplo.Za + ((count(momChrGeno,"Z"))/2)*((count(momAllGeno,"a"))/2);

            % fill in each column (for females)
            % XXAA = XA + XA
            female_genoArray(i,1) =  (fatherHaplo.XA)*(motherHaplo.XA);
            % XXAa = XA + Xa, Xa + XA
            female_genoArray(i,2) =  (fatherHaplo.XA)*(motherHaplo.Xa) + (fatherHaplo.Xa)*(motherHaplo.XA);
            % XXaa = Xa + Xa
            female_genoArray(i,3) =  (fatherHaplo.Xa)*(motherHaplo.Xa);
            % ZXAA = ZA + XA, XA + ZA
            female_genoArray(i,4) =  (fatherHaplo.ZA)*(motherHaplo.XA) + (fatherHaplo.XA)*(motherHaplo.ZA);
            % ZXAa = ZA + Xa, Za + XA, XA + Za, Xa + ZA
            female_genoArray(i,5) =  (fatherHaplo.ZA)*(motherHaplo.Xa) + (fatherHaplo.Za)*(motherHaplo.XA) + ... 
                (fatherHaplo.XA)*(motherHaplo.Za) + (fatherHaplo.Xa)*(motherHaplo.ZA);
            % ZXaa = Za + Xa, Xa + Za
            female_genoArray(i,6) =  (fatherHaplo.Za)*(motherHaplo.Xa) + (fatherHaplo.Xa)*(motherHaplo.Za);
            % ZZAA = ZA + ZA
            female_genoArray(i,7) =  (fatherHaplo.ZA)*(motherHaplo.ZA);
            % ZZAa = ZA + Za, Za + ZA
            female_genoArray(i,8) =  (fatherHaplo.ZA)*(motherHaplo.Za) + (fatherHaplo.Za)*(motherHaplo.ZA);
            % ZZaa = Za + Za
            female_genoArray(i,9) =  (fatherHaplo.Za)*(motherHaplo.Za);
            
        else
            % parental sterility (no progeny are produced)
            female_genoArray(i,:) = 0;
        end    

    end % end of for numPairings loop

elseif (driveType == 2)
    % recessive sterile--hemizygotes fertile, homozygotes sterile
    
    % for each possible pairing, determine zygote frequencies...
    for i = 1:numPairings
        % initialize haplotypes @ 0 (for male progeny, the father never
        % contributes a Y chromosome)
        fatherHaplo = struct();    
        fatherHaplo.XA = 0;
        fatherHaplo.Xa = 0;
        fatherHaplo.ZA = 0;
        fatherHaplo.Za = 0;    

        motherHaplo = struct();
        motherHaplo.XA = 0;
        motherHaplo.Xa = 0;
        motherHaplo.ZA = 0;
        motherHaplo.Za = 0;

        dadGenotype = string(tmpArray(1+i,1));
        momGenotype = string(tmpArray(1+i,2));

        % extract locus genotypes for each parent
        dadChrGeno = extractBetween(dadGenotype,1,2);
        dadAllGeno = extractBetween(dadGenotype,3,4);
        momChrGeno = extractBetween(momGenotype,1,2);
        momAllGeno = extractBetween(momGenotype,3,4);

        % first make sure that offspring are actually produced
        nonviablePairing = (dadAllGeno == "AA" || momAllGeno == "AA");
        if (~nonviablePairing)
            % progeny are produced, check for homing in either parent
            homingDadBool = (dadAllGeno == "Aa") && (contains(dadChrGeno,"Z"));
            homingMomBool = (momAllGeno == "Aa") && (contains(momChrGeno,"Z"));

            if (homingDadBool)
                % with MALE_CONVERSION_PROB, then father's genotype is ZYAa;
                % note that male progeny are guaranteeed to receive a Y...
                fatherHaplo.ZA = MALE_CONVERSION_PROB + fatherHaplo.ZA;
                % homing fails...
                fatherHaplo.ZA = (1-MALE_CONVERSION_PROB)*0.5 + fatherHaplo.ZA;
                fatherHaplo.Za = (1-MALE_CONVERSION_PROB)*0.5 + fatherHaplo.Za;
            else 
                % homing is not possible, genotype is XYAA, XYAa, XYaa,
                % ZYaa, ZYAA
                fatherHaplo.XA = fatherHaplo.XA + count(dadChrGeno,"X")*count(dadAllGeno,"A")/2;
                fatherHaplo.Xa = fatherHaplo.Xa + count(dadChrGeno,"X")*count(dadAllGeno,"a")/2;
                fatherHaplo.ZA = fatherHaplo.ZA + count(dadChrGeno,"Z")*count(dadAllGeno,"A")/2;
                fatherHaplo.Za = fatherHaplo.Za + count(dadChrGeno,"Z")*count(dadAllGeno,"a")/2;
            end

            if (homingMomBool)
                % with FEMALE_CONVERSION_PROB, then mother's genotype is
                % ZZAa or ZXAa
                motherHaplo.XA = FEMALE_CONVERSION_PROB*(count(momChrGeno,"X")/2) + motherHaplo.XA;
                motherHaplo.ZA = FEMALE_CONVERSION_PROB*(count(momChrGeno,"Z")/2) + motherHaplo.ZA;
                % homing fails...
                motherHaplo.XA = (1-FEMALE_CONVERSION_PROB)*(count(momChrGeno,"X")/2)*0.5 + motherHaplo.XA;
                motherHaplo.Xa = (1-FEMALE_CONVERSION_PROB)*(count(momChrGeno,"X")/2)*0.5 + motherHaplo.Xa;
                motherHaplo.ZA = (1-FEMALE_CONVERSION_PROB)*(count(momChrGeno,"Z")/2)*0.5 + motherHaplo.ZA;
                motherHaplo.Za = (1-FEMALE_CONVERSION_PROB)*(count(momChrGeno,"Z")/2)*0.5 + motherHaplo.Za;                
            else 
                % homing is not possible, genotype is XXAA, XXAa, XXaa,
                % ZZaa, ZZAA, ZXAA, ZXaa, etc.
                motherHaplo.XA = motherHaplo.XA + ((count(momChrGeno,"X"))/2)*((count(momAllGeno,"A"))/2);
                motherHaplo.Xa = motherHaplo.Xa + ((count(momChrGeno,"X"))/2)*((count(momAllGeno,"a"))/2);
                motherHaplo.ZA = motherHaplo.ZA + ((count(momChrGeno,"Z"))/2)*((count(momAllGeno,"A"))/2);
                motherHaplo.Za = motherHaplo.Za + ((count(momChrGeno,"Z"))/2)*((count(momAllGeno,"a"))/2);
            end            

            % fill in each column (for females)
            % XXAA = XA + XA
            female_genoArray(i,1) =  (fatherHaplo.XA)*(motherHaplo.XA);
            % XXAa = XA + Xa, Xa + XA
            female_genoArray(i,2) =  (fatherHaplo.XA)*(motherHaplo.Xa) + (fatherHaplo.Xa)*(motherHaplo.XA);
            % XXaa = Xa + Xa
            female_genoArray(i,3) =  (fatherHaplo.Xa)*(motherHaplo.Xa);
            % ZXAA = ZA + XA, XA + ZA
            female_genoArray(i,4) =  (fatherHaplo.ZA)*(motherHaplo.XA) + (fatherHaplo.XA)*(motherHaplo.ZA);
            % ZXAa = ZA + Xa, Za + XA, XA + Za, Xa + ZA
            female_genoArray(i,5) =  (fatherHaplo.ZA)*(motherHaplo.Xa) + (fatherHaplo.Za)*(motherHaplo.XA) + ... 
                (fatherHaplo.XA)*(motherHaplo.Za) + (fatherHaplo.Xa)*(motherHaplo.ZA);
            % ZXaa = Za + Xa, Xa + Za
            female_genoArray(i,6) =  (fatherHaplo.Za)*(motherHaplo.Xa) + (fatherHaplo.Xa)*(motherHaplo.Za);
            % ZZAA = ZA + ZA
            female_genoArray(i,7) =  (fatherHaplo.ZA)*(motherHaplo.ZA);
            % ZZAa = ZA + Za, Za + ZA
            female_genoArray(i,8) =  (fatherHaplo.ZA)*(motherHaplo.Za) + (fatherHaplo.Za)*(motherHaplo.ZA);
            % ZZaa = Za + Za
            female_genoArray(i,9) =  (fatherHaplo.Za)*(motherHaplo.Za);        

            % is the mother hemizygous? 
            if (momAllGeno == "Aa")
                % reduced fecundity
                female_genoArray(i,:) = female_genoArray(i,:) * RELATIVE_FEMALE_FECUND;
            end
            
            % check for maternal deposition of Cas9, occurring when mother
            % has Cas9/gRNA complex (both Cas9 and gRNA present)
            if (count(momChrGeno,"Z") > 0) && (count(momAllGeno,"A") > 0)
                motherDepoInd(i) = 1; 
            end

        else
            % parental sterility (no progeny are produced)
            female_genoArray(i,:) = 0;
        end    


    end % end of for numPairings loop    
    
else
    error("Incorrect driveType. Must either be 1 (dominant female sterile) or 2 (recessive sterile).");    
end

% append a column for death probability
female_genoArray = [female_genoArray, 1-sum(female_genoArray,2)];


%% return
breedingTables = struct();
breedingTables.male_genoArray = male_genoArray;
breedingTables.female_genoArray = female_genoArray;
breedingTables.motherDepoInd = motherDepoInd;

end








