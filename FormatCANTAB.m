function FormatCANTAB(inFilename,outFilename)
%
% Usage: FormatCANTAB(input_filename, output_filename)
%
% e.g. FormatCANTAB('raw_test_data.xlsx','processed_test_data.xslx')
%
% This function uses MATLAB's built-in xlsread and xlswrite functions to
% import data from a raw CANTAB output file (the 'input_filename' argument)
% and export it in processed form (the 'output_filename' argument). 
%
% A couple of notes: (1) Both input and output files must be Microsoft
% Excel spreadsheets. (2) Because of limitations in the MATLAB functions 
% that handle reading and writing from excel files, this function will only
% work in Windows, and not on Macs. We wrote and tested it using MATLAB
% version R2015b and Microsoft Excel 2010 for Windows. While we would not
% expect performance to be substantially different in other versions of
% these programs, we have not tested these and make no guarantees regarding
% the software.
%
% Contact: Daniel Bennett, danielbrianbennett (at) gmail (dot) com
%
% This code is distributed under a GNU GPL license as documented at
% http://www.gnu.org/licenses/gpl.txt. Software is distributed as is,
% completely without warranty or service support. The authors of the
% software are not liable for the condition or performance of the code, and
% hereby disclaim all implied warranties.

%% Import the data and reformat it
[~, ~, raw] = xlsread(inFilename);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
cellVectors = raw(:,[2,3,7,10,11,12,14,16]);
raw = raw(:,[1,4,5,6,8,9,13,15,17,18]);

% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

% Create output variable
data = reshape([raw{:}],size(raw));

% Allocate imported array to column variable names
blocknumber = data(:,1);
assessed = cellVectors(:,1);
wincondition = cellVectors(:,2);
trialnumber = data(:,4);
passed = cellVectors(:,3);
colourchosen = cellVectors(:,5);
majoritycolourchosen = cellVectors(:,6);
sequence = cellVectors(:,7);
openedbox = data(:,8);
colourrevealed = cellVectors(:,8);

% Clear temporary variables
clearvars data raw cellVectors R;

%% work out block start locations
blocks = unique(blocknumber(~isnan(blocknumber))); % find different block numbers
for i = 1:numel(blocks)
    blockLocs(i) = find(blocknumber == blocks(i)); % find row indices corresponding to start of each block
end

fixedBlockLoc = blockLocs(strcmp(assessed(blockLocs),'yes') & strcmp(wincondition(blockLocs),'FIXED')); % row index of start of fixed block
decreasingBlockLoc = blockLocs(strcmp(assessed(blockLocs),'yes') & strcmp(wincondition(blockLocs),'DECREASING')); % row index of start of decreasing block

%% work out where trials start in each block
if numel(fixedBlockLoc) > 1 || numel(decreasingBlockLoc) > 1
    error('Too many blocks found.')
elseif fixedBlockLoc < decreasingBlockLoc
    firstCondition = 'fixed';
    fixedTrialLocs = fixedBlockLoc + find(~isnan(trialnumber(fixedBlockLoc:decreasingBlockLoc-1))) - 1;
    decreasingTrialLocs = decreasingBlockLoc + find(~isnan(trialnumber(decreasingBlockLoc:end))) - 1;
elseif decreasingBlockLoc < fixedBlockLoc
    firstCondition = 'decreasing';
    decreasingTrialLocs = decreasingBlockLoc + find(~isnan(trialnumber(decreasingBlockLoc:fixedBlockLoc-1))) - 1;
    fixedTrialLocs = fixedBlockLoc + find(~isnan(trialnumber(fixedBlockLoc:end))) - 1;
end

%% work out how many boxes were opened in each trial
for i = 1:numel(fixedTrialLocs)   
    temp = find(openedbox(fixedTrialLocs(i):numel(openedbox)) == 1,2);
    if numel(temp) > 1
        nBoxesFixed(i) = temp(2)-1; 
    elseif numel(temp) == 1
        nBoxesFixed(i) = numel(openedbox(~isnan(openedbox))) - fixedTrialLocs(i) + 1;
    end    
    fixedBoxIndices(i,:) = [fixedTrialLocs(i), fixedTrialLocs(i) + nBoxesFixed(i) - 1];    
    clear temp
end

for i = 1:numel(decreasingTrialLocs)   
    temp = find(openedbox(decreasingTrialLocs(i):numel(openedbox)) == 1,2);
    if numel(temp) > 1
        nBoxesDecreasing(i) = temp(2)-1; 
    elseif numel(temp) == 1
        nBoxesDecreasing(i) = numel(openedbox(~isnan(openedbox))) - decreasingTrialLocs(i) + 2;
    end   
    decreasingBoxIndices(i,:) = [decreasingTrialLocs(i), decreasingTrialLocs(i) + nBoxesDecreasing(i) - 1];
    clear temp    
end

%% work out number of each colour
for i = 1:size(fixedBoxIndices,1)
    S = regexp(sequence{fixedTrialLocs(i)},';[\s-]','split');
    fixedColours{i} = unique(S);
    nColour1_fixed(i) = sum(strcmp({colourrevealed{fixedBoxIndices(i,1):fixedBoxIndices(i,2)}},fixedColours{i}(1)));
    nColour2_fixed(i) = sum(strcmp({colourrevealed{fixedBoxIndices(i,1):fixedBoxIndices(i,2)}},fixedColours{i}(2)));
    colourChosen_fixed(i) = find(ismember(fixedColours{i},colourchosen{fixedTrialLocs(i)}));
end

for i = 1:size(decreasingBoxIndices,1)
    S = regexp(sequence{decreasingTrialLocs(i)},';[\s-]','split');
    decreasingColours{i} = unique(S);
    nColour1_decreasing(i) = sum(strcmp({colourrevealed{decreasingBoxIndices(i,1):decreasingBoxIndices(i,2)}},decreasingColours{i}(1)));
    nColour2_decreasing(i) = sum(strcmp({colourrevealed{decreasingBoxIndices(i,1):decreasingBoxIndices(i,2)}},decreasingColours{i}(2)));
    colourChosen_decreasing(i) = find(ismember(decreasingColours{i},colourchosen{decreasingTrialLocs(i)}));
end

%% write data to output spreadsheet
columnNames = {'trial','condition','block number', 'number of boxes opened','n colour 1','n colour 2','colour 1 name', 'colour 2 name', 'chosen colour', 'correct?','majority colour chosen?','old P(correct)','new P(correct)'};
xlswrite(outFilename, columnNames);
counter = 2;
% write for fixed condition first
for i = 1:numel(fixedTrialLocs)
   
    if strcmp(firstCondition,'fixed')
        blockNumber = 1;
    elseif strcmp(firstCondition,'decreasing')
        blockNumber = 2;
    else
        error()
    end
    
    [newProb_fixed(i), oldProb_fixed(i)] = ISTProb(nColour1_fixed(i), nColour2_fixed(i), colourChosen_fixed(i) , 25);
    
    writeData = {...
        trialnumber(fixedTrialLocs(i))              % trial number
        'fixed'                                     % win condition
        blockNumber                                 % whether first or second block
        nColour1_fixed(i) + nColour2_fixed(i)       % n boxes opened in total
        nColour1_fixed(i)                           % n opened of colour 1
        nColour2_fixed(i)                           % n opened of colour 2
        fixedColours{i}{1}                          % name of colour 1
        fixedColours{i}{2}                          % name of colour 2
        colourChosen_fixed(i)                       % which colour was chosen?
        passed{fixedTrialLocs(i)}                   % correct? (yes = correct, no = incorrect)
        majoritycolourchosen{fixedTrialLocs(i)}     % did the participant choose the colour in the majority of open boxes?
        oldProb_fixed(i)                            % P(correct) calculated as per Clark et al., 2006
        newProb_fixed(i)                            % P(correct) calculated as per modified formula
        }';
    
    xlswrite(outFilename,writeData,1,sprintf('A%.0f',counter))
    counter = counter + 1;
end

% write for decreasing condition second
for i = 1:numel(decreasingTrialLocs)
   
    if strcmp(firstCondition,'decreasing')
        blockNumber = 1;
    elseif strcmp(firstCondition,'fixed')
        blockNumber = 2;
    else
        error()
    end
    
    % Calculate different P(correct) measures using the ISTProb helper function (included below)
    [newProb_decreasing(i), oldProb_decreasing(i)] = ISTProb(nColour1_decreasing(i), nColour2_decreasing(i), colourChosen_decreasing(i) , 25);
    
    writeData = {...
        trialnumber(decreasingTrialLocs(i))              % trial number
        'decreasing'                                     % win condition
        blockNumber                                      % whether first or second block
        nColour1_decreasing(i) + nColour2_decreasing(i)  % n boxes opened in total
        nColour1_decreasing(i)                           % n opened of colour 1
        nColour2_decreasing(i)                           % n opened of colour 2
        decreasingColours{i}{1}                          % name of colour 1
        decreasingColours{i}{2}                          % name of colour 2
        colourChosen_decreasing(i)                       % which colour was chosen?
        passed{decreasingTrialLocs(i)}                   % correct? (yes = correct, no = incorrect)
        majoritycolourchosen{decreasingTrialLocs(i)}     % did the participant choose the colour in the majority of open boxes?
        oldProb_decreasing(i)                            % P(correct) calculated as per Clark et al., 2006
        newProb_decreasing(i)                            % P(correct) calculated as per modified formula
        }';
    
    xlswrite(outFilename,writeData,1,sprintf('A%.0f',counter))
    counter = counter + 1;
end

xlswrite(outFilename,{'Old Mean of P(correct) - fixed condition:'},1,'O2');
xlswrite(outFilename,{'New Mean of P(correct) - fixed condition:'},1,'O3');
xlswrite(outFilename,{'Old Mean of P(correct) - decreasing condition:'},1,'O5');
xlswrite(outFilename,{'New Mean of P(correct) - decreasing condition:'},1,'O6');
xlswrite(outFilename,mean(oldProb_fixed),1,'P2');
xlswrite(outFilename,mean(newProb_fixed),1,'P3');
xlswrite(outFilename,mean(oldProb_decreasing),1,'P5');
xlswrite(outFilename,mean(newProb_decreasing),1,'P6');

%% helper function
function [new_PCorrect, old_PCorrect] = ISTProb(nColour1, nColour2, chosenColour, N)
%
% Usage: [new_PCorrect, old_PCorrect] = ISTProb(nColour1, nColour2, chosenColour, N)
% Check inputs

    if isnan(nColour1) || isnan(nColour2) || isnan(chosenColour)
        new_PCorrect = NaN;
        old_PCorrect = NaN;
    else

    % Derive some useful variables
    nRequiredForMajority = ceil(N/2);
    probIndices = 0:N;

    if chosenColour == 1
        nChosenColour = nColour1;
        nUnchosenColour = nColour2;
    elseif chosenColour == 2
        nChosenColour = nColour2;
        nUnchosenColour = nColour1;
    else
        error('Something went wrong.');
    end

    %% Calculate New P(correct) Measure

    % Calculate the new P(correct) measure using a ratio of hypergeometric
    % distributions. This is equivalent to the equation presented in the
    % commentary, but substantially faster. To verify the equivalence of
    % these measures, you can uncomment the block of code below titled
    % 'Alternative Calculation of new P(Correct) measure'.
    outputPMF = hygepdf(nChosenColour,N,probIndices,nColour1+nColour2) ./ sum(hygepdf(nChosenColour,N,probIndices,nColour1+nColour2)  );
    new_PCorrect = sum(outputPMF(probIndices >= nRequiredForMajority));
    
    % %% Alternative calculation of new P(Correct) measure
    % % This code explicitly implements the formula in the Bennett et al. commentary.   
    %   
    %     M = 13:25;
    %     denominatorIndices = nChosenColour:(25-nUnchosenColour);
    %     Pr = nan(size(M));     
    %     for i = 1:numel(M)        
    %         denominator = nan(size(denominatorIndices));       
    %         for j = 1:numel(denominatorIndices)         
    %             denominator(j) = nchoosek(denominatorIndices(j), nChosenColour) * nchoosek(25-denominatorIndices(j),nUnchosenColour);          
    %         end        
    %         try
    %             Pr(i) = ( nchoosek(M(i),nChosenColour) * nchoosek(25-M(i),nUnchosenColour) ) / sum(denominator);
    %         catch
    %             Pr(i) = 0;
    %         end
    %     end   
    %     alternative_PCorrect = sum(Pr);

    %% Calculate P(correct) as per Clark et al. (2006)

    Z = N - (nColour1 + nColour2);
    A = nRequiredForMajority - nChosenColour;

    if nChosenColour >= nRequiredForMajority
        old_PCorrect = 1;
    elseif A > Z
        old_PCorrect = 0;
    else
        for ii = A:Z       
            aCoeff(ii) = nchoosek(Z,ii);       
        end    
        old_PCorrect = sum(aCoeff) ./ (2 ^ Z);
    end

end
end
end