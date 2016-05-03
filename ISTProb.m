function [new_PCorrect, old_PCorrect] = ISTProb(nColour1, nColour2, chosenColour, N)
%
% Usage: [new_PCorrect, old_PCorrect] = ISTProb(nColour1, nColour2, chosenColour, N)
%
% where 
%   nColour1 is the number of open boxes of colour 1 at time of choice
%   nColour2 is the number of open boxes of colour 2 at time of choice
%   chosenColour is the colour chosen by the participant (either 1 or 2)
%   N is the total number of boxes in the grid (typically 25)
%
% This function is identical to the helper function contained in
% FormatCANTAB.m, and is included here for convenience.
%
% Contact: Daniel Bennett, danielbrianbennett (at) gmail (dot) com
%
% This code is distributed under a GNU GPL license as documented at
% http://www.gnu.org/licenses/gpl.txt. Software is distributed as is,
% completely without warranty or service support. The authors of the
% software are not liable for the condition or performance of the code, and
% hereby disclaim all implied warranties.

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
        for i = A:Z       
            aCoeff(i) = nchoosek(Z,i);       
        end    
        old_PCorrect = sum(aCoeff) ./ (2 ^ Z);
    end

end
end