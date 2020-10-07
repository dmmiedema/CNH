%% MATLAB function to calculate copy number intra-tumor heterogeneity (CNH) from a single copy number measurement.
% In this method, relative segmented copy numbers are transformed to absolute copy numbers and the distance
% to integer values is measured for each segment. This procedure is
% performed for either a range of ploidies/purities or for fixed
% ploidy/purity, in case either or both quantities are known from previous
% measurements.
% Copy number heterogeneity (CNH) is defined as the minimum average distance of
% segments to the closest integer value, weighted by segment length. 
% In addition to CNH, this function returns the ploidy and purity
% corresponding to the inferred CNH.

%% Input
% 1st input argument:   (seg_val)       values of relative copy numbers per segment, vector of size Nx1
% 2nd input argument:   (seg_len)       segment lengths, vector of size Nx1
% 3th input argument:   (ploidy)        tumour ploidy, use empty vector ( = [] ) for grid search 
% 4th input argument:   (purity)        sample purity, use empty vector ( = [] ) for grid search     

%% Output
% 1st output argument:  (CNH_out)       inferred CNH 
% 2nd output argument:  (ploidy_out)    inferred ploidy for empty input ploidy, otherwise same as input ploidy. 
% 3th output argument:  (purity_out)    inferred purity for empty input purity, otherwise same as input purity.

%% Main function
function [ CNH_out, ploidy_out, purity_out ]= CNH( seg_val ,seg_len, ploidy, purity ) 
    % check if input seg_val and seg_len are vectors of equal size
    if ~( isvector(seg_val) && isvector(seg_len) && numel(seg_val) == numel(seg_len) && size(seg_val,2) == 1)
        error('Segment values (1st input argument) and segment lengths (second input argument) appear not to be column vectors of equal length');
    end    
    % check if ploidy is scalar argument or empty
    if ~( (isscalar(ploidy) && ploidy > 0) || isempty(ploidy) )
        error('Ploidy is not a positive scalar or empty');
    end
    % check if purity is a scalar between 0 and 1 argument or empty
    if ~( (isscalar(purity) && purity > 0 && purity <=1 ) || isempty(purity) )
        error('Purity is not a positive scalar or empty');
    end
    
    % specify default range of ploidy for grid search, if input ploidy is empty
    if isempty(ploidy)
        ploidy = 1.5:0.01:5;   % tumor ploidy
    end
    % specify default range of purity purity for grid search, if input purity is empty
    if isempty(purity)
        purity = 0.2:0.01:1;  % tumor purity: standard range 
    end
    
    % get number of ploidies and purities for grid search
    Nploidy = numel(ploidy);
    Npurity = numel(purity);
        
    % initialize vectors a1 and a2 from all combinations of ploidy and purity, for the transformation of 
    % measured relative copy number profile (seg_val) to absolute values (q) using   
    %  q = seg_val*a1+a2. 
    a1 = zeros(Nploidy*Npurity,1);
    a2 = zeros(Nploidy*Npurity,1);
    purity_all = zeros(Nploidy*Npurity,1); % vector that contains all purities used in 2D grid search
    ploidy_all = zeros(Nploidy*Npurity,1); % vector that contains all ploidies used in 2D grid search
    for i = 1:Nploidy
        a1((i-1)*Npurity+1:i*Npurity) = (purity*ploidy(i)+2*(1-purity))./purity;
        a2((i-1)*Npurity+1:i*Npurity) = -2*(1-purity)./purity;
        
        purity_all( (i-1)*Npurity+1:i*Npurity) = purity;
        ploidy_all( (i-1)*Npurity+1:i*Npurity) = ploidy(i)*ones(Npurity,1);
    end    
        
    % iniatilize output: CNH_out, ploidy_out and purity_out
    CNH_out = 1;
    purity_out = 0;
    ploidy_out = 0;
    
    % grid search over all ploidies and purities to infer CNH
    for i = 1:Nploidy*Npurity
        % transform relative copy numbers to absolute numbers
        q = a1(i)*seg_val+a2(i);
        
        % measure distance to closest integer value of each segment
        q_dist_down = mod(q,1); 
        q_dist_up = 1-mod(q,1);
        q_dist_min = min(q_dist_up,q_dist_down);
        
        % calculate the mean distance of segments to integer values,
        % weighted by segment length
        CNHnew = sum(q_dist_min.*seg_len) / (sum(seg_len));
        
        % if the new weighted distance to integers is smaller than any
        % previously calculated, replace CNH_out, ploidy_out and purity_out with the new values.
        if CNHnew < CNH_out
            CNH_out = CNHnew;
            purity_out = purity_all(i);
            ploidy_out = ploidy_all(i);
        end
    end
end