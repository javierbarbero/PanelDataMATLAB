function [ est ] = defaultVarNames( est )
%DEFAULTVARNAMES Internal Function
%   Internal Function
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%

    % Dependent variable
    if ~iscell(est.ynames)
        est.ynames = cellstr(sprintf('deptvar'));
    end
    
    % Independent variables
    korig = size(est.X,2);
    k = est.k;
    lxnames = length(est.xnames);
    
    if ~iscell(est.xnames) || lxnames < korig              
        
        % Convert to cell if it is not cell
        if ~iscell(est.xnames)            
            est.xnames = cell(k,1);
            begin = 0;
        else
            % Variables that already have names
            begin = lxnames;
        end
        
        % Fill missing varnames        
        for i=begin+1:1:korig
            est.xnames(i) = cellstr(sprintf('var%d',i));
        end
    end
    
    % Spatial    
    if est.isSpatial
        klast = korig;
        if est.slagy            
            klast = klast + 1;
            est.xnames(klast) = cellstr(sprintf('W*%s',char(est.ynames(1))));
        end
        if est.slagX
            lslagX = length(est.slagX);
            for i=1:lslagX
                klast = klast + 1;
                est.xnames(klast) = cellstr(sprintf('W*%s',char(est.xnames(est.slagX(i)))));
            end
        end
    end
    
    % Constant term
    if est.hasConstant
        if lxnames < k
            est.xnames(k) = cellstr('CONST');
        end
    end
    
    
    % ---------------
    %   Instruments
    % ---------------
    if est.isInstrumental
        lneworig = size(est.Z,2);
        lnew = est.lnew;
        lznames = length(est.znames);
        
        % Convert to cell if it is not cell
        if ~iscell(est.znames)            
            est.znames = cell(lnew,1);
            begin = 0;
        else
            % Variables that already have names
            begin = lznames;
        end
        
        % Fill missing varnames        
        for i=begin+1:1:lneworig
            est.znames(i) = cellstr(sprintf('iv%d',i));
        end
    end
    
    % OLD
    % Create cell to store names
    %{
    est.xnames = cell(est.k,est.g);
    pos = 1;
    for i=1:est.g
        if est.hasConstant
            last = est.keq(i)-1;
        else
            last = est.keq(i);
        end

        % Write varX for each variable
        for i=1:last
            est.xnames(pos) = cellstr(sprintf('var%d',i));
            pos = pos+1;
        end

        % Add constant name if not specified
        if est.hasConstant       
            est.xnames(pos) = cellstr('CONSTANT');
        end
        pos = pos+1;
    end
    
    % TODO: Add multi-equation to instruments
    % Names of instruments if IV estimation
    if est.isInstrumental()
        est.znames = cell(est.lnew,1);
        
        for i=1:est.lnew
            est.znames(i) = cellstr(sprintf('iv%d',i));
        end
    end
    %}
    

end