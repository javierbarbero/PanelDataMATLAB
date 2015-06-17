function [ Xrep ] = replicate_elements( X, Tid )
%REPLICATE_ELEMENTS Internal Function
%   Internal Function
%
%   Copyright 2013-2015 Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 17, June, 2015
%
	if ~verLessThan('matlab', '8.5')
        % Use the new built in function 'repelem'
        Xrep = repelem(X,Tid,1);
    else
        % Before R2015a (no 'repelem' function')
        k = size(X,2);
        Xrep = nan(sum(Tid),k);
        pos = 1;
        for i=1:length(Tid)
            Xrep(pos:pos+Tid(i)-1,1:k) = repmat(X(i,1:k),Tid(i),1);
            pos = pos+Tid(i);
        end
    end

end

