function [ p ] = normalden( x, mu, sigma )
%NORMALDEN Return the probability of a normal function
%   Return the probability of a normal function
%
%   AUTHORS: Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 9, June, 2015
%
    
    % Normal(0,1)
    if nargin < 2
        mu = 0;
        sigma = 1;
    end
    
    % Normal(mu,1)
    if nargin < 3
        sigma = 1;
    end

    p = 1/(sqrt(2*pi)*sigma) * exp(-1/2 * ((x-mu)./sigma).^2);

end

