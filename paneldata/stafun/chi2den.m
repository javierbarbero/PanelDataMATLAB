function [ p ] = chi2den( x, v )
%CHI2DEN Return the probability of a Chi-Square
%   Return the probability of a Chi-Square with v degrees of freedom
%
%   AUTHORS: Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 9, June, 2015
%

    p = (1/2)^(1/2*v) *(x.^(1/2*v-1))./(gamma(1/2*v)) .* exp(-1/2*x);

end

