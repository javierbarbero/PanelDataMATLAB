function [ p ] = tden( t, v )
%TDEN Return the probability of a student-t function
%   Return the probability of a student-t function with v degrees of
%   freedom
%
%   AUTHORS: Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 9, June, 2015
%

    p = gamma(1/2*(v+1))/(gamma(1/2*v)*sqrt(v*pi)) * (1 + 1/v*t.^2).^(-1/2*(v+1));

end

