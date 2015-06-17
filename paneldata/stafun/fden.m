function [ p ] = fden( f, v1, v2 )
%FDEN Return the probability of a F
%   Return the probability of a F with v1 and v2 degrees of freedom
%
%   AUTHORS: Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 9, June, 2015
%

    p = (v1^(1/2*v1) * v2^(1/2*v2) * f.^(1/2*v1 - 1))./(beta(1/2*v1, 1/2*v2) * (v2+v1*f).^((v1+v2)/2));

end

