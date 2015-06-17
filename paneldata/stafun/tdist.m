function [ p ] = tdist( t, v )
%TDIST Return the cumulative distribution function of a student-t
%   Return the cumulative distribution function of a student-t with v
%   degrees of freedom
%
%   AUTHORS: Inmaculada C. Álvarez, Javier Barbero, José L. Zofío
%   http://www.paneldatatoolbox.com
%
%   Version: 2.0
%   LAST UPDATE: 9, June, 2015
%

    x = v./(v+t.^2);
    p(t <= 0) = 1/2 .* betainc(x(t<=0), 1/2*v, 1/2);
    p(t > 0) = 1 - 1/2 .* betainc(x(t>0), 1/2*v, 1/2);

end

