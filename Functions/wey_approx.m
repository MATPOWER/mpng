function [f,df,d2f] = wey_approx(K,pi_i,pi_j,PI_STAR)
% WEY_APPROX Computes the gas flow function and its first and second order derivatives 
%   for a set of pipelines according to the Weymouth's equation.
%   [F,DF,D2F] = WEY_APPROX(K,PI,PJ,P_APPROX)
% 
%   This function also calculates a non-linear approximation of the Weymouth's 
%   equation within an interval around zero. Inputs are as follow:
%   
%   K - Weymouth constant vector for all pipelines.
% 
%   PI - Quadratic pressure vector at node from for all pipelines.
% 
%   PJ - Quadratic pressure vector at node for of all pipelines.
%   
%   PI_STAR (optional) - Vector with the values of PI* for each pipeline
%                               
%
%
%   For more details on the derivations behind the derivative code used
%   in MPNG information, see the user's manual.
%
%   See also MPNG
%
%   Main author: Wilson González Vanegas M.Sc

%   MPNG: MATPOWER - Natural Gas
%   Copyright (c) 2019-2022 - v0.99beta
%   Sergio García-Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González-Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo-Sánchez - Universidad Nacional de Colombia - Sede Manizales
% 
%   This file is part of MPNG.
%   Covered by the 3-clause BSD License (see LICENSE file for details).

%% 

%% default margins for the weymouth polynomial approximation

if nargin < 4                               % 'Do not' approximate
    PI_STAR = 1e-15*ones(length(K),1);      %  Set a very small approx. interval to simulate the 'no approx. status'
end
%% inizialization
f = zeros(numel(K),1);
df = zeros(numel(K),1);
if nargout < 4 
    d2f = zeros(numel(K),1);
end

%% finding bounds for the approximation
p = pi_i-pi_j;
id_l = p < -PI_STAR;
id_u = p > PI_STAR;
id_lu = p >= -PI_STAR & p <= PI_STAR;
%% calculation outside the bounds
f(id_l) = -sqrt(-K(id_l).*p(id_l));
df(id_l) = K(id_l)./(2*sqrt(-K(id_l).*p(id_l)));
if nargout < 4
    d2f(id_l) = (K(id_l).^2)./(4*sqrt((-K(id_l).*p(id_l)).^3));
end

f(id_u) = sqrt(K(id_u).*p(id_u));
df(id_u) = K(id_u)./(2*sqrt(K(id_u).*p(id_u)));
if nargout < 4
    d2f(id_u) = -(K(id_u).^2)./(4*sqrt((K(id_u).*p(id_u)).^3));
end

%% calculation of the approximation inside the bounds
if any(id_lu)
    k1 = sqrt(K(id_lu).*PI_STAR(id_lu));
    k2 = K(id_lu)./(2*sqrt(K(id_lu).*PI_STAR(id_lu)));
    k3 = -((K(id_lu).^2)./(4*(sqrt((K(id_lu).*PI_STAR(id_lu)).^3))));
    c = (k1 - PI_STAR(id_lu).*k2 + (k3/3).*(PI_STAR(id_lu).^2))./((8/3)*(PI_STAR(id_lu).^5));
    b = (k3 - (20*c).*(PI_STAR(id_lu).^3))./(6*PI_STAR(id_lu));
    a = k2 - 0.5*PI_STAR(id_lu).*(k3 - (20*c).*(PI_STAR(id_lu).^3)) - (5*c).*(PI_STAR(id_lu).^4);

    f(id_lu)    = a.*p(id_lu) + b.*(p(id_lu).^3) + c.*(p(id_lu).^5);
    df(id_lu)   = a + (3*b).*(p(id_lu).^2) + (5*c).*(p(id_lu).^4);
    if nargout < 4
        d2f(id_lu)  = (6*b).*p(id_lu) + (20*c).*(p(id_lu).^3);
    end
end
