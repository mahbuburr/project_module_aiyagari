% Package on Inequality Metrics such as the Gini Coefficient associated to the Lorenz Curve, the Theil and the Atkinson Indexes.
%
% For more references:
% http://en.wikipedia.org/wiki/Lorenz_curve
% http://en.wikipedia.org/wiki/Gini_coefficient
% http://en.wikipedia.org/wiki/Theil_index
% http://en.wikipedia.org/wiki/Atkinson_index

%-*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-*%
%                                                                                               %
%            Author: Liber Eleutherios                                             %
%            E-Mail: libereleutherios@gmail.com                             %
%            Date: 19 May 2008                                                     %
%                                                                                               %
%-*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-* -*-*%

%-*-* -*-* -*-* -*-* -*-* -*-* -*-*

% The zip file contains:
% ginicoeff.m
% lorenzcurve.m
% theilt.m
% theill.m
% atkinsonineq.m

%-*-* -*-* -*-* -*-* -*-* -*-* -*-*

% Example for ginicoeff.m:
N = 1000;
p = rand(N, 1); w = rand(N, 1);
y = ginicoeff(p, w)

%-*-* -*-* -*-* -*-* -*-* -*-* -*-*

% Example for lorenzcurve.m:
N = 1000;
p = rand(N, 1);
w = exp(randn(N, 1));
h = lorenzcurve(p, w);

%-*-* -*-* -*-* -*-* -*-* -*-* -*-*

% Example for theilt.m:
N = 1000;
w = rand(N, 1);
y = theilt(w)

%-*-* -*-* -*-* -*-* -*-* -*-* -*-*

% Example for theill.m:
N = 1000;
w = rand(N, 1);
y = theill(w)

%-*-* -*-* -*-* -*-* -*-* -*-* -*-*

% Example for atkinsonineq.m:
N = 1000;
w = rand(N, 1);
epsilon = rand;
y = atkinsonineq(w, epsilon)
