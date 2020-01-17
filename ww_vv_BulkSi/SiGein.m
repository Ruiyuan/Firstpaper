function [eps,isigma,leps,epssig, imass] = SiGein(lambdaSi, lambdaGe, epsSi, epsGe, sigmaSi, sigmaGe, massSi, massGe)
% ! S-W potential parameters for mixtures
eps(1,1) = epsSi;
eps(1,2) = 0.5*(epsSi+epsGe);
eps(2,1) = 0.5*(epsSi+epsGe);
eps(2,2) = epsGe;

isigma(1,1) = 1.0/sigmaSi;
isigma(1,2) = 2.0/(sigmaSi+sigmaGe);
isigma(2,1) = 2.0/(sigmaSi+sigmaGe);
isigma(2,2) = 1.0/sigmaGe;

imass(1,1) = 1.0; %/massSi;
imass(1,2) = 1.0; %2.0/(massSi+massGe);
imass(2,1) = 1.0;% 2.0/(massGe+massSi);
imass(2,2) = 1.0; %/massGe;

leps(1,1,1) = lambdaSi*epsSi;
leps(1,1,2) = 2.0*lambdaSi*epsSi/3.0 + lambdaGe*epsGe/3.0;
leps(1,2,1) = 2.0*lambdaSi*epsSi/3.0 + lambdaGe*epsGe/3.0;
leps(2,1,1) = 2.0*lambdaSi*epsSi/3.0 + lambdaGe*epsGe/3.0;
leps(1,2,2) = lambdaSi*epsSi/3.0 + 2.0*lambdaGe*epsGe/3.0;
leps(2,1,2) = lambdaSi*epsSi/3.0 + 2.0*lambdaGe*epsGe/3.0;
leps(2,2,1) = lambdaSi*epsSi/3.0 + 2.0*lambdaGe*epsGe/3.0;
leps(2,2,2) = lambdaGe*epsGe;

epssig(1,1) = epsSi/sigmaSi;
epssig(1,2) = (epsSi+epsGe)/(sigmaSi+sigmaGe);
epssig(2,1) = (epsSi+epsGe)/(sigmaSi+sigmaGe);
epssig(2,2) = epsGe/sigmaGe;

end