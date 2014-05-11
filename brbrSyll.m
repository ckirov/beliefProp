%Berber Syllabification

%Inputs:
%phons: a list of phonological segment sonority levels. Each sonority may be
%valued between 1 and 8. No two consecutive segments may have the same
%sonority level.

%Outputs:
%out: a list indicating which positions should be treated as syllable
%nuclei

function out = brbrSyll(phons)


%string
slen = length(phons);

%edges
E = [1:slen-1;2:slen];
E = [E [E(2,:);E(1,:)]]';
elen = size(E,1);

%univariates
U = zeros(slen,2);
U(:,2) = (2.^phons-1)';
U = exp(U);

%bivariates
B = zeros(elen,2,2);
B(:,2,2) = (-2^8);
B = exp(B);

%run beliefprop in maximization mode
Prbs = beliefProp(E,U,B,2);

%read out most probable states
out = Prbs(:,2) > Prbs(:,1);
out = out';

