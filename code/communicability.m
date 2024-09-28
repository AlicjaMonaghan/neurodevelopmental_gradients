% create a function to calculate communicability
function communicability_matrix = communicability(a)
% 'a' is an nroi x nroi connectivity matrix
% find the strength matrix
s = diag(sum(a,2));
% normalise matrix by strength
pow = (s^-.5)*a*(s^-.5);
% take the normalised weighted communicability
c = expm(pow);
communicability_matrix = c;
end

