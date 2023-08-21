########################################################################################
# Kalle Kossio et al. 2018, Phys. Rev. Lett. 121, 058301,
# https://doi.org/10.1103/PhysRevLett.121.058301
########################################################################################

function [t,n,n0] = durhist(durations, bin)

## Given list of avalanche durations (e.g. from outpoot of aval() function)
## and bin size, returns normalized histogram, t is time and n is fraction of
## durations within interval (t, t+bin]


ind = durations > 0;
nonzero_dur_fraction = sum(ind)/length(durations);
t = [0:bin:max(durations(ind))];
m = histc(durations(ind), t);
n = nonzero_dur_fraction * m ./ trapz(t,m');
t(n==0) = [];
n(n==0) = [];

n0 = 1-nonzero_dur_fraction;

endfunction
