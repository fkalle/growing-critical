########################################################################################
# Octave script to replicate Fig. 3 of Kalle Kossio et al. 2018, 
# Phys. Rev. Lett. 121, 058301, https://doi.org/10.1103/PhysRevLett.121.058301
########################################################################################

# Near-critical (NC)
FNAME_NC = "../data/nc/spiketime";	# File with spike times [milliseconds]
TBIN_NC = 45;				# Time bin for detecting avalanches [milliseconds]
SIGMA_NC = 0.995;			# Branching parameter
TAU_NC = 10; 				# Time constant [milliseconds]

# Subcritical (SC)
FNAME_SC = "../data/sc/spiketime";	# File contining spike times [milliseconds]
TBIN_SC = 30;				# Time bin for detecting avalanches [milliseconds]
SIGMA_SC = 0.75; 			# Branching parameter
TAU_SC = 10;		 		# Time constant [milliseconds]

DBIN = 0.2;				# Bin size for durations, both NC and SC [milliseconds]
STARTTIME = 0.9e9;			# From what time to do analysis [milliseconds], allows network to equilibrate
	
graphics_toolkit gnuplot;

########################################################################################

spiketime_file = fopen(FNAME_NC);
spiketime = fread(spiketime_file, 'double');
fclose(spiketime_file);
spiketime = spiketime(spiketime>STARTTIME);

[N_NC,D_NC] = aval(spiketime, TBIN_NC);

spiketime_file = fopen(FNAME_SC);
spiketime = fread(spiketime_file, 'double');
fclose(spiketime_file);
spiketime = spiketime(spiketime>STARTTIME);

[N_SC,D_SC] = aval(spiketime, TBIN_SC);

########################################################################################
subplot(1,2,1);

s1 = unique(N_NC);
Ps1 = histc(N_NC,s1);
Ps1 = Ps1/sum(Ps1);

s2 = unique(N_SC);
Ps2 = histc(N_SC,s2);
Ps2 = Ps2/sum(Ps2);

# Stirling's approximation Eq. (5)
As1 = s1.^(-3/2).*exp((-SIGMA_NC + log(SIGMA_NC) +1)*s1)/(sqrt(2*pi())*SIGMA_NC);
As2 = s2.^(-3/2).*exp((-SIGMA_SC + log(SIGMA_SC) +1)*s2)/(sqrt(2*pi())*SIGMA_SC);

loglog(
   s2,Ps2,
  'linestyle', 'none',
  'color', [128,128,128]./256,
  'markersize', 2,
  'marker', 's',
   s2, As2, 
  'color', [133,149,225]./256,
  'linewidth', 1.5,
   s1,Ps1,
  'linestyle', 'none',
  'color','k',
  'markersize', 2,
  'marker', 's',
   s1, As1, 
  'color', [211,63,106]./256,
  'linewidth', 1.5
   );
   
axis([1e+0 10^4.3 10^(-5.5) 1e+0]);
axis('square');

ylabel("Probability");
xlabel("Avalanche size");

set(gca(), 
  'linewidth', 2,
  'tickdir', 'out'
);

box off;

set(gcf,
  'PaperUnits', 'normalized',
  'PaperSize', [1,1],
  'PaperPosition', [ 0, 0, 1, 1]
);

########################################################################################
subplot(1,2,2);

# Binning durations
[t, n, n0] = durhist(D_NC, DBIN);

# Analitycal solution Eq. (10)
adot = @(a,t) (-a/TAU_NC + exp(SIGMA_NC*a/TAU_NC) -1);
jac = @(a,t) (-1 + SIGMA_NC*exp(SIGMA_NC*a/TAU_NC))/TAU_NC;

T = 0: 0.1: 10;
T = [T 11:1:2000];
A = lsode({adot,jac}, [-TAU_NC,0], T); A = A(:,1);
dur = SIGMA_NC*adot(A,T).*exp(SIGMA_NC*A/TAU_NC)/TAU_NC;

# Approximation Eq. (12)
a = 2*TAU_NC^2./(-2*TAU_NC - T);
dadt = a.^2/(2*TAU_NC^2);
f = exp(a/TAU_NC).*dadt/TAU_NC;

loglog(
   t,n,
  'linestyle', 'none',
  'marker', 's',
  'color', 'k',
  'markersize', 2,
   T,f,
  'linewidth', 1.5,
  'color', [211,63,106]./256,
   T, dur,
  'linewidth', 1.5,
  'color', [239,151,8]./256
   );
      
hold on;

# Binning durations
[t, n, n0] = durhist(D_SC, DBIN);

# Analitycal solution Eq. (10)
adot = @(a,t) (-a/TAU_SC + exp(SIGMA_SC*a/TAU_SC) -1);
jac = @(a,t) (-1 + SIGMA_SC*exp(SIGMA_SC*a/TAU_SC))/TAU_SC;

T = 0: 0.1: 10;
T = [T 11:1:2000];
A = lsode({adot,jac}, [-TAU_SC,0], T); A = A(:,1);
dur = SIGMA_SC*adot(A,T).*exp(SIGMA_SC*A/TAU_SC)/TAU_SC;

loglog(
   t,n, 
  'linestyle', 'none', 
  'marker', 's', 
  'color', [128,128,128]./256,
  'markersize', 2,
   T, dur, 
  'linewidth', 1.5, 
  'color', [17,198,56]./256
   );
      
hold off;

xlabel("Duration (ms)");

axis([1e-1 1e+3 1e-4 0.02]);
axis('square');

set(gca(),
  'linewidth', 1,
  'tickdir', 'out',
  'ticklength', [0.02, 0.02]
);

box off;

set(gcf,
  'PaperUnits', 'normalized',
  'PaperSize', [1,1],
  'PaperPosition', [ 0, 0, 1, 1]
);

print(
  'figure3.eps', 
  '-S600,600',
  '-color',
  '-F:Helvetica:14',
  '-tight'
);
