########################################################################################
# Octave script to replicate Fig. 2 of Kalle Kossio et al. 2018, 
# Phys. Rev. Lett. 121, 058301, https://doi.org/10.1103/PhysRevLett.121.058301
########################################################################################

FNAME_SP = "../data/nc/spiketime";		# File with spike times [milliseconds]
FNAME_N = "../data/nc/neuron";			# File with neuron indices
FNAME_R = "../data/nc/radii";			# File with neuron radii
FNAME_INIT = "../data/nc/initstate";		# File with initial state

TIMEPOINTS = [0, 0.8e5, 8e5];			# Time points to plot [seconds]

C1 = [0.55686 0.02352 0.23137];
C2 = [0.00784 0.24705 0.64705];			# Colors
C3 = [0.06666 0.77647 0.21960];
COLS = [C1; C2; C3];

# Panel (a)
BOXLIMITS = [-0.21, 1.25 , -0.21, 1.25];	# Size of the box to plot [unitless]

# Panel (b)
NPLOT = 25;					# Number of neurons to plot
TPLOT = 100.0; 					# Duration of time interval to plot [seconds]

graphics_toolkit gnuplot;

########################################################################################

function h = drawcircle(x, y, radius, npoints)
	theta=linspace(0,2*pi,npoints);
	rho=ones(1,npoints)*radius;
	[xs,ys] = pol2cart(theta,rho);
	xs=xs+x;
	ys=ys+y;
	h=plot(xs,ys);
endfunction

########################################################################################
bfile = fopen(FNAME_INIT);
N = fread(bfile, 1, 'int'); #Number of neurons
start_time = fread(bfile, 1, 'double');
net = fread(bfile, 11*N, 'double');
fclose(bfile);
net = reshape(net, 11, []);
xs = net(2,:);
ys = net(3,:);

bfile = fopen(FNAME_R);
radii = fread(bfile, 'double');
fclose(bfile);

radii = reshape(radii, N+1, []);
radii_times = radii(1,:)./1e3; # Convert from ms to s
radii = radii(2:end,:);

bfile = fopen(FNAME_SP);
spiketime = fread(bfile, 'double');
fclose(bfile);
spiketime = spiketime./1e3; # Convert from ms to s

bfile = fopen(FNAME_N);
neuron = fread(bfile, 'int');
fclose(bfile);
neuron = neuron.+1; # Start neuron index from 1 instead of 0

########################################################################################
for i=1:3
  time = TIMEPOINTS(i);
  ind = (radii_times > time & radii_times < time + TPLOT);
  rs = radii(:,ind);
  
  subplot(3, 3, i);
  for n=1:N
    c=drawcircle(xs(n), ys(n), rs(n,1), min(round(rs(n,1)*1000), 100));
    set(c, 'linewidth', 1, 'color', COLS(i,:));
    hold on;
  endfor
  hold off;

  axis(BOXLIMITS);
  axis('square');
  axis('equal');
  a=gca;
  set(a, 'box','off', 'color','none');
  axis('off');
  b=axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
  set(b, 'linewidth', 1)
  axis('square');
  axes(a);  
endfor

########################################################################################

for i=1:3
  subplot(3, 3, i+3);
  time = TIMEPOINTS(i); # Start time of the plot

  ind=(neuron < NPLOT & spiketime > time  & spiketime < time  + TPLOT);
  plot(spiketime(ind), neuron(ind),
    'linestyle', 'none',
    'linewidth', 0.5,
    'color', COLS(i,:),
    'markersize', 1,
    'marker','s'
     );
  axis([time, time+TPLOT, 0, NPLOT+1]);
  pbaspect([1 0.5 1]);

  set(gca()
    ,'linewidth', 1
    ,'tickdir', 'out'
    ,'ytick', [1, NPLOT]
    ,'xtick', [time]
    );
endfor

########################################################################################

subplot(3,3,7:9)

print('figure2.eps', 
  '-S300,300',
  '-color',
  '-F:Helvetica:7',
  '-tight'
);
