########################################################################################
# Octave script to replicate Fig. 2 of Kalle Kossio et al. 2018, 
# Phys. Rev. Lett. 121, 058301, https://doi.org/10.1103/PhysRevLett.121.058301
########################################################################################

FNAME_SP = "../data/nc/spiketime";		# File with spike times [milliseconds]
FNAME_N = "../data/nc/neuron";			# File with neuron indices
FNAME_R = "../data/nc/radii";			# File with neuron radii
FNAME_O = "../data/nc/overlap";			# File with neuron total overlap
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
N = fread(bfile, 1, 'int'); # Number of neurons
start_time = fread(bfile, 1, 'double');
net = fread(bfile, 11*N, 'double');
fclose(bfile);
net = reshape(net, 11, []);
xs = net(2,:); 
ys = net(3,:);
g = net(10,1); # Coupling proportionality constant [1/milliseconds]
tau = net(8,1); # Time constant [milliseconds]

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

bfile = fopen(FNAME_O);
overlap = fread(bfile, 'double');
fclose(bfile);

overlap = transpose(reshape(overlap, N+1, []));

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

average_total_overlap = mean(overlap(:,2:end),2);

interval = [0:1e6:2e8 (2e8+1e7):1e7:1e9];
bin_ov = [];
bin_av_ov =[];
for i=1:1:length(interval)-1
  ind = overlap(:,1) > interval(i) & overlap(:,1) < interval(i+1);
  if(sum(ind) > 10)
    bin_ov = [bin_ov ; mean(overlap(ind,:)) ];
    bin_av_ov = [bin_av_ov; mean(average_total_overlap(ind))];
  else 
    bin_ov = [bin_ov ; overlap(ind,:) ];
    bin_av_ov = [bin_av_ov; average_total_overlap(ind)];
  endif
endfor
  
hold off;

for i=1:1:NPLOT
  plot(bin_ov(:,1)/1e3, g*tau*bin_ov(:,i),
    'color', [190,193,212]./255,
    'linewidth', 1
     );
  hold on;
endfor

plot(bin_ov(:,1)/1e3, g*tau*bin_av_ov,'k', 'linewidth', 2);
ylabel("g\\tau A");
xlabel("Time");

axis([0.0, 1e6,0.0, 1.55]);

for i=1:1:3
  plot([TIMEPOINTS(i),TIMEPOINTS(i)], [0,1.55], 'color', COLS(i,:));
endfor

set(gca(), 
  'linewidth', 1,
  'tickdir', 'out',
  'ytick', [0, 0.5,1,1.5]
);

hold off;
box off;

print('figure2.eps', 
  '-S300,300',
  '-color',
  '-F:Helvetica:7',
  '-tight'
);
