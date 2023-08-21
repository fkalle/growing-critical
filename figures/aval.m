########################################################################################
# Kalle Kossio et al. 2018, Phys. Rev. Lett. 121, 058301,
# https://doi.org/10.1103/PhysRevLett.121.058301
########################################################################################

function [N, D] = aval(spike_times, bin_size)

  # usage: [N, D] = aval(spike_times, bin_size)
  # This function accepts list of event times
  # and size of time bin, it returns sizes of 
  # detected avalanches and their durations

  bin_num = fix((spike_times(end)-spike_times(1))/bin_size) + 1;

  # Make bin edges, starting from first spike to last
  edges = linspace(spike_times(1), spike_times(1)+bin_num*bin_size, bin_num + 1);
  n = histc(spike_times, edges); # Bin the spike times

  aval_edges = [edges(1) edges(n == 0)]; # Edges separating avalanches
   
  N = histc(spike_times, aval_edges);
  N(end) = [];
  N(N==0) = [];

  # Calculation of avalanche durations
  end_element = cumsum(N);
  begin_element = (end_element - N) .+ 1;
  D = spike_times(end_element) - spike_times(begin_element);
  
endfunction
