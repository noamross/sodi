match_times = function(sodi, times) {
  sodi_times = unique(sodi$time)
  matched_times = sodi_times[approx(x = sodi_times, y = 1:length(sodi_times),
                            xout = times, method = "constant", rule = 2)$y]
  return(matched_times)
}