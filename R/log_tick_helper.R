log_tick_helper <- function(min = -1, max = 4, base = 10) {
  stop_flag = F
  curr = min
  rtn <- numeric(0)
  while(!stop_flag) {
    nxt <- curr + 1
    rtn <- c(
      rtn,
      seq(base^curr, base^nxt-base^curr, by = base^curr)
    )
    curr <- nxt
    if (nxt > max) {
      stop_flag = T
    }
  }
  
  rtn <- rtn[base^max >= rtn]
  return(rtn)
}
