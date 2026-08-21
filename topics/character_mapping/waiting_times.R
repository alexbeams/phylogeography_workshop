

Q <- matrix(c(-1,10,1,-10), nrow=2, byrow=T)


events <- 0
states <- 1

Tmax <- 60

tm <- 0
state <- 1

while(tm < Tmax){
    tm <- tm + rexp(1, rate= -Q[state, state])
    events = c(events, tm)

    state <- if(state==1){2}else{1}
    states <- c(states, state )
}


dat <- cbind(events, states)

