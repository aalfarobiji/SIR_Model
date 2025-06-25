#'Wherever you see the word *MISSING*, you need to specify it!
# For mac users, you may not be able to run the ODE solver.
# In this case, please read the comment at the bottom of this file.

# Overview ----------------------------------------------------------------


# 1. Loading the necessary package
# 2. Set up the parameters
# 3. Set up the initial conditions
# 4. Define the model's rates of change
# 5. Setting the start and end time for the simulation
# 6. Run the model
# 7. Plot the output, and post-process the model output 
# 8. Explore! 


# 1. Loading package from library -----------------------------------------

# We will use the ode function in the deSolve package.
# If you have not installed deSolve, first you need to do that.
# Then you need to load this package from the library.
library(deSolve)


# 2. Set up parameters ----------------------------------------------------

 
# We will set the parameters for the model. 
# R0 is the basic reproduction number
# gamma is the recovery rate per day
# The average number of days of being infectious is 5, fill in the gamma parameter. 
parameters <- c(R0=1, gamma= 1/5)


# 3. Set up initial conditions --------------------------------------------

# We want to model a population size of 5.7 million, and will seed with one infected individual. 
# Everyone else is initially susceptible. Fill in the initial conditions and run this chunk of code.

initN <- 5700000 # population size
initI <- 1 # Infectious
initR <- 0 # Immune
initS <- initN-initI-initR # Susceptible 

# We need the next line of code to bring together these into initial conditions for the deSolve.

state <- c(S = initS, I = initI,R = initR)

rm(initN,initI,initS,initR) # Tidy up after using things for the last time

# 4. Create model function ------------------------------------------------

# Below is a function to run the COVID model. Fill in the dI equation. 
#' *Q1 On paper, draw the flow diagram for this model. Check if it differs from what we did in class.*

# Version 1 (the one usually found in examples online)
SIRmodel <- function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         
         # define variables
         N    <- S+I+R
         beta <- R0*gamma/N
         
         # rate of change
         dS <- -beta*I*S
         dI <- (beta*I*S)-(gamma*I)
         dR <- gamma*I 
         
         # return the rate of change
         list(c(dS, dI, dR))
       }
  ) 
}

# Version 2 which does the same thing
SIRmodel <- function(t, state, parameters)
{
  # This is an alternative function structure that some people prefer
  # because they do not like the "with" approach as it is conceptually
  # more complex. However this version requires more typing
  S = state['S']
  I = state['I']
  R = state['R']
  
  R0 = parameters['R0']
  gamma = parameters['gamma']
  
  # define variables
  N    <- S+I+R
  beta <- R0*gamma/N
         
  # rate of change
  dS <- -beta*I*S
  dI <- (beta*I*S)-(gamma*I)
  dR <- gamma*I 
         
  # return the rate of change
  return(list(c(dS, dI, dR)))
}


# 5. Specify time range ---------------------------------------------------

# We want to run the model for 6 months. Specify the start and end times that
# meet this expectation and then the output times we want from the model as below.

day_start <- 1
day_stop <- 365/2
times <- seq(day_start, day_stop)

rm(day_start,day_stop) # tidy tidy


# 6. Run the model --------------------------------------------------------



# Add the missing arguments to the ode function (which runs the model). 
# Results are stored in the "out" variable
# Hint: use ?deSolve if you want to get more information on the package
# and ?ode for more information on the function. 

out <- ode(y = state, times = times, func = SIRmodel,parms = parameters)

out2 = as.data.frame(out) # Some people prefer to have the output in a data frame
# so we will show alternatives that use this format

# 7. Make plots -----------------------------------------------------------


# The code below plots simple output. Run it to see the S, I and R compartment over time. 
plot(out)

# Alex alternative:
plot(out2$t,out2$S,type='l',xlab="Time (d)",ylab="Number of susceptibles",ylim=c(0,10))
plot(out2$t,out2$I,type='l',xlab="Time (d)",ylab="Number of infectives",ylim=c(0,10))
plot(out2$t,out2$R,type='l',xlab="Time (d)",ylab="Number removed",ylim=c(0,10))


# Often what we want to plot is the incidence of reported cases. 
# For this we need to have the reporting rate. 
# Fill in the reportprop parameter to represent that 25% of infections are reported. 

reportprop<-0.25
N <- out[,"S"]+out[,"I"]+out[,"R"]

# Alex alternative:
N = out2$S+out2$I+out2$R

# daily incidence
incidence <- reportprop*parameters["gamma"]*out[,"I"]
plot(incidence)

# Alex alternative
incidence <- reportprop*parameters["gamma"]*out2$I
plot(out2$t,incidence,type='l',xlab="Time (d)",ylab="Incidence of cases",ylim=c(0,1))


#' *Q2 Why are we plotting this and not reportprop times I?* 


# 8. Vary parameters ------------------------------------------------------


# You aren't sure whether the R0 is 3, 4 or 5. Run the model with all these values.
#' *Q3 How does the output differ?*  



















# For mac users:
# If you are having trouble installing the ODE solver package on a mac computer,
# it may be because you do not have a prerequisite software installed, namely a 
# fortran compiler. If this is indeed a problem for you, then to remedy this, 
# you should open up the terminal (probably 
# Launchpad -> Other [or Utilities?] -> Terminal 
# In the terminal, type the following to first download fortran to the current 
# directory and then to install it (you need to know the password for your 
# computer to do this):
# curl -OL http://r.research.att.com/libs/gfortran-4.8.2-darwin13.tar.bz2
# sudo tar fvxj gfortran-4.8.2-darwin13.tar.bz2 -C /
# Please note that this does actually involve downloading and then running
# "random" stuff off the internet, and no guarantees are offered by us that this
# will work or that it # will not destroy your computer. 
# Though it worked last year so Insha'Allah it will work this year too.




