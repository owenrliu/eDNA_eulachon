# Eulachon allometry and eDNA
# Try to understand the difference between potential biomass and age structure
# weight-length relationship from Harvey 1987: W=0.00667L^3.12 (cm,gr)
# max length - ~25cm
# larvae - ~4-8mm = 0.4-0.8cm

library(tidyverse)

lw <- function(l){
  0.00667*l^(3.12)
}
wl <- function(w){
  (w/0.00667)^(1/3.12)
}

lwdat <- tibble(lcm=seq(0.4,25,length.out=100)) %>% 
  mutate(wg=lw(lcm))

glimpse(lwdat)
lwdat %>% 
  ggplot(aes(lcm,wg))+
  geom_line()+
  labs(x="Length (cm)",y="Weight (g)")

# eDNA shedding varies by life history stage/size
# We adopt the values directly from Maruyama et al. 2014, who studied juvenile and adult sunfish,
# which are actually similar length to eulachon. The release rate of juveniles (0.5-2g) was about
# 4 times that of adults (30-75g). What if we assume a release rate multiplier that scales down as
# the animal grows, from a max of 4 down to a minimum of 1

lwdat <- lwdat %>% 
  mutate(rscl=seq(4,1,length.out=100))
lwdat %>% 
  ggplot(aes(lcm,rscl))+
  geom_line()+
  labs(x="Length (cm)",y="eDNA release multiplier")
lwdat %>% 
  ggplot(aes(wg,rscl))+
  geom_line()+
  labs(x="Weight (g)",y="eDNA release multiplier")

# finally, if we pin this to some known eDNA release rate, we can get a sense for how different the final
# estimated biomass might be, depending on the age structure of our samples
# Again, following Maruyama et al., adults were shedding ~5*10^5 copies/g. use this to scale

lwdat <- lwdat %>% 
  mutate(copies_per_g=rscl*5e5) %>%
  mutate(total_copies=copies_per_g*wg)
lwdat %>% 
  ggplot(aes(wg,copies_per_g))+
  geom_line()+
  labs(x="Weight (g)",y="Copies per g")

lwdat %>% 
  ggplot(aes(wg,total_copies))+
  geom_line()+
  labs(x="Weight (g)",y="Copies")

copies_fxn <- function(jscl=4,minw,maxw,wg,sa=5e5){
  l <- wl(wg)
  minl <- wl(minw)
  maxl <- wl(maxw)
  scl=(1-(l-minl)/(maxl-minl))
  scl2 = (jscl-1)*scl+1
  cpg = scl2*sa
  tc = cpg*wg
  tc
  # scl2
}

x <- lwdat %>% 
  mutate(y=copies_fxn(wg=wg,minw=min(wg),maxw=max(wg)))

glimpse(x)
x %>% ggplot(aes(wg,y))+geom_line()

# allometric scaling of shedding?
shed <- function(a,b=0.7,wg){
  a*wg^b
}

