# Utils for plotting output of sdmTMB models
# see https://github.com/ericward-noaa/sdmTMB-stdcurve/blob/main/R/fit.R

library(contoureR)

plot_theme <- theme_minimal()+theme(panel.border = element_rect(color="black",fill=NA))
theme_black <- theme_minimal()+theme(panel.border = element_rect(color="black",fill=NA),text=element_text(color='white'))
theme_set(plot_theme)

make_pred_obs_plots <- function(modelobj,stands, model_name="",saveplots=T,savetype=c('pdf','png')){
  
  r <- modelobj$sd_report$par.random
  
  df <- data.frame(id = modelobj$plates$id,
                   plate = modelobj$plates$plate,
                   std_xi_0 = r[grep("std_xi_0", names(r))],
                   std_xi_1 = r[grep("std_xi_1", names(r))],
                   std_xi_2 = r[grep("std_xi_2", names(r))],
                   std_xi_3 = r[grep("std_xi_3", names(r))]
  )
  
  pred.samp <- predict(modelobj) %>% 
    # join hierarchical params for std curves
    left_join(df,by=join_by(plate)) %>% 
    mutate(kappa_samp=std_xi_0+std_xi_1*est,
           theta_samp=std_xi_2+std_xi_3*exp(est)) %>% 
    mutate(Ct_bin=ifelse(Ct>0,1,0))
  
  pred.stand <- stands %>% 
    mutate(theta_stand=modelobj$sd_report$value[grep("theta_stand", names(modelobj$sd_report$value))],
                              kappa_stand=NA) %>% 
    mutate(Ct_bin=ifelse(Ct>0,1,0))
  pred.stand$kappa_stand[which(stands$Ct > 0)] <- modelobj$sd_report$value[grep("kappa_stand", names(modelobj$sd_report$value))]
  
  
  
  # dplyr::filter(standards, Ct > 0) |>
  #   ggplot(aes(kappa_stand, Ct)) + geom_point(alpha=0.5) + 
  #   geom_abline(aes(intercept=0,slope=1),col="red") +
  #   theme_bw()
  
  ### MAKE SOME STANDARDS PLOTS
  STAND.BREAKS = c(1,5,10,100,1000,10000,100000) %>% as.integer()
  
  p_bin_stand <- ggplot(pred.stand) +
    geom_jitter(aes(y=Ct_bin,x=known_conc_ul),width=0,alpha=0.5,height=0.05) +
    geom_point(aes(y=plogis(theta_stand),x=known_conc_ul),color="red") +
    geom_line(aes(y=plogis(theta_stand),x=known_conc_ul,group=plate),color="red") +
    scale_x_continuous(trans="log10",breaks = STAND.BREAKS)+
    facet_wrap(~year) + 
    labs(x="Known Concentration (copies/uL)",y="Ct",title="Standards: Binary")
  
  p_pos_stand <- ggplot(pred.stand %>% filter(Ct_bin==1)) +
    geom_jitter(aes(y=Ct,x=known_conc_ul),width=0,alpha=0.5,height=0.05) +
    #geom_point(aes(y=kappa_stand,x=known_conc_ul),color="red") +
    geom_line(aes(y=kappa_stand,x=known_conc_ul,group=plate,color=plate)) +
    # geom_line(aes(y=kappa_plus,x=known_conc_ul),color="red",linetype="dashed") +
    # geom_line(aes(y=kappa_minus,x=known_conc_ul),color="red",linetype="dashed") +
    scale_x_continuous(trans="log10",breaks = STAND.BREAKS)+
    guides(color='none')+
    facet_wrap(~year) +
    labs(x="Known Concentration (copies/uL)",y="Ct",title="Standards: Positive")
  # p_pos_stand
  
  # Unknown samples
  p_bin_unk <- ggplot(pred.samp) +
    geom_jitter(aes(y=Ct_bin,x=plogis(theta_samp)),width=0,alpha=0.5,height=0.05) +
    stat_smooth(aes(y=Ct_bin,x=plogis(theta_samp))) +
    geom_abline(intercept=0,slope = 1,linetype="dashed",color="red") +
    scale_x_continuous(limits=c(0,1)) +
    facet_grid(year~depth_cat)  +
    labs(title="Unknown Samples: Pres/Abs",x="Predicted Probability of Amplification (theta)",
         y="Observed Ct>0")
  # p_bin_unk
  
#   
  p_pos_unk <- ggplot(pred.samp %>% filter(Ct_bin ==1)) +
    geom_jitter(aes(y=Ct,x=kappa_samp),width=0,alpha=0.5,height=0.05) +
    stat_smooth(aes(y=Ct,x=kappa_samp)) +
    geom_abline(intercept=0,slope = 1,linetype="dashed",color="red") +
    facet_grid(depth_cat ~ year)  +
    labs(title="Unknown Samples: Positive",x="Predicted Ct (kappa)",
         y="Observed Ct")
  
  if(saveplots){
    if(savetype=="png"){
      ggsave(here('plots',paste(model_name,"Standards Binary.png")),p_bin_stand,w=6,h=4.5)
      ggsave(here('plots',paste(model_name,"Standards Positive.png")),p_pos_stand,w=6,h=4.5)
      ggsave(here('plots',paste(model_name,"Unknowns Binary.png")),p_bin_unk,w=6,h=4.5)
      ggsave(here('plots',paste(model_name,"Unknown Pred vs Obs.png")),p_pos_unk,w=6,h=4.5)
      
    } else{
      pdf(file = here("plots",paste(model_name,"Fits Standards and Pred-Obs.pdf")),  
          onefile=T,width = 11,height=8.5)
      print(p_bin_stand)
      print(p_pos_stand)
      print(p_bin_unk)
      print(p_pos_unk)
      dev.off()    
      
    }
  }
  # p_pos_unk
  # p_bin_unk2
  allp<-list(
    p_bin_stand=p_bin_stand,
    p_pos_stand=p_pos_stand,
    p_bin_unk=p_bin_unk,
    p_pos_unk=p_pos_unk
  )
  return(allp)
}

# MAKE MAP OF PREDICTIONS

# coastline
coast <- ne_states(country='United States of America',returnclass = 'sf') %>%
  filter(name %in% c('California','Oregon','Washington','Nevada')) %>%
  st_transform(crs = pred.crs)

# here's the grid to predict with
grid.pred <- read_rds(here('data','prediction_grid_5km_sdmTMB_with_covars.rds')) %>% 
  mutate(across(c(thetao,so,bathy.bottom.depth,bottomT,k2,k4,k5),list(ln=function(x){
    x[x==0]<-1
    log(x)
  }))) %>% 
  # dummy covariate for 'washed'
  mutate(washed=0) %>% 
  # factors for year and depth (for some models that need these)
  mutate(depth_cat=as.factor(depth_cat),
         yr_fct=as.factor(year)) %>% 
  mutate(d1=ifelse(depth_cat==0,1,0),
         d2=ifelse(depth_cat==50,1,0),
         d3=ifelse(depth_cat==150,1,0))

# bathymetry
limits.for.map <- d_obs_filt %>% 
  st_as_sf(coords=c('lon','lat'),crs=4326) %>% 
  st_bbox()+c(-1,-1,1,1) #add an extra degree on each side for buffer


b = getNOAA.bathy(lon1 = limits.for.map["xmin"],
                  lon2 = limits.for.map[ "xmax" ],
                  lat1 = limits.for.map["ymin"],
                  lat2 = limits.for.map["ymax"],
                  resolution = 1,keep = TRUE)
# plot(b, image=TRUE, deep=-1000, shallow=0, step=250)
# make into a dataframe of contours for ggplotting
bdf <- fortify(b) %>% 
  contoureR::getContourLines(levels=c(-50, -150, -1000,-1500)) %>% 
  st_as_sf(coords=c("x","y"),crs=4326) %>% 
  st_transform(pred.crs) %>%
  mutate(x=st_coordinates(.)[,1],y=st_coordinates(.)[,2])

# bathy.only.map <- ggplot(bdf,aes(x,y,group=Group,colour=factor(z))) + geom_path()
# bathy.only.map+coord_equal()

# bounding box
bbox=st_bbox(grid.pred %>% st_as_sf(coords=c('x','y'),crs=pred.crs))
coastcrop = st_crop(coast,bbox)
bathycrop <- st_crop(bdf,bbox) %>% st_set_geometry(NULL)

make_map <- function(dat, column) {
  ggplot(coastcrop) + geom_sf()+
    geom_raster(data=dat, aes(x, y, fill = {{ column }})) +
    facet_grid(year~depth_cat)+
    theme_minimal()+
    labs(x="",y="")+
    theme(panel.border = element_rect(color='black',fill=NA),
          axis.text.x=element_blank())
}

make_map_bathy <- function(dat, column) {
  ggplot() + geom_sf(data=coastcrop)+
    geom_raster(data=dat, aes(x, y, fill = {{ column }})) +
    geom_path(data=bathycrop,aes(x,y,group=Group),color='gray30',linewidth=0.5)+
    facet_grid(year~depth_cat)+
    theme_minimal()+
    theme(panel.border = element_rect(color='black',fill=NA),
          axis.text.x=element_blank())
}

# MAKE PRESENCE-ABSENCE MAP- I THINK WE HAVE TO DO THIS MANUALLY

make_presence_absence_map <- function(modelobj,predgrid=grid.pred){
  x <- predict(modelobj,newdata=predgrid,nsim=30,model=1)
  y <- rowMeans(plogis(x))
  predgrid$pa <- y
  
  p1 <- ggplot(coastcrop) + geom_sf()+
    geom_raster(data=predgrid, aes(x, y, fill = pa)) +
    facet_grid(year~depth_cat)+
    scale_fill_viridis(option="C")+
    theme_minimal()+
    labs(fill="p(presence)")+
    theme(panel.border = element_rect(color='black',fill=NA),
          axis.text.x=element_blank())
  
  p2 <- ggplot(coastcrop) + geom_sf()+
    geom_raster(data=predgrid %>% 
                  mutate(top5=ifelse(pa>=0.95,T,F)), aes(x, y, fill = top5)) +
    facet_grid(year~depth_cat)+
    labs(fill="95% p(presence)")+
    scale_fill_manual(values=c('gray20','#BD3786FF'))+
    theme_minimal()+
    theme(panel.border = element_rect(color='black',fill=NA),
          axis.text.x=element_blank())
    
  return(list(p1,p2))
}
# MAKE CONDITIONAL PLOTS WITH PREDICT()
plot_marginal <- function(m,varname="thetao"){
  
  # conditional on means
  meanlat=mean(m$data$utm.lat.km,na.rm=T)
  meanlon=mean(m$data$utm.lon.km,na.rm=T)
  mean_thetao <- mean(m$data$thetao_ln,na.rm=T)
  mean_k4_ln <- mean(m$data$k4_ln,na.rm=T)
  #####
  
  newdata <- m$data %>% 
    mutate(varname=varname) %>% 
    mutate(thetao_ln=if_else(varname!="thetao",mean_thetao,thetao_ln)) %>% 
    mutate(k4_ln=if_else(varname!="k4",mean_k4_ln,k4_ln)) %>% 
    mutate(washed=1) %>% 
    mutate(offsets_all=mean(m$data$offsets_all)) %>% 
    mutate(yr_fct=as.factor(2019)) %>% 
    mutate(utm.lon.km=meanlon,utm.lat.km=meanlat)
  
  # run 100 sims from the joint precision matrix
  sims <- predict(m,newdata,return_tmb_object=F,nsim=100)
  mean_est <- apply(sims,1,mean)
  sd_est <- apply(sims,1,sd)
  
  # add sims to data to plot
  df_out <- newdata %>% 
    mutate(est=mean_est,upper=mean_est+sd_est,lower=mean_est-sd_est)
  
  df_out %>% 
    ggplot(aes(.data[[varname]],est,ymax=upper,ymin=lower))+
    geom_ribbon(color='gray50',fill='gray50',alpha=0.5)+geom_line(color='black',size=1.5)+
    labs(y="Marginal Effect")+
    # ylim(-50,NA)+
    theme_minimal()+
    theme(panel.background = element_rect(color='black'))
  
}
# make_cond_plot <- function(modelobj,v,exp_var=F,saveplot=T,return_dat=F){
#   d <- modelobj$data %>% 
#     # what if we used conditional means based on presence only...
#     filter(Ct>0)
#   var_range <- d %>% dplyr::select({{v}}) %>% pull(1) %>% range()
#   varseq <- seq(var_range[1],var_range[2],length.out=100)
#   
#   
#   dmeans <- d %>%
#     # what if we used conditional means based on presence only...
#     filter(Ct>0) %>%
#     summarise(across(where(is.numeric),mean)) %>% 
#     # and some factor levels (assume surface, 2019)
#     mutate(depth_cat=factor("0",levels=c("0","50","150")),
#            yr_fct=factor("2019",levels=c("2019","2021")))
#     
#   dummydat <- map_df(seq_len(length(varseq)),~dmeans) %>% mutate({{v}} := varseq)
#   pr <- predict(modelobj,newdata=dummydat,se_fit=TRUE,re_form=NA)
#   pr <- pr %>% mutate(
#     ymin = exp(est - 1.96 * est_se),
#     ymax = exp(est + 1.96 * est_se)
#   ) %>% 
#     dplyr::select({{v}},est,ymin,ymax)
#   
#   if(exp_var) pr <- pr %>% mutate({{v}} := exp({{v}}))
#   
#   vquo <- enquo(v)
#   
#   fit.name <- deparse(substitute(modelobj))
#   var.name <- deparse(substitute(v))
#   
#   plot.out <- pr %>% 
#     ggplot(aes(!!vquo,exp(est),ymin=ymin,ymax=ymax))+
#     geom_ribbon(alpha=0.4)+geom_line()+
#     theme_minimal()+
#     labs(x=var.name,y="eDNA copies",title=paste("Fit:",fit.name))+
#     theme(panel.border = element_rect(color='black',fill=NA))
#   
#   if(saveplot){
#     plot.name <- paste(fit.name,var.name,"conditional effect.png")
#     ggsave(here('plots',plot.name),plot.out,h=6,w=6)
#   }
#   if(return_dat) return(pr)
#   else return(plot.out)
# }

# SUM ABUNDANCE WITH PREDICTIONS

make_abun_index <- function(modelobj){
  preds <- predict(modelobj,newdata=grid.pred) %>% 
    mutate(expest=exp(est)) %>% 
    mutate(depth_cat=factor(depth_cat,levels=c("150","50","0"))) %>% 
    # make coastline segments
    mutate(state=case_when(
      lat_glorys>46.2 ~ "Washington (>46N)",
      lat_glorys<=46.2 & lat_glorys>42 ~ "Oregon (42-46N)",
      lat_glorys<=42 ~ "California (<42N)"
    )) %>% 
    mutate(state=factor(state,levels=c("Washington (>46N)","Oregon (42-46N)","California (<42N)")))
  abun_yr <- preds %>% 
    group_by(year,state,depth_cat) %>% 
    summarise(abun_idx=sum(expest))
  

  # abun_index_plot <- abun_yr %>% 
  #   ggplot(aes(year,abun_idx,fill=factor(depth_cat)))+
  #   geom_col(position='stack')+
  #   scale_fill_manual(values=viridis_pal(option="C")(3))+
  #   labs(x="Year",y="Index of Abundance\n(total eulachon eDNA)",fill="Depth Category")
  abun_index_plot2 <- abun_yr %>% 
    ggplot(aes(depth_cat,abun_idx,fill=factor(year)))+
    scale_fill_manual(values=c('gray20','#5DC863FF'))+
    geom_col(position='dodge')+
    # scale_y_reverse()+
    coord_flip()+
    facet_wrap(~state,nrow=3)+
    labs(x="Depth Category",y="Index of Abundance\n(total eulachon eDNA)",fill="Year")
  return(abun_index_plot2)
}

make_lat_abundance <- function(modelobj,depth_integrate=F){
  preds <- predict(modelobj,newdata=grid.pred) %>% 
    mutate(expest=exp(est))
  
  if(depth_integrate){
    abun_yr_lat <- preds %>% 
      group_by(year,lat_glorys) %>% 
      summarise(abun_idx=sum(expest)/1e05) # just to make plotting simpler
    abun_index_plot <- abun_yr_lat %>% 
      ggplot(aes(lat_glorys,abun_idx))+
      geom_line()+
      coord_flip()+
      facet_wrap(~year)+
      labs(y="Total eDNA Index",x="Latitude")+
      theme(text=element_text(size=14))
  } else{
    abun_yr_lat <- preds %>% 
      group_by(year,lat_glorys,depth_cat) %>% 
      summarise(abun_idx=sum(expest))
    abun_index_plot <- abun_yr_lat %>% 
      ggplot(aes(lat_glorys,abun_idx,color=factor(depth_cat)))+
      geom_line()+
      coord_flip()+
      scale_color_manual(values=PNWColors::pnw_palette("Starfish",3))+
      facet_wrap(~year)+
      labs(y="Total eDNA Index",x="Latitude",color="Depth")+
      theme(text=element_text(size=14))
  }
  
  return(abun_index_plot)
}
