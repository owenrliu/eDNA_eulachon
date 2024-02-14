# Utils for plotting output of sdmTMB models

make_pred_obs_plots <- function(modelobj,stands){
  tmbdat <- modelobj$tmb_data
  
  sdr <- modelobj$sd_report$value
  sdr <- tibble(param=names(sdr),val=sdr)
  # theta_stand(i) = std_xi_2(pcr_stand_bin_idx(i)) + std_xi_3(pcr_stand_bin_idx(i)) * (D_bin_stand(i));// - stand_offset);
  std_xi_2 <- sdr %>% filter(param=='std_xi_2') %>% pull(val)
  std_xi_3 <- sdr %>% filter(param=='std_xi_3') %>% pull(val)
  pcr_stand_bin_idx <- tmbdat$pcr_stand_bin_idx+1
  D_bin_stand <- tmbdat$D_bin_stand
  theta_stand <- std_xi_2[pcr_stand_bin_idx]+std_xi_3[pcr_stand_bin_idx]*D_bin_stand %>% as.numeric()
  theta_stand <- tibble(theta_stand=theta_stand) %>% 
    bind_cols(stands) %>% 
    mutate(Ct_bin=ifelse(Ct>0,1,0))
  
  #kappa_stand(i) = std_xi_0(pcr_stand_pos_idx(i)) + std_xi_1(pcr_stand_pos_idx(i)) * (D_pos_stand(i));// - stand_offset);
  std_xi_0 <- sdr %>% filter(param=='std_xi_0') %>% pull(val)
  std_xi_1 <- sdr %>% filter(param=='std_xi_1') %>% pull(val)
  pcr_stand_pos_idx <- tmbdat$pcr_stand_pos_idx+1
  D_pos_stand <- tmbdat$D_pos_stand
  kappa_stand <- std_xi_0[pcr_stand_pos_idx]+std_xi_1[pcr_stand_pos_idx]*D_pos_stand
  kappa_stand <- tibble(kappa_stand=kappa_stand,plate=modelobj$plates$plate[pcr_stand_pos_idx])
  
  # bring kappas and thetas together with plates
  pred.stand <- theta_stand %>% mutate(kappa_stand=NA)
  pred.stand$kappa_stand[pred.stand$Ct_bin==1] <- kappa_stand$kappa_stand

    #3 Pred- Obs Plots
  # pred.samp <- data.frame(theta_samp =report$theta_samp,
  #                         kappa_samp = report$kappa_samp)
  # pred.stand <- data.frame(theta_stand = report$theta_stand,
  #                          kappa_stand = report$kappa_stand)
  # pred.samp <- cbind(dat.samp,pred.samp)
  # pred.stand <- cbind(dat.stand,pred.stand)
  # 
  # pred.stand$sd_pred <- exp(report$ln_std_tau[1] + report$ln_std_tau[2] * pred.stand$ln_copies)
  # pred.stand <- pred.stand %>% mutate(kappa_plus = kappa_stand+sd_pred,
  #                                     kappa_minus = kappa_stand-sd_pred)
  
  ### MAKE SOME STANDARDS PLOTS
  STAND.BREAKS = c(1,5,10,100,1000,10000,100000)*2
  
  p_bin_stand <- ggplot(pred.stand) +
    geom_jitter(aes(y=Ct_bin,x=known_conc_ul),width=0,alpha=0.5,height=0.05) +
    geom_point(aes(y=plogis(theta_stand),x=known_conc_ul),color="red") +
    geom_line(aes(y=plogis(theta_stand),x=known_conc_ul),color="red") +
    scale_x_continuous(trans="log",breaks = STAND.BREAKS)+
    facet_wrap(~year) + 
    theme_bw()
  
  p_bin_stand2 <- ggplot(pred.stand) +
    geom_jitter(aes(y=Ct_bin,x=known_conc_ul),width=0,alpha=0.5,height=0.05) +
    geom_point(aes(y=plogis(theta_stand),x=known_conc_ul,color=plate,group=plate)) +
    geom_line(aes(y=plogis(theta_stand),x=known_conc_ul,color=plate,group=plate)) +
    scale_x_continuous(trans="log",breaks = STAND.BREAKS)+
    theme_bw() +
    theme(legend.position = "none")
  p_bin_stand2
  
  p_pos_stand <- ggplot(pred.stand %>% filter(Ct_bin ==1)) +
    geom_jitter(aes(y=Ct,x=known_conc_ul),width=0,alpha=0.5,height=0.05) +
    #geom_point(aes(y=kappa_stand,x=known_conc_ul),color="red") +
    geom_line(aes(y=kappa_stand,x=known_conc_ul,group=plate,color=plate)) +
    # geom_line(aes(y=kappa_plus,x=known_conc_ul),color="red",linetype="dashed") +
    # geom_line(aes(y=kappa_minus,x=known_conc_ul),color="red",linetype="dashed") +
    scale_x_continuous(trans="log",breaks = STAND.BREAKS)+
    guides(color='none')+
    facet_wrap(~year) +
    theme_bw()
  p_pos_stand
#   
#   p_pos_stand2 <- ggplot(pred.stand %>% filter(Ct_bin ==1)) +
#     geom_jitter(aes(y=Ct,x=known_conc_ul,group=year,color=year),width=0,alpha=0.5,height=0.05) +
#     #geom_point(aes(y=kappa_stand,x=known_conc_ul),color="red") +
#     geom_line(aes(y=kappa_stand,x=known_conc_ul,group=year,color=year)) +
#     #geom_line(aes(y=kappa_plus,x=known_conc_ul),color="red",linetype="dashed") +
#     #geom_line(aes(y=kappa_minus,x=known_conc_ul),color="red",linetype="dashed") +
#     scale_x_continuous(trans="log",breaks = STAND.BREAKS)+
#     facet_wrap(~year) + 
#     theme_bw() +
#     theme(legend.position = "none")
#   
#   # Unknown samples.
#   p_bin_unk <- ggplot(pred.samp) +
#     geom_jitter(aes(y=Ct_bin,x=plogis(theta_samp)),width=0,alpha=0.5,height=0.05) +
#     stat_smooth(aes(y=Ct_bin,x=plogis(theta_samp))) +
#     geom_abline(intercept=0,slope = 1,linetype="dashed",color="red") +
#     scale_x_continuous(limits=c(0,1)) +
#     facet_grid(year~depth_cat)  +
#     theme_bw()
#   
#   p_bin_unk2 <- ggplot(pred.samp) +
#     geom_jitter(aes(y=Ct_bin,x=plogis(theta_samp)),width=0,alpha=0.5,height=0.05) +
#     stat_smooth(aes(y=Ct_bin,x=plogis(theta_samp))) +
#     geom_abline(intercept=0,slope = 1,linetype="dashed",color="red") +
#     scale_x_continuous(limits=c(0,1)) +
#     facet_grid(depth_cat ~ year)  +
#     theme_bw()
#   
#   p_pos_unk2 <- ggplot(pred.samp %>% filter(Ct_bin ==1)) +
#     geom_jitter(aes(y=Ct,x=kappa_samp),width=0,alpha=0.5,height=0.05) +
#     stat_smooth(aes(y=Ct,x=kappa_samp)) +
#     geom_abline(intercept=0,slope = 1,linetype="dashed",color="red") +
#     facet_grid(depth_cat ~ year)  +
#     theme_bw()
#   #  scale_x_continuous(limits=c(0,1))
#   
#   
#   pdf(file = here("Plots and figures","Fits Standards and Pred-Obs.pdf"),  
#       onefile=T,width = 11,height=8.5)
#   print(p_bin_stand)
#   print(p_pos_stand)
#   print(p_bin_unk)
#   print(p_pos_unk)
#   dev.off()
#   
# }
