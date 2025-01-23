# Make ms figures for eulachon

# Figure 1 samples map and amplification
# samples with presence/absence (from simple_samples_map.R)
samples_pa_plot
numsamps_plot

f1 <-plot_grid(samples_pa_plot,numsamps_plot,ncol=2,rel_widths= c(1.25,1),labels='auto')
ggsave(here('plots','figure1.png'),f1,h=8,w=8)

# Figure 2 WA and OR pred maps (from main sdms script)
ppreds <-make_map_bathy(m.preds,est)+
  scale_fill_viridis(option="D")+
  labs(title="",fill="log(eDNA)",x="",y="")+
  ylim(wobbox[c(2,4)])+xlim(wobbox[c(1,3)])+
  coord_sf(datum=NA)+
  theme(legend.position = 'bottom')
ppreds
m.pa.map <- make_presence_absence_map(m)[[2]]+labs(x="",y="")+
  ylim(wobbox[c(2,4)])+xlim(wobbox[c(1,3)])+
  theme(legend.position = 'bottom')+
  coord_sf(datum=NA)
m.pa.map
f2 <-plot_grid(ppreds,m.pa.map,ncol=2,rel_widths= c(1,1),labels='auto')
ggsave(here('plots','figure2.png'),f2,h=8,w=8)

# Figure 3 index of abundance (from main sdms script)
m_abun
ggsave(here('plots','figure3.png'),m_abun,h=6,w=4)

# Figure 4 index of abundance and river discharge (from main sdms script)
lat_abun_river_disc_p
ggsave(here('plots','figure4.png'),lat_abun_river_disc_p,w=8,h=6,bg = "white")

# Figure 5 conditional covariate effects (from main sdms script)
mod1_covar_cond_effs
ggsave(here('plots','figure5.png'),mod1_covar_cond_effs,w=4,h=6,bg = "white")

## Supplementary Figures

## S1 cumulative samples
p_cumulative_samps
ggsave(here('plots','figureS1.png'),p_cumulative_samps,w=6,h=5,bg = "white")

## S2 and S3 covariates
thetao_plot <- make_map(grid.pred,thetao)+
  scale_fill_viridis(option="C")+
  labs(fill="Temperature")
ggsave(here('plots','figureS2.png'),thetao_plot,w=6,h=6,bg = "white")
k4_plot <- make_map(grid.pred,log(k4))+
  scale_fill_viridis(option="C",na.value="grey90")+
  labs(fill="Krill\nIndex")
ggsave(here('plots','figureS3.png'),k4_plot,w=6,h=6,bg = "white")

## S4-S7 Diagnostic Curves (from main sdms script)
diags <- make_pred_obs_plots(m,d,model_name="",saveplots = F)
pbin_st <- diags$p_bin_stand
ppos_st <- diags$p_pos_stand
pbin_unk <- diags$p_bin_unk
ppos_unk <- diags$p_pos_unk
ggsave(here('plots','figureS4.png'),pbin_st,w=6,h=5,bg = "white")
ggsave(here('plots','figureS5.png'),ppos_st,w=6,h=5,bg = "white")
ggsave(here('plots','figureS6.png'),pbin_unk,w=6,h=5,bg = "white")
ggsave(here('plots','figureS7.png'),ppos_unk,w=6,h=5,bg = "white")


## S8 and S9 INLA Mesh (from Test Mesh Complexity part of main sdms script)
s8 <- meshes.p+
  annotate("rect",xmin=0.16,xmax=0.20,ymin=0.32,ymax=0.39,col="red",fill=NA)+
  annotate("rect",xmin=0.71,xmax=0.75,ymin=0.41,ymax=0.48,col="red",fill=NA)
s8
ggsave(here('plots','figureS8.png'),s8,w=7,h=5,bg = "white")

mesh <- make_mesh(d_obs_filt,xy_cols=c('utm.lon.km','utm.lat.km'),type='cutoff',cutoff=40)
png(here('plots','figureS9.png'),w=480,h=800,bg="white")
plot(mesh)
dev.off()