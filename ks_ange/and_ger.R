# 1. Setting up workspace and data----------------------------------------------------------------------------------------

#Project ks_ange
#By Diana Spurite
#diana.spurite@posteo.de


library(MuMIn) 
#library( bblme ) # ver 1.0.25
library( ggplot2 ) # ver 3.4.4
library( htmlTable ) # ver 2.4.2
library( ipmr ) # ver 0.0.7
library( lme4 ) # ver 1.1-33
library( patchwork ) # ver 1.1.2
#library( pdbDigitUtils ) # ver 0.0.0.90
library( readxl ) # ver 1.4.2
library( tidyverse ) # ver 2.0.0
library( writexl ) # ver 1.4.2
library(boot)

and_ger            <- read.csv( "ks_ange/ks_ange.csv" )

#and_ger$logsize_t1 <- log( and_ger$size_tplus1 )

and_ger            <- and_ger %>% 
  mutate(logsize_t0 = log(basalArea_genet)) %>%
  mutate(logsize_t1 = log(size_tplus1)) 

surv               <- subset( and_ger, !is.na( survives_tplus1 ) ) %>%
  subset( basalArea_genet != 0 ) %>%
  select( Quad, Year, trackID,
          basalArea_genet, logsize_t0,logsize_t1,
          survives_tplus1, size_tplus1 )


grow              <- and_ger %>% 
  subset( basalArea_genet != 0) %>%
  subset( size_tplus1 != 0) %>% 
  select( Quad, Year, trackID,
          basalArea_genet, logsize_t0,logsize_t1,
          survives_tplus1, size_tplus1 )

#prep for recruitment

quad_df           <- and_ger %>%
  group_by( Species, Quad, Year ) %>%
  summarise( totPsize = sum( basalArea_genet ) ) %>%
  ungroup

group_df          <- quad_df %>%
  group_by( Species, Year ) %>%
  summarise( Gcov = mean( totPsize ) ) %>%
  ungroup

cover_df         <- left_join( quad_df, group_df ) %>%
  mutate( year = Year + 1 ) %>%
  mutate( year = as.integer( year ) ) %>%
  drop_na()
#recruitment

recr_df          <- and_ger %>%
  group_by( Species, Quad, Year ) %>%
  summarise( NRquad   = sum( recruit, na.rm=T ) ) %>%
  ungroup

recr             <- left_join( cover_df, recr_df ) %>%
  drop_na

write.csv( surv,  "ks_ange/data/survival_df.csv" )
write.csv( grow, "ks_ange/data/growth_df.csv" )
write.csv( recr, "ks_ange/data/recruitment_df.csv" )


# 2. Plotting data by year #---------------------------------------------


and_ger_long       <- pivot_longer( and_ger, cols = c( logsize_t0, logsize_t1 ), names_to = "size", values_to = "size_value" )

and_ger_long$size     <- as.factor( and_ger_long$size )
and_ger_long$Year_fac <- as.factor( and_ger_long$Year )

size_labs <- c( "at time t0", "at time t1" )
names(size_labs) <- c( "logsize_t0", "logsize_t1" )

print(head(and_ger_long))

plot <- and_ger_long %>% ggplot( aes( x = size_value ) ) +
  geom_histogram( binwidth = 1 ) +
  facet_wrap( ~ Year_fac,
              ncol = 4, 
            #  ~ size, 
              scales = "free_y",
              labeller = labeller( size = size_labs ) ) +
  labs( x = "log( size )",
        y = "Frequency" ) 

ggsave(filename = "ks_ange/results/histograms.png", plot = plot)
# given code " png( 'results/and_ger_yr/histograms.png', width = 6, height = 13, units = "in", res = 150 ) "
# didn't work, so I recycled code from another script


# 2. Survival data by year #---------------------------------------------

df_binned_prop_year <- function( ii, df_in, n_bins, siz_var, rsp_var, years ){
  
  # make sub-selection of data
  df   <- subset( df_in, Year == years$Year[ii] )
  
  if( nrow( df ) == 0 ) return( NULL )
  
  size_var <- deparse( substitute( siz_var ) )
  resp_var <- deparse( substitute( rsp_var ) )
  
  # binned survival probabilities
  h    <- ( max(df[,size_var], na.rm = T ) - min( df[,size_var], na.rm = T ) ) / n_bins
  lwr  <- min( df[,size_var], na.rm = T ) + ( h * c( 0:( n_bins - 1 ) ) )
  upr  <- lwr + h
  mid  <- lwr + ( 1/2 * h )
  
  binned_prop <- function( lwr_x, upr_x, response ){
    
    id  <- which( df[,size_var] > lwr_x & df[,size_var] < upr_x ) 
    tmp <- df[id,]
    
    if( response == 'prob' ){   return( sum( tmp[,resp_var], na.rm = T ) / nrow( tmp ) ) }
    if( response == 'n_size' ){ return( nrow( tmp ) ) }
    
  }
  
  y_binned <- Map( binned_prop, lwr, upr, 'prob' ) %>% unlist
  x_binned <- mid
  y_n_size <- Map( binned_prop, lwr, upr, 'n_size' ) %>% unlist
  
  # output data frame
  data.frame( xx  = x_binned, 
              yy  = y_binned,
              nn  = y_n_size ) %>% 
    setNames( c( size_var, resp_var, 'n_size' ) ) %>% 
    mutate( Year      = years$Year[ii] )
  
}

surv_yrs       <- data.frame( Year = surv$Year %>% unique %>% sort )
surv_bin_yrs   <- lapply( 1:39, df_binned_prop_year, and_ger, 20, 
                          logsize_t0, survives_tplus1, surv_yrs )

surv_yr_pan_df <- bind_rows( surv_bin_yrs ) %>% 
  mutate( transition = paste( paste0( Year ),
                              substr( paste0( Year + 1 ), 3, 4 ),
                              sep = '-' ) ) %>% 
  mutate( year       = as.integer( Year - 31  ) )


surv_yr_plot <- ggplot( data   = surv_yr_pan_df, 
                aes( x = logsize_t0, 
                     y = survives_tplus1 ) ) +
  geom_point( alpha = 0.5,
              pch   = 16,
              size  = 1,
              color = 'red' ) +
  scale_y_continuous( breaks = c( 0.1, 0.5, 0.9 ) ) +
  # split in panels
  facet_wrap( .~ transition, 
              nrow = 4 ) +
              theme_bw( ) +
              theme( axis.text = element_text( size = 8 ),
              title         = element_text( size = 10 ),
              strip.text.y  = element_text( size   = 5,
                                       margin = margin( 0.5, 0.5, 0.5, 0.5,
                                                        'mm' ) ),
         strip.text.x  = element_text( size   = 5,
                                       margin = margin( 0.5, 0.5, 0.5, 0.5,
                                                        'mm' ) ),
         strip.switch.pad.wrap = unit( '0.5', unit = 'mm' ),
         panel.spacing         = unit( '0.5', unit = 'mm' ) ) +
  labs( x = expression( 'log(size)'[t0] ),
        y = expression( 'Survival to time t1' ) )


# reoccurring warning message (in all the exercises):
# "Removed 95 rows containing missing values or values outside the scale range (`geom_point()`)."


ggsave(filename = "ks_ange/results/survival_binned_yr.png", plot = surv_yr_plot,                            # The plot object
       width = 10,                              
       height = 6,                              
       dpi = 300)                           



# 3. Growth data by year #---------------------------------------------

grow_yr_pan_df <- grow %>%
  mutate( transition = paste( paste0( Year ),
                              substr( paste0( Year + 1 ), 3, 4 ),
                              sep = '-' ) ) %>% 
  mutate( year       = as.integer( Year - 31 ) )

# png( 'results/and_ger_yr/growth_yr.png', width = 10, height = 6, units = "in", res = 150 )

plot <-  ggplot(data  = grow_yr_pan_df, aes( x = logsize_t0, y = log( size_tplus1 ) ) ) +
  geom_point( alpha = 0.5,
              pch   = 16,
              size  = 0.7,
              color = 'red' ) +
  # split in panels
  facet_wrap( .~ transition,
            nrow = 4 ) +
            theme_bw( ) +
            theme( axis.text     = element_text( size   = 8 ),
            title         = element_text( size   = 10 ),
            strip.text.y  = element_text( size   = 8,
                                       margin = margin( 0.5, 0.5, 0.5, 0.5,
                                                        'mm' ) ),
            strip.text.x  = element_text( size   = 8,
                                       margin = margin( 0.5, 0.5, 0.5, 0.5,
                                                        'mm' ) ),
         strip.switch.pad.wrap = unit( '0.5', unit = 'mm' ),
         panel.spacing         = unit( '0.5', unit = 'mm' ) ) +
  labs( x = expression( 'log( size )'[t0] ),
        y = expression( 'log( size )'[t1] ) )

ggsave(filename = "ks_ange/results/growth_yr.png", plot = plot)


# 4. Recruitment data by year #---------------------------------------------

indiv_qd <- surv %>%
  group_by( Quad ) %>%
  count( Year ) %>% 
  rename( n_adults = n ) %>% 
  mutate( Year = Year + 1 )

repr_yr <- indiv_qd %>% 
  left_join( recr ) %>%
  mutate( repr_pc    = NRquad / n_adults ) %>% 
  mutate( Year = Year - 1 ) %>% 
  drop_na

## Joining with `by = join_by(Quad, Year)`

# png( 'results/and_ger_yr/recruit_yr.png', width = 10, height = 6, units = "in", res = 150 )

plot <- repr_yr %>%
  filter( NRquad != max( repr_yr$NRquad ) ) %>% 
  filter( n_adults != max( repr_yr$n_adults ) ) %>% 
  ggplot( aes( x = n_adults, y = NRquad ) ) +
  geom_point( alpha = 1,
              pch   = 16,
              size  = 1,
              color = 'red' ) +
  facet_wrap( .~ Year, 
              nrow = 4 ) +
  theme_bw( ) +
  theme( axis.text     = element_text( size   = 8 ),
         title         = element_text( size   = 10 ),
         strip.text.y  = element_text( size   = 8,
                                       margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm' ) ),
         strip.text.x  = element_text( size   = 8,
                                       margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm' ) ),
         strip.switch.pad.wrap    = unit( '0.5', unit = 'mm' ),
         panel.spacing            = unit( '0.5', unit = 'mm' ) ) +
  labs( x = expression( 'Number of adults '[ t0] ),
        y = expression( 'Number of recruits '[ t1] ) )


ggsave(filename = "ks_ange/results/recruit_yr.png", plot = plot)


# 5. Recruitment histograms by year #---------------------------------------------

recSize <- and_ger %>% subset( recruit == 1)

recSize$Year_fac <- as.factor( recSize$Year )

# png( 'results/and_ger_yr/recr_histograms.png', width = 10, height = 6, units = "in", res = 150 )

plot <- recSize %>% ggplot( aes( x = logsize_t0 ) ) +
  geom_histogram( ) +
  facet_wrap( Year_fac ~ ., 
              scales = "free_y",
              ncol = 4 ) +
  labs( x = expression('log( size )'[t0]),
        y = "Frequency" )



ggsave(filename = "ks_ange/results/recruit_histograms.png", plot = plot)

# Check if there are any records for Year 33
table(recSize$Year)

# 6. Fitting vital rate models
##### Survival model####-----------------------

surv_df      <- surv %>% 
  mutate( logsize_t0 = log( basalArea_genet ) )

su_mod_yr <- glmer( survives_tplus1 ~ logsize_t0 + ( logsize_t0 | Year ), data = surv_df, family = binomial )

ranef_su <- data.frame( coef( su_mod_yr )[1] )
years_v  <- c( 32:70 )

surv_yr_plots <- function( i ){
  surv_temp   <- as.data.frame( surv_bin_yrs[[i]] )
  x_temp      <- seq( min( surv_temp$logsize_t0, na.rm = T ), 
                      max( surv_temp$logsize_t0, na.rm = T ), 
                      length.out = 100)
  pred_temp   <- boot::inv.logit( ranef_su[i,1] + ranef_su[i,2] * x_temp ) 
  pred_temp_df <- data.frame( logarea = x_temp, survives_tplus1 = pred_temp )
  temp_plot <- surv_temp %>% ggplot( ) +
    geom_point( aes( x = logsize_t0, y = survives_tplus1 ) ) +
    geom_line( data = pred_temp_df, aes( x     = logarea,
                                         y     = survives_tplus1 ),
               color = 'red',
               lwd   = 1  ) +
    labs( title = paste0( years_v[i] ),
          x = expression( 'log( size )'[t0] ),
          y = expression( 'Survival probability  '[ t1] ) )
  if( i %in% c( 33:41, 43:51, 53:61, 63:70 ) ){
    temp_plot <- temp_plot + theme( axis.title.y = element_blank( ) )
  }
  
  return(temp_plot)
}

surv_yrs <- lapply(1:39, surv_yr_plots)
surv_years <- wrap_plots(surv_yrs) + plot_layout(nrow = 4)

# png( 'results/Bou_gra_yr/survival_pred.png', width = 10, height = 8, units = "in", res = 150 )

surv_years


#ggsave(filename = "ks_ange/results/survival_pred.png", plot = surv_years)

ggsave(filename = "ks_ange/results/survival_pred.png", 
       plot = surv_years, 
       width = 30,
       height = 24,
       dpi = 300)  



#### Growth model####------------------------------------------------

grow_df      <- grow %>% 
  mutate( logsize_t0 = log( basalArea_genet ),
          logsize_t1 = log( size_tplus1 ) )

gr_mod_yr <- lmer( logsize_t1 ~ logsize_t0 + ( logsize_t0 | Year ), data = grow_df )


ranef_gr <- data.frame( coef( gr_mod_yr )[1] )

grow_yr_plots <- function( i ){
  temp_plot <- grow_df %>% filter( Year == i ) %>% ggplot( ) +
    geom_point( aes( x = logsize_t0, 
                     y = logsize_t1 ) ) +
    geom_abline( aes( intercept = ranef_gr[which(rownames( ranef_gr ) == i ),1],
                      slope     = ranef_gr[which(rownames( ranef_gr ) == i ),2] ),
                 color = "red",
                 lwd   = 1 ) +
    labs( title = paste0( i ),
          x = expression( 'log( size ) '[ t0] ),
          y = expression( 'log( size ) '[ t1] ) )
  if( i %in% c( 33:41, 43:51, 53:61, 61:70 ) ){
    temp_plot <- temp_plot + theme( axis.title.y = element_blank( ) )
  }
  
  return(temp_plot)
}
grow_yrs <- lapply( 32:70, grow_yr_plots )
grow_years <- wrap_plots( grow_yrs ) + plot_layout( nrow = 4 )

# png( 'results/Bou_gra_yr/grow_pred.png', width = 10, height = 6, units = "in", res = 150 )

grow_years

ggsave(filename = "ks_ange/results/grow_pred.png", plot = grow_years,
          width = 30,
         height = 20,
            dpi = 300) 

# 6. Quadratic & cubic term---------------------------------------------------------

grow_df$logsize_t0_2 <- grow_df$logsize_t0^2
grow_df$logsize_t0_3 <- grow_df$logsize_t0^3

gr_mod_yr2 <- lmer( logsize_t1 ~ logsize_t0 + logsize_t0_2 + ( logsize_t0 | Year ), data = grow_df )
gr_mod_yr3 <- lmer( logsize_t1 ~ logsize_t0 + logsize_t0_2 + logsize_t0_3 + ( logsize_t0 | Year ), data = grow_df )

g_mods <- c( gr_mod_yr, gr_mod_yr2, gr_mod_yr3 )

AICtab( g_mods, weights = T )

# Current R version does not allow this function


# 7. Quadratic term growth modeling--------------------------------------------------------------

ranef_gr2 <- data.frame( coef( gr_mod_yr2 )[1] )

grow_yr_plots2 <- function( i ){
  temp_f <- function( x ) ranef_gr2[which(rownames( ranef_gr2 ) == i ),1] + ranef_gr2[which(rownames( ranef_gr2 ) == i ),2] * x + ranef_gr2[which(rownames( ranef_gr2 ) == i ),3] * x^2 
  temp_plot <- grow_df %>% filter( Year == i ) %>% ggplot( ) +
    geom_point( aes( x = logsize_t0, 
                     y = logsize_t1 ) ) +
    geom_function( fun = temp_f,
                   color = "orange",
                   lwd   = 1 ) +
    labs( title = paste0( i ),
          x = expression( 'log( size ) '[ t0] ),
          y = expression( 'log( size ) '[ t1] ) )
  if( i %in% c( 33:41, 43:51, 53:61, 63:70 ) ){
    temp_plot <- temp_plot + theme( axis.title.y = element_blank( ) )
  }
  
  return(temp_plot)
}
grow_yrs2 <- lapply( 32:70, grow_yr_plots2 )
grow_years2 <- wrap_plots( grow_yrs2 ) + plot_layout( nrow = 4 )

# png( 'results/Bou_gra_yr/grow_pred_2.png', width = 10, height = 6, units = "in", res = 150 )

ggsave(filename = "ks_ange/results/grow_pred_quadratic.png", plot = grow_years2,
       width = 30,
       height = 20,
       dpi = 300) 

grow_years2


# 8. Cubic term growth modeling-------------------------------------------------

ranef_gr3 <- data.frame( coef( gr_mod_yr3 )[1] )

summary(gr_mod_yr3)

grow_yr_plots3 <- function( i ){
  temp_f <- function( x ) ranef_gr3[which(rownames( ranef_gr3 ) == i ),1] + ranef_gr3[which(rownames( ranef_gr3 ) == i ),2] * x + ranef_gr2[which(rownames( ranef_gr2 ) == i ),3] * x^2 
  temp_plot <- grow_df %>% filter( Year == i ) %>% ggplot( ) +
    geom_point( aes( x = logsize_t0, 
                     y = logsize_t1 ) ) +
    geom_function( fun = temp_f,
                   color = "turquoise",
                   lwd   = 1 ) +
    labs( title = paste0( i ),
          x = expression( 'log( size ) '[ t0] ),
          y = expression( 'log( size ) '[ t1] ) )
  if( i %in% c( 33:41, 43:51, 53:61, 63:70 ) ){
    temp_plot <- temp_plot + theme( axis.title.y = element_blank( ) )
  }
  
  return(temp_plot)
}
grow_yrs3 <- lapply( 32:70, grow_yr_plots3 )
grow_years3 <- wrap_plots( grow_yrs3 ) + plot_layout( nrow = 4 )

# png( 'results/Bou_gra_yr/grow_pred_3.png', width = 10, height = 6, units = "in", res = 150 )

ggsave(filename = "ks_ange/results/grow_pred_cubic.png", plot = grow_years3,
       width = 30,
       height = 20,
       dpi = 300) 

grow_years3

# 9. Growth variance model------------------------------------------------------------------

x <- fitted( gr_mod_yr2 )
y <- resid( gr_mod_yr2 )^2

gr_var <- nls( y ~ a * exp( b * x ), start = list( a = 1, b = 0 ) )


# 10. Recruitment model-----------------------------------------------------------------------

rec_mod <- glmer.nb( NRquad ~ ( 1 | Year ), data = recr )

recr_df        <- recr %>% 
  mutate( pred_mod      = predict( rec_mod, type = 'response' ) ) 

rec_sums_df <- recr_df %>% 
  group_by( Year ) %>% 
  summarise( NRquad    = sum( NRquad ),
             pred_mod  = sum( pred_mod ) ) %>% 
  ungroup


indiv_yr <- surv_df %>%
  count( Year ) %>% 
  rename( n_adults = n ) %>% 
  mutate( Year = Year + 1 )

repr_pc_yr <- indiv_yr %>% 
  left_join( rec_sums_df ) %>%
  mutate( repr_percapita = pred_mod / n_adults ) %>% 
  mutate( repr_pc_obs    = NRquad / n_adults ) %>% 
  mutate( Year = Year - 1 ) %>% 
  drop_na

recruit_pred <- repr_pc_yr %>% 
                ggplot() +
                geom_point( aes( x = repr_pc_obs,
                   y = repr_percapita ) ) +
                geom_abline( aes( intercept = 0,
                    slope     = 1 ),
                    color     = "red",
                    lwd       = 2,
                    alpha     = 0.5 ) +
               labs( x = "Observed per capita recruitment",
                     y = "Predicted per capita recruitment" )


## png( 'results/Bou_gra_yr/recruit_pred.png', width = 6, height = 4, units = "in", res = 150 )

ggsave(filename = "ks_ange/results/recruit_pred.png", plot = recruit_pred,
       width = 10,
       height = 8,
       dpi = 300) 


# 11. Exporting parameter estimates--------------------------------------------------

#survival####

su_yr_r <- data.frame( coefficient = paste0( "year_", rownames( coef( su_mod_yr )$Year ) ), 
                       value       = coef( su_mod_yr )$Year[,"(Intercept)"] )
su_la_r <- data.frame( coefficient = paste0( "logarea", rownames( coef( su_mod_yr )$Year ) ), 
                       value       = coef( su_mod_yr )$Year[,"logsize_t0"] )

surv_out_yr <- Reduce( function(...) rbind(...), list( su_la_r, su_yr_r ) ) %>%
  mutate( coefficient = as.character( coefficient ) )

write.csv( surv_out_yr, "ks_ange/data//surv_pars.csv", row.names = F )

#growth####


var_fe  <- data.frame( coefficient = names( coef( gr_var ) ),
                       value       = coef( gr_var ) )

year_re <- data.frame( coefficient = paste0( "year_", rownames( coef( gr_mod_yr2 )$Year ) ), 
                       value       = coef( gr_mod_yr2 )$Year[,"(Intercept)"] )

la_re   <- data.frame( coefficient = paste0( "logsize_t0", rownames( coef( gr_mod_yr2 )$Year ) ), 
                       value       = coef( gr_mod_yr2 )$Year[,"logsize_t0"] )

la2_re  <- data.frame( coefficient = paste0( "logsize_t0_2", rownames( coef( gr_mod_yr2 )$Year ) ), 
                       value       = coef( gr_mod_yr2 )$Year[,"logsize_t0_2"] )


grow_out_yr <- Reduce( function(...) rbind(...), list( var_fe, la_re, la2_re, year_re ) ) %>%
  mutate( coefficient = as.character( coefficient ) )

write.csv( grow_out_yr, "ks_ange/data/grow_pars.csv", row.names = F )

#recruitment####

rc_pc <- data.frame( coefficient = paste0( "rec_pc_", repr_pc_yr$Year ),
                     value = repr_pc_yr$repr_percapita )

rc_sz <- data.frame( coefficient = c( "rec_siz", "rec_sd" ),
                     value = c( mean( recSize$logsize_t0 ),
                                sd( recSize$logsize_t0 ) ) )

recr_out_yr <- Reduce( function(...) rbind(...), list( rc_pc, rc_sz ) ) %>%
  mutate( coefficient = as.character( coefficient ) )

write.csv( recr_out_yr, "ks_ange/data/recr_pars.csv", row.names = F )

#constant pars####

constants <- data.frame( coefficient = c( "recr_sz",
                                          "recr_sd",
                                          "a",
                                          "b",
                                          "L",
                                          "U",
                                          "mat_size" ),
                         value = c( mean( recSize$logsize_t0 ),  # are these four "logsize_t0" the same 
                                    sd( recSize$logsize_t0 ),    # are these four "logsize_t0" the same 
                                    as.numeric(coef(gr_var)[1]),
                                    as.numeric(coef(gr_var)[2]),
                                    grow_df$logsize_t0 %>% min,  # are these four "logsize_t0" the same 
                                    grow_df$logsize_t0 %>% max,  # are these four "logsize_t0" the same 
                                    200 ) )

surv_fe <- data.frame( coefficient = c( "surv_b0",
                                        "surv_b1" ),
                       value       = fixef( su_mod_yr ) )

grow_fe <- data.frame( coefficient = c( "grow_b0",
                                        "grow_b1",
                                        "grow_b2" ),
                       value       = fixef( gr_mod_yr2 ) )

rec_fe  <- data.frame( coefficient = "fecu_b0",
                       value       = mean( repr_pc_yr$repr_percapita ) )

pars_cons <- Reduce(function(...) rbind(...), list( surv_fe, grow_fe, rec_fe, constants ) ) %>%
  mutate(coefficient = as.character( coefficient))
          

rownames( pars_cons ) <- 1:13

pars_cons_wide <- as.list( pivot_wider( pars_cons, names_from = "coefficient", values_from = "value" ) )

write.csv( pars_cons_wide, "ks_ange/data/pars_cons.csv", row.names = F )

#varying pars####
su_b0 <- data.frame( coefficient = paste0( "surv_b0_", rownames( coef( su_mod_yr )$Year ) ), 
                       value       = coef( su_mod_yr )$Year[,"(Intercept)"] )
su_b1 <- data.frame( coefficient = paste0( "surv_b1_", rownames( coef( su_mod_yr )$Year ) ), 
                       value       = coef( su_mod_yr )$Year[,"logsize_t0"] )
grow_b0 <- data.frame( coefficient = paste0( "grow_b0_", rownames( coef( gr_mod_yr2 )$Year ) ), 
                       value       = coef( gr_mod_yr2 )$Year[,"(Intercept)"] )
grow_b1   <- data.frame( coefficient = paste0( "grow_b1_", rownames( coef( gr_mod_yr2 )$Year ) ), 
                       value       = coef( gr_mod_yr2 )$Year[,"logsize_t0"] )
grow_b2   <- data.frame( coefficient = paste0( "grow_b2_", rownames( coef( gr_mod_yr2 )$Year ) ), 
                       value       = coef( gr_mod_yr2 )$Year[,"logsize_t0_2"] )
fecu_b0 <- data.frame( coefficient = paste0( "fecu_b0_", repr_pc_yr$Year ),
                     value = repr_pc_yr$repr_percapita )

pars_var <- Reduce(function(...) rbind(...), list( su_b0, su_b1, grow_b0, grow_b1, grow_b2, fecu_b0 ) )

pars_var_wide <- as.list( pivot_wider( pars_var, names_from = "coefficient", values_from = "value" ) )

write.csv( pars_var_wide, "ks_ange/data/pars_var.csv", row.names = F )

