# 1. Setting up workspace and data----------------------------------------------------------------------------------------

#Project ks_bogr
#By Diana Spurite
#diana.spurite@posteo.de

library( bbmle ) # ver 1.0.25
library( ggplot2 ) # ver 3.4.4
library( htmlTable ) # ver 2.4.2
library( ipmr ) # ver 0.0.7
library( lme4 ) # ver 1.1-33
library(remotes)
library( patchwork ) # ver 1.1.2
#library( pdbDigitUtils ) # ver 0.0.0.90
library( readxl ) # ver 1.4.2
library( tidyverse ) # ver 2.0.0
library( writexl ) # ver 1.4.2
library(boot)

bou_gra  <- read.csv( "ks_bogr/ks_bogr.csv" )

#bou_gra$logsize_t1 <- log( bou_gra$size_tplus1 )

bou_gra <- bou_gra %>% 
           mutate( logsize_t0 = log(basalArea_genet)) %>%
           mutate( logsize_t1 = log(size_tplus1)) 

surv    <- subset( bou_gra, !is.na( survives_tplus1 ) ) %>%
           subset( basalArea_genet != 0 ) %>%
           select( Quad, Year, trackID,
                   basalArea_genet,
                   logsize_t0,logsize_t1,
                   survives_tplus1, size_tplus1 )

grow    <- bou_gra %>% 
           subset( basalArea_genet != 0) %>%
           subset( size_tplus1 != 0) %>% 
           select( Quad, Year, trackID,
                   basalArea_genet,
                   logsize_t0,logsize_t1,
                   survives_tplus1, size_tplus1 )

#prep for recruitment

quad_df  <- bou_gra %>%
            group_by( Species, Quad, Year ) %>%
            summarise( totPsize = sum( basalArea_genet ) ) %>%
            ungroup

group_df <- quad_df %>%
            group_by( Species, Year ) %>%
            summarise( Gcov = mean( totPsize ) ) %>%
            ungroup

cover_df <- left_join( quad_df, group_df ) %>%
            mutate( year = Year + 1 ) %>%
            mutate( year = as.integer( year ) ) %>%
            drop_na()

#recruitment

recr_df <- bou_gra %>%
           group_by( Species, Quad, Year ) %>%
           summarise( NRquad = sum( 
           recruit, na.rm=T ) ) %>%
           ungroup

recr    <- left_join( cover_df, recr_df ) %>%
           drop_na


write.csv ( surv, "ks_bogr/data/survival_df.csv" )
write.csv ( grow, "ks_bogr/data/growth_df.csv" )
write.csv ( recr, "ks_bogr/data/recruitment_df.csv" )


# 2. Plotting data by year #---------------------------------------------


bou_gra_long          <- pivot_longer( bou_gra, cols = 
                         c( logsize_t0, logsize_t1 ), 
                         names_to = "size", 
                         values_to = "size_value" )

bou_gra_long$size     <- as.factor( bou_gra_long$size )
bou_gra_long$Year_fac <- as.factor( bou_gra_long$Year )

size_labs             <- c( "at time t0", "at time t1" )
names(size_labs)      <- c( "logsize_t0", "logsize_t1" )

print(head(bou_gra_long))

plot <- bou_gra_long %>% 
        ggplot( aes( x = size_value ) ) +
        geom_histogram( binwidth = 1 ) +
        facet_wrap( ~ Year_fac, 
                  ncol = 4, 
                scales = "free_y",
              labeller = labeller( 
                  size = size_labs ) ) +
               labs( x = "log( size )",
                     y = "Frequency" ) + theme( 
          strip.text.y = element_text( 
                  size = 8,
                margin = margin( 0.5, 0.5, 0.5, 0.5,'mm' ) ),
          strip.text.x = element_text( 
                  size = 8,
                margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm' ) ),
 strip.switch.pad.wrap = unit( '0.5', unit = 'mm' ),
         panel.spacing = unit( '0.5', unit = 'mm' ) ) +
                         theme_bw()


ggsave(filename = "ks_bogr/results/histograms.png",
           plot = plot,
          width = 6.3,
         height = 9,
            dpi = 300)


# 2. Survival data by year #---------------------------------------------

df_binned_prop_year <- function( ii, df_in, n_bins, siz_var, rsp_var, years ){ 
    
                 df <- subset( df_in, Year == years$Year[ii] ) # make sub-selection of data
                       if( nrow( df ) == 0 ) return( NULL )
                 
           size_var <- deparse( substitute( siz_var ) )
           resp_var <- deparse( substitute( rsp_var ) )  # binned survival probabilities 
           
                  h <- ( max ( df [,size_var], na.rm = T ) - 
                         min ( df [,size_var], na.rm = T ) ) / n_bins
                  
                lwr <-   min ( df [,size_var],  na.rm = T ) + ( h * c( 0:( n_bins - 1 ) ) )
                upr <- lwr + h
                mid <- lwr + ( 1/2 * h )
                
        binned_prop <- function( lwr_x, upr_x, response ){
                 id <- which( df[,size_var] > lwr_x & df[,size_var] < upr_x ) 
                 
                tmp <- df[id,] 
                       if( response == 'prob' )  { 
                           return( sum( tmp [,resp_var], na.rm = T ) / nrow( tmp ) ) }
                       if( response == 'n_size' ){ 
                           return( nrow( tmp ) ) }
  }
           y_binned <- Map( binned_prop, lwr, upr, 'prob' ) %>% unlist
           x_binned <- mid
           y_n_size <- Map( binned_prop, lwr, upr, 'n_size' ) %>% unlist
                       data.frame( xx  = x_binned,
                                   yy  = y_binned,
                                   nn  = y_n_size ) %>%
                      setNames( c ( size_var, resp_var, 'n_size' ) ) %>% 
                       mutate( Year = years$Year[ii] ) 
}

surv_yrs       <- data.frame( Year = surv$Year %>% unique %>% sort )

surv_bin_yrs   <- lapply( 1:38, df_binned_prop_year, bou_gra, 20, 
                          logsize_t0, survives_tplus1, surv_yrs )

surv_yr_pan_df <- bind_rows( surv_bin_yrs ) %>% 
                  mutate( transition = 
                  paste ( paste0( Year ),
                  substr( paste0( Year + 1 ), 3, 4 ), 
                  sep = '-' ) ) %>% 
                  mutate( year = as.integer( Year - 31  ) )

surv_yr_plot   <- ggplot( data = surv_yr_pan_df, 
                       aes ( x = logsize_t0, 
                             y = survives_tplus1 ) ) +
                  geom_point(alpha = 0.5,
                           pch = 16,
                          size = 1,
                         color = 'red' ) +
                  scale_y_continuous(
                        breaks = c( 0.1, 0.5, 0.9 ) ) +
                  facet_wrap( .~ transition, 
                          ncol = 4 ) + 
                  theme_bw( ) + theme(   
                     axis.text = element_text( size = 8 ),
                         title = element_text( size = 10 ),
                  strip.text.y = element_text( size = 5,
                        margin = margin( 0.5, 0.5, 0.5, 0.5,  'mm' ) ),
                  strip.text.x = element_text( size = 5,
                        margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm' ) ),
         strip.switch.pad.wrap = unit ( '0.5', unit = 'mm' ),
                 panel.spacing = unit ( '0.5', unit = 'mm' ) ) +
                       labs( x = expression( 'log(size)'[t0] ),
                             y = expression( 'Survival to time t1' ) )


ggsave(filename = "ks_bogr/results/survival_binned_yr.png",           
           plot = surv_yr_plot,
          width = 6.3,                              
         height = 9,                              
            dpi = 300)                         


# 3. Growth data by year #---------------------------------------------

grow_yr_pan_df <- grow %>%
                  mutate( transition = 
                  paste ( paste0( Year ),
                  substr( paste0( Year + 1 ), 3, 4 ), 
                          sep = '-' ) ) %>% 
                  mutate(year = as.integer( Year - 31 ) )

plot          <-  ggplot(data = grow_yr_pan_df, 
                        aes(x = logsize_t0, 
                            y = log( size_tplus1 ) ) ) +
                  geom_point( alpha = 0.5,
                          pch = 16,
                         size = 0.7,
                        color = 'red' ) +
                  facet_wrap( .~ transition, # split in panels
                         ncol = 4 ) + 
                  theme_bw( ) + theme( 
                    axis.text = element_text( size = 8 ),
                        title = element_text( size = 10 ),
                 strip.text.y = element_text( size = 8,
                       margin = margin( 0.5, 0.5, 0.5, 0.5,'mm' ) ),
                 strip.text.x = element_text( size = 8,
                       margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm' ) ),
        strip.switch.pad.wrap = unit( '0.5', unit = 'mm' ),
                panel.spacing = unit( '0.5', unit = 'mm' ) ) +
                      labs( x = expression( 'log( size )'[t0] ),
                            y = expression( 'log( size )'[t1] ) )


ggsave(filename = "ks_bogr/results/growth_yr.png",
           plot = plot,
         width  = 6.3,                              
         height = 9,                              
            dpi = 300 )


# 4. Recruitment data by year #---------------------------------------------

indiv_qd <- surv %>%
            group_by( Quad ) %>%
            count ( Year ) %>% 
            rename( n_adults = n ) %>% 
            mutate( Year = Year + 1 )

repr_yr  <- indiv_qd %>% 
            left_join( recr ) %>%
            mutate( repr_pc = NRquad / n_adults ) %>% 
            mutate( Year = Year - 1 ) %>% 
            drop_na

plot     <- repr_yr %>%
            filter( NRquad  != max( repr_yr$NRquad ) ) %>% 
            filter(n_adults != max( repr_yr$n_adults ) ) %>% 
            ggplot( aes ( x  = n_adults, y = NRquad ) ) +
            geom_point(alpha = 1,
                         pch = 16,
                        size = 1,
                       color = 'red' ) + 
            facet_wrap(.~ Year, 
                        ncol = 4 )     + 
            theme_bw( ) + theme(   
                   axis.text = element_text( size = 8 ),
                       title = element_text( size = 10 ),
                strip.text.y = element_text( size = 8,
                      margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm' ) ),
                strip.text.x = element_text( size = 8,
                      margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm' ) ),
       strip.switch.pad.wrap = unit ( '0.5', unit = 'mm' ),
               panel.spacing = unit ( '0.5', unit = 'mm' ) ) +
                     labs( x = expression( 'Number of adults '[ t0] ),
                           y = expression( 'Number of recruits '[ t1] ) )


ggsave(filename = "ks_bogr/results/recruit_yr.png",
           plot = plot,
          width = 6.3,
         height = 9, 
            dpi = 300)

# 5. Recruitment histograms by year #---------------------------------------------

recSize          <- bou_gra %>% subset( recruit == 1)

recSize$Year_fac <- as.factor( recSize$Year )

plot            <- recSize %>% 
                   ggplot( aes( x = logsize_t0 ) ) +
                   geom_histogram( ) +
                   facet_wrap( Year_fac ~ ., 
                   scales = "free_y",
                     ncol = 4 ) +
                   labs(x = expression('log( size )'[t0]),
                        y = "Frequency" )

table(recSize$Year) # no record for the years 32 and 33


ggsave(filename = "ks_bogr/results/recruit_histograms.png", 
           plot = plot,
          width = 6.3,
         height = 9, 
            dpi = 300)

# 6. Fitting vital rate models
##### * Survival model####-----------------------

surv_df        <- surv %>% mutate( logsize_t0 = log( basalArea_genet ) )

su_mod_yr      <- glmer( survives_tplus1 ~ logsize_t0 + ( logsize_t0 | Year ), 
                         data = surv_df, family = binomial )

ranef_su       <- data.frame( coef( su_mod_yr )[1] )
years_v        <- c( 34:71 )

surv_yr_plots  <- function( i ){
  surv_temp    <- as.data.frame( surv_bin_yrs[[i]] )
  x_temp       <- seq( min( surv_temp$logsize_t0, na.rm = T ), 
                       max( surv_temp$logsize_t0, na.rm = T ), 
                                             length.out = 100)
  
  pred_temp     <- boot::inv.logit( ranef_su[i,1] + ranef_su[i,2] * x_temp ) 
  pred_temp_df  <- data.frame( logarea = x_temp, survives_tplus1 = pred_temp )
  
  temp_plot     <- surv_temp %>% 
                   ggplot( ) +
                   geom_point(aes( x = logsize_t0, 
                                   y = survives_tplus1 ) ) +
                   geom_line (  data = pred_temp_df, 
                              aes( x = logarea,
                                   y = survives_tplus1 ),
                               color = 'red',
                                 lwd = 1  ) +
                   theme_bw() + theme( 
                           axis.text = element_text( size = 5 ),
                               title = element_text( size = 5),
                          plot.title = element_text( size = 10,
                               hjust = 0.5,
                              margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm' ) ),
                         plot.margin = margin( 0,0,0,0) ) +
                         labs( title = paste0( years_v[i] ),
                                   x = expression( 'log( size )'[t0] ),
                                   y = expression( 'Survival probability'[ t1] ) )
                                       if( i %in% c(1:38)[!(c(1:38) %% 4) == 1] ){
                                         
     temp_plot <- temp_plot + theme( 
                  axis.title.y = element_blank( ) ) }
                  return(temp_plot)
}

surv_yrs   <- lapply(1:38, surv_yr_plots)
surv_years <- wrap_plots(surv_yrs) + plot_layout(ncol = 4)


ggsave(filename = "ks_bogr/results/survival_pred.png", 
           plot = surv_years, 
          width = 6.3,
         height = 9,
            dpi = 300)    


#### Growth model####------------------------------------------------

grow_df   <- grow %>% 
             mutate( logsize_t0 = log( basalArea_genet ),
                     logsize_t1 = log( size_tplus1 ) )

gr_mod_yr <- lmer(logsize_t1 ~ logsize_t0 + (logsize_t0 | Year), data = grow_df)

ranef_gr  <- data.frame( coef( gr_mod_yr )[1] )

grow_yr_plots <- function( i ){
    temp_plot <- grow_df %>% 
                 filter( Year == i ) %>% 
                 ggplot( ) +
                 geom_point(
                 aes( x = logsize_t0, 
                      y = logsize_t1 ) ) +
                 geom_abline( aes( 
              intercept = ranef_gr[which(rownames( ranef_gr ) == i ),1],
                  slope = ranef_gr[which(rownames( ranef_gr ) == i ),2] ),
                  color = "red",
                    lwd = 1 ) +
                 theme_bw() + theme( 
              axis.text = element_text( size = 5 ),
                  title = element_text( size = 8),
             plot.title = element_text( size = 10,
                  hjust = 0.5,
                 margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm' ) ),
            plot.margin = margin(0,0,0,0,'mm') ) +
            labs( title = paste0( i ),
                      x = expression( 'log( size ) '[ t0 ] ),
                      y = expression( 'log( size ) '[ t1 ] ) )
                          if( i %in% c(34:71)[!(c(34:71) %% 4) == 0] ){
                            
    temp_plot <- temp_plot + 
                 theme( axis.title.y = element_blank( ) )  }
                 return(temp_plot)
}

grow_yrs    <- lapply( 34:71, grow_yr_plots )
grow_years  <- wrap_plots( grow_yrs ) + plot_layout( ncol = 4 )


ggsave(filename = "ks_bogr/results/grow_pred.png", 
           plot = grow_years,
          width = 6.3,
         height = 9,
            dpi = 300) 


# 6. Quadratic & cubic term---------------------------------------------------------

#grow_df    <- grow_df %>% filter(Year>33)

grow_df$logsize_t0_2 <- grow_df$logsize_t0^2
grow_df$logsize_t0_3 <- grow_df$logsize_t0^3

gr_mod_yr2 <- lmer( logsize_t1 ~ logsize_t0 + logsize_t0_2 + 
                  ( logsize_t0 | Year ), data = grow_df )

gr_mod_yr3 <- lmer( logsize_t1 ~ logsize_t0 + logsize_t0_2 + 
                    logsize_t0_3 + ( logsize_t0 | Year ), 
                    data = grow_df )

table(grow_df$Year)

g_mods     <- c( gr_mod_yr, 
                 gr_mod_yr2, 
                 gr_mod_yr3 )

AICtab( g_mods, weights = T )

#dAIC   df weight
#model3    0.0 8  0.9945
#model2   10.4 7  0.0055
#model1 1018.5 6  <0.001

# 7. Quadratic term growth modeling--------------------

ranef_gr2      <- data.frame( coef( gr_mod_yr2 )[1] )

grow_yr_plots2 <- function( i ){
        temp_f <- function( x ) 
                  ranef_gr2[which(rownames( ranef_gr2 ) == i ),1] + 
                  ranef_gr2[which(rownames( ranef_gr2 ) == i ),2] * x + 
                  ranef_gr2[which(rownames( ranef_gr2 ) == i ),3] * x^2 
     temp_plot <- grow_df %>% filter( Year == i ) %>% ggplot( ) +
                  geom_point( aes( x = logsize_t0, 
                                   y = logsize_t1 ) ) +
                  geom_function( fun = temp_f,
                               color = "orange",
                               lwd   = 1 ) +
                  theme_bw() + theme( 
                           axis.text = element_text( size = 3 ),
                               title = element_text( size = 5),
                          plot.title = element_text( size = 10,
                               hjust = 0.5,
                              margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm' ) ),
                         plot.margin = margin(0,0,0,0) ) +
                         labs( title = paste0( i ),
                                   x = expression( 'log( size ) '[ t0] ),
                                   y = expression( 'log( size ) '[ t1] ) )
                                       if( i %in% c(34:71)[!(c(34:71) %% 4) == 0] ){
   temp_plot <- temp_plot + 
                theme ( axis.title.y = 
                element_blank( ) )}
                return(temp_plot)
}

grow_yrs2    <- lapply( 34:71, grow_yr_plots2 )
grow_years2  <- wrap_plots( grow_yrs2 ) + 
                plot_layout( ncol = 4 )


ggsave(filename = "ks_bogr/results/grow_pred_quadratic.png", 
           plot = grow_years2,
          width = 6.3,
         height = 9,
            dpi = 300)


# 8. Cubic term growth modeling-------------------------------------------------

ranef_gr3      <- data.frame( coef( gr_mod_yr3 )[1] )

summary(gr_mod_yr3)

grow_yr_plots3 <- function( i ){
        temp_f <- function( x ) 
                  ranef_gr3[which(rownames( ranef_gr3 ) == i ),1] + 
                  ranef_gr3[which(rownames( ranef_gr3 ) == i ),2] * x + 
                  ranef_gr2[which(rownames( ranef_gr2 ) == i ),3] * x^2 
     temp_plot <- grow_df %>% 
                  filter( Year == i ) %>% 
                  ggplot( ) +
                  geom_point( aes( x = logsize_t0, 
                                   y = logsize_t1 ) ) +
                  geom_function( fun = temp_f,
                               color = "turquoise",
                                 lwd = 1 ) +
                  theme_bw( ) + theme( 
                          axis.text  = element_text( size = 3 ),
                               title = element_text( size = 5),
                          plot.title = element_text( size = 10,
                               hjust = 0.5,
                              margin = margin( 0.5, 0.5, 0.5, 0.5, 'mm' ) ),
                         plot.margin = margin(0,0,0,0) ) +
                         labs( title = paste0( i ),
                                   x = expression( 'log( size ) '[ t0] ),
                                   y = expression( 'log( size ) '[ t1] ) )
                                       if( i %in% c(34:71)[!(c(34:71) %% 4) == 0] ){
                                        
    temp_plot <- temp_plot + 
                 theme( axis.title.y = 
                 element_blank( ) ) }
                 return(temp_plot)
}


grow_yrs3   <- lapply( 34:71, grow_yr_plots3 )
grow_years3 <- wrap_plots( grow_yrs3 ) +
               plot_layout( ncol = 4 )


ggsave(filename = "ks_bogr/results/grow_pred_cubic.png" ,
           plot = grow_years3,
          width = 6.3,
         height = 9,
            dpi = 300) 


# 9. Growth variance model------------------------------------------------------------------

x <- fitted( gr_mod_yr2 )
y <- resid( gr_mod_yr2 )^2

gr_var <- nls( y ~ a * exp( b * x ), start = list( a = 1, b = 0 ) )

# 10. Recruitment model-----------------------------------------------------------------------

rec_mod       <- glmer.nb( NRquad ~ ( 1 | Year ), data = recr )

recr_df       <- recr %>% 
                 mutate ( pred_mod = predict( rec_mod, 
                              type = 'response' ) ) 

rec_sums_df   <- recr_df %>% 
                 group_by ( Year ) %>% 
                 summarise( NRquad = sum( NRquad ),
                          pred_mod = sum( pred_mod ) ) %>% 
                          ungroup  


indiv_yr     <- surv_df %>%
                count( Year ) %>% 
                rename( n_adults = n ) %>% 
                mutate( Year = Year + 1 )

repr_pc_yr   <- indiv_yr %>% 
                left_join( rec_sums_df ) %>%
                mutate( repr_percapita = pred_mod / n_adults ) %>% 
                mutate( repr_pc_obs = NRquad / n_adults ) %>% 
                mutate( Year = Year - 1 ) %>% 
                drop_na

recruit_pred <- repr_pc_yr %>% 
                ggplot() +
                geom_point( aes( x = repr_pc_obs,
                                 y = repr_percapita ) ) +
                geom_abline( aes( intercept = 0,
                                      slope = 1 ),
                                      color = "red",
                                        lwd = 2,
                                      alpha = 0.5 ) +
                                    labs( x = "Observed per capita recruitment",
                                          y = "Predicted per capita recruitment" )


ggsave(filename = "ks_bogr/results/recruit_pred.png", 
           plot = recruit_pred,
          width = 6.3,
         height = 9,
            dpi = 300) 


# 11. Exporting parameter estimates--------------------------------------------------

# * survival####

su_yr_r   <- data.frame( coefficient = paste0( "year", 
             rownames  ( coef ( su_mod_yr )$Year ) ), 
                 value = coef ( su_mod_yr )$Year[,"(Intercept)"] )

su_la_r   <- data.frame( coefficient = paste0( "logsize_t0", 
             rownames  ( coef ( su_mod_yr )$Year ) ), 
                 value = coef ( su_mod_yr )$Year[,"logsize_t0"] )

surv_out_yr <- Reduce ( function(...) 
               rbind  (...), 
               list ( su_la_r, 
                      su_yr_r ) ) %>%
               mutate ( coefficient = as.character( coefficient ) )


write.csv( surv_out_yr, "ks_bogr/data//surv_pars.csv", row.names = F )


# * growth####

var_fe  <- data.frame( coefficient = names( 
           coef( gr_var ) ), value = coef( gr_var ) )

year_re <- data.frame( coefficient = paste0( "year", 
           rownames( coef( gr_mod_yr2 )$Year ) ), 
             value = coef( gr_mod_yr2 )$Year[,"(Intercept)"] )

la_re   <- data.frame( coefficient = paste0( "logsize_t0", 
           rownames( coef( gr_mod_yr2 )$Year ) ), 
             value = coef( gr_mod_yr2 )$Year[,"logsize_t0"] )

la2_re  <- data.frame( coefficient = paste0( "logsize_t0_2", 
           rownames( coef( gr_mod_yr2 )$Year ) ), 
             value = coef( gr_mod_yr2 )$Year[,"logsize_t0_2"] )


grow_out_yr <- Reduce ( function(...) 
               rbind  (...), 
               list ( var_fe, 
                      la_re, 
                      la2_re, 
                      year_re ) ) %>%
               mutate ( coefficient = as.character( coefficient ) )


write.csv( grow_out_yr, "ks_bogr/data/grow_pars.csv", row.names = F )


# * recruitment####

rc_pc       <- data.frame( coefficient = 
               paste0( "rec_pc_", repr_pc_yr$Year ),
                          value = repr_pc_yr$repr_percapita )

rc_sz       <- data.frame( coefficient = c( "rec_siz", "rec_sd" ),
                                 value = c( mean( recSize$logsize_t0 ),
                                              sd( recSize$logsize_t0 ) ) )

recr_out_yr <- Reduce ( function(...) 
               rbind(...), 
               list( rc_pc, 
                     rc_sz ) ) %>%
               mutate ( coefficient = as.character( coefficient ) )


write.csv( recr_out_yr, "ks_bogr/data/recr_pars.csv", row.names = F )


# * constant pars####

constants <- data.frame( coefficient = 
                     c ( "recr_sz",
                         "recr_sd",
                         "a",
                         "b",
                         "L",
                         "U",
                         "mat_size" ),
              value = c ( mean( recSize$logsize_t0 ),   
                      sd( recSize$logsize_t0 ), 
                      as.numeric(coef(gr_var)[1]),
                      as.numeric(coef(gr_var)[2]),
                      grow_df$logsize_t0 %>% min,  
                      grow_df$logsize_t0 %>% max, 200 ) )

surv_fe <- data.frame( coefficient = c( "surv_b0",
                                        "surv_b1" ),
                             value = fixef( su_mod_yr ) )

grow_fe <- data.frame( coefficient = c( "grow_b0",
                                        "grow_b1",
                                        "grow_b2" ),
                             value = fixef( gr_mod_yr2 ) )

rec_fe  <- data.frame( coefficient = "fecu_b0",
                             value = mean( repr_pc_yr$repr_percapita ) )

pars_cons <- Reduce( function(...) 
             rbind (...), 
             list  ( surv_fe, 
                     grow_fe, 
                     rec_fe, 
                     constants ) ) %>%
             mutate ( coefficient = as.character( coefficient))
             rownames( pars_cons ) <- 1:13

pars_cons_wide <- as.list( 
                  pivot_wider( pars_cons, 
                  names_from = "coefficient", 
                  values_from = "value" ) )


write.csv( pars_cons_wide, "ks_bogr/data/pars_cons.csv", row.names = F )


# * varying pars####

su_b0  <- data.frame( coefficient = 
          paste0( "surv_b0_", 
          rownames( coef( su_mod_yr ) $Year ) ), 
          value = coef  ( su_mod_yr ) $Year[,"(Intercept)"] )

su_b1   <- data.frame( coefficient = 
           paste0( "surv_b1_", 
           rownames( coef( su_mod_yr ) $Year ) ), 
           value = coef  ( su_mod_yr ) $Year[,"logsize_t0"] )

grow_b0 <- data.frame( coefficient = paste0( "grow_b0_", 
           rownames( coef( gr_mod_yr2 ) $Year ) ), 
           value = coef  ( gr_mod_yr2 ) $Year[,"(Intercept)"] )

grow_b1 <- data.frame( coefficient = 
           paste0( "grow_b1_", 
           rownames( coef( gr_mod_yr2 ) $Year ) ), 
           value = coef  ( gr_mod_yr2 ) $Year[,"logsize_t0"] )

grow_b2 <- data.frame( coefficient = 
           paste0( "grow_b2_", 
           rownames( coef( gr_mod_yr2 ) $Year ) ), 
           value = coef  ( gr_mod_yr2 ) $Year[,"logsize_t0_2"] )

fecu_b0   <- data.frame( coefficient = 
             paste0( "fecu_b0_", repr_pc_yr$Year ),
                         value = repr_pc_yr$repr_percapita )

pars_var  <- Reduce(function(...) 
             rbind(...), list( su_b0, 
                               su_b1, 
                               grow_b0, 
                               grow_b1, 
                               grow_b2, 
                               fecu_b0 ) )

pars_var_wide <- as.list( 
                 pivot_wider( pars_var, 
                 names_from = "coefficient", 
                 values_from = "value" ) )


write.csv( pars_var_wide, "ks_bogr/data/pars_var.csv", row.names = F )


# 12.  Building year specific IPM------------------------------------------------

lmer( logsize_t1 ~ logsize_t0 + logsize_t0_2 + ( logsize_t0 | Year ), data = grow_df )

grow_sd   <- function( x, pars ) {
             pars$a * ( exp( pars$b * x ) ) %>% sqrt 
}            # Standard deviation of growth model

gxy       <- function( x, y, pars ) {
             return( dnorm( y,  
             mean = pars$grow_b0 + pars$grow_b1*x + pars$grow_b2*x^2, 
             sd = grow_sd( x, pars ) ) )
}            # Growth from size x to size y

inv_logit <- function( x ) { exp( x ) / ( 1 + exp( x ) ) } # Inverse logit

sx        <- function( x, pars ) {
             return( inv_logit( pars$surv_b0 + pars$surv_b1 * x ) ) 
}            # Survival of x-sized individual to time t1


pxy <- function( x, y, pars ) {
       return( sx( x, pars ) * gxy( x, y, pars ) ) 
}      # Transition of x-sized individual to y-sized individual at time t1


fy       <- function( y, pars, h ){
  n_recr <- pars$fecu_b0
  recr_y <- dnorm( y, pars$recr_sz, pars$recr_sd ) * h 
  recr_y <- recr_y / sum( recr_y )
  f      <- n_recr * recr_y
            return( f )
  
}          # Per-capita production of y-sized recruits


kernel  <- function( pars ) {
  
  n  <- pars$mat_siz
  L  <- pars$L
  U  <- pars$U
  h  <- ( U - L ) / n
  b  <- L + c( 0:n ) * h
  y  <- 0.5 * ( b[1:n] + b[2:( n + 1 )] )
  
Fmat   <- matrix( 0, n, n )
Fmat[] <- matrix( fy( y, pars, h ), n, n )
  
Smat   <- c( )
Smat   <- sx( y, pars )
  
Gmat   <- matrix( 0, n, n )
Gmat[] <- t( outer( y, y, gxy, pars ) ) * h
  
Tmat   <- matrix( 0, n, n )
  
          for( i in 1:( n / 2 ) ) {
            
Gmat[1,i] <- Gmat[1,i] + 1 - sum( Gmat[,i] )
Tmat[,i]  <- Gmat[,i] * Smat[i] }
  
             for( i in ( n / 2 + 1 ):n ) {
               
Gmat[n,i] <- Gmat[n,i] + 1 - sum( Gmat[,i] )
Tmat[,i]  <- Gmat[,i] * Smat[i] }
  
     k_yx <- Fmat + Tmat
  
             return( list( k_yx = k_yx,
                Fmat = Fmat,
                Tmat = Tmat,
                Gmat = Gmat,
                meshpts = y ) ) }

pars_mean  <- pars_cons_wide

lambda_ipm <- function( i ) {
              return( Re( eigen( kernel( i )$k_yx )$value[1] ) )
}

lam_mean   <- lambda_ipm( pars_mean )
lam_mean
#1.113085

bou_gra_yr      <- c( 34:71 )
pars_yr         <- vector( mode = "list", length = length( bou_gra_yr ) )
extr_value_list <- function( x, field ) {
                   return( as.numeric( x[paste0( field )] %>% unlist( ) ) )
}

prep_pars   <- function( i ) {
  yr_now    <- bou_gra_yr[i]
  pars_year <- list( surv_b0 = extr_value_list( pars_var_wide, paste( "surv_b0", yr_now, sep = "_" ) ),
                     surv_b1 = extr_value_list( pars_var_wide, paste( "surv_b1", yr_now, sep = "_" ) ),
                     grow_b0 = extr_value_list( pars_var_wide, paste( "grow_b0", yr_now, sep = "_" ) ),
                     grow_b1 = extr_value_list( pars_var_wide, paste( "grow_b1", yr_now, sep = "_" ) ),
                     grow_b2 = extr_value_list( pars_var_wide, paste( "grow_b2", yr_now, sep = "_" ) ),
                     a       = extr_value_list( pars_cons_wide, "a" ),
                     b       = extr_value_list( pars_cons_wide, "b" ),
                     fecu_b0 = extr_value_list( pars_var_wide, paste( "fecu_b0", yr_now, sep = "_" ) ),
                     recr_sz = extr_value_list( pars_cons_wide, "recr_sz" ),
                     recr_sd = extr_value_list( pars_cons_wide, "recr_sd" ),
                     L       = extr_value_list( pars_cons_wide, "L" ),
                     U       = extr_value_list( pars_cons_wide, "U" ),
                     mat_siz = 200 )
             return( pars_year )
  
}

pars_yr     <- lapply( 1:length( bou_gra_yr ), prep_pars )

calc_lambda <- function( i ) {
        lam <- Re( eigen( kernel( pars_yr[[i]] )$k_yx )$value[1] )
               return( lam )
}

lambdas_yr  <- lapply( 1:length( bou_gra_yr ), calc_lambda )
               names( lambdas_yr ) <- bou_gra_yr

#comparing lambdas

year_kern  <- function( i ) {
              return( kernel( pars_yr[[i]] )$k_yx )
}

kern_yr    <- lapply( 1:length( bou_gra_yr ), year_kern )

all_mat    <- array( dim = c( 200, 200, 38 ) ) #38 years, not 13 rows

              for( i in 1:length( bou_gra_yr ) ) {
                
all_mat[,,i]  <- as.matrix( kern_yr[[i]] )
}

mean_kern     <- apply( all_mat, c( 1, 2 ), mean )

lam_mean_kern <- Re( eigen( mean_kern )$value[1] )
lam_mean_kern
#1.091588

# Population counts at time t0
pop_counts_t0 <- bou_gra %>%
                 group_by( Year, Quad ) %>%
                 summarize( n_t0 = n( ) ) %>% 
                 ungroup %>% 
                 mutate( Year = Year + 1 )

# Population counts at time t1
pop_counts_t1 <- bou_gra %>%
                 group_by( Year, Quad ) %>%
                 summarize( n_t1 = n( ) ) %>% 
                 ungroup 


# Calculate observed population growth rates, 
#  accounting for discontinued sampling!
pop_counts <- left_join( pop_counts_t0, 
                         pop_counts_t1 ) %>% # by dropping NAs, we remove gaps in sampling!
              drop_na %>% 
              group_by (  Year ) %>% 
              summarise( n_t0 = sum( n_t0 ),
                         n_t1 = sum( n_t1 ) ) %>% 
              ungroup %>% 
              mutate( obs_pgr = n_t1 / n_t0 ) %>% 
              mutate(  lambda = lambdas_yr %>% unlist )

lam_mean_yr      <- mean( pop_counts$lambda, na.rm = T )
lam_mean_count   <- mean( pop_counts$obs_pgr, na.rm = T )

lam_mean_geom    <- exp( mean( log( pop_counts$obs_pgr ), na.rm = T ) )
lam_mean_geom
#1.042462

lam_mean_overall <- sum( pop_counts$n_t1 ) / sum( pop_counts$n_t0 )
lam_mean_overall
# 0.99627

count_indivs_by_size <- function( size_vector, 
                                  lower_size, 
                                  upper_size, 
                                  matrix_size ){
  size_vector %>% 
    cut( breaks = seq( lower_size - 0.00001, 
                       upper_size + 0.00001, 
                       length.out = matrix_size + 1 ) ) %>% 
    table %>% 
    as.vector 
}

yr_pop_vec <- function( i ) {
  vec_temp <- surv_df %>% 
              filter( Year == i ) %>% 
              select( logsize_t0 ) %>% 
              unlist( ) 
  min_sz  <- pars_mean$L
  max_sz  <- pars_mean$U
  pop_vec <- count_indivs_by_size( vec_temp, min_sz, max_sz, 200 )
            return( pop_vec )
}

year_pop <- lapply( 34:71, yr_pop_vec )

proj_pop <- function( i ) {
             sum( all_mat[,,i] %*% year_pop[[i]] )
}

projected_pop_ns  <- sapply( 1:38, proj_pop ) #double check if 1:13 or 1:38 when pulled

pop_counts_update <- pop_counts %>% 
                     mutate ( proj_n_t1 = projected_pop_ns ) %>% 
                     mutate ( proj_pgr  = proj_n_t1/n_t0 ) 
                     ggplot ( pop_counts_update ) +
                     geom_point( aes( x = lambda,
                                      y = obs_pgr ),
                                  color = 'green' ) +
                    geom_point ( aes( x = proj_pgr,
                                      y = obs_pgr ),
                                  color = 'red' ) +
             geom_abline( aes(intercept = 0,
                                  slope = 1) ) +
                                labs( x = "Modeled lambda",
                                      y = "Observed population growth rate" )


# Building year specific IPM with 'ipmr'

all_pars <- c( pars_cons_wide, pars_var_wide )


write.csv( all_pars, "ks_bogr/data/all_pars.csv", row.names = F )


proto_ipm_yr <- init_ipm( sim_gen = "simple",
                            di_dd = "di",
                        det_stoch = "det" ) %>% 
               define_kernel(name = "P_yr",
                           family = "CC",
                          formula = s_yr * g_yr,
                             s_yr = plogis( surv_b0_yr + surv_b1_yr * size_1), 
                             g_yr = dnorm ( size_2, mu_g_yr, grow_sig ),
                          mu_g_yr = grow_b0_yr + grow_b1_yr * size_1 + grow_b2_yr * size_1^2,
                         grow_sig = sqrt( a * exp( b * size_1 ) ),
                        data_list = all_pars,
                           states = list( c( 'size' ) ),
                    uses_par_sets = TRUE,
                  par_set_indices = list( yr = 34:71 ),
                        evict_cor = TRUE,
                        evict_fun = truncated_distributions( 
                              fun = 'norm',
                           target = 'g_yr' ) ) %>% 
               define_kernel(name = 'F_yr',
                           family = 'CC',
                          formula = fecu_b0_yr * r_d,
                              r_d = dnorm( size_2, recr_sz, recr_sd ),
                        data_list = all_pars,
                           states = list( c( 'size' ) ),
                    uses_par_sets = TRUE,
                  par_set_indices = list( yr = 34:71 ),
                        evict_cor = TRUE,
                        evict_fun = truncated_distributions( "norm", "r_d" ) ) %>% 
               define_impl(
               make_impl_args_list(
                     kernel_names = c( "P_yr", "F_yr" ),
                         int_rule = rep( "midpoint", 2 ),
                     state_start  = rep( "size", 2 ),
                        state_end = rep( "size", 2 ) ) ) %>% 
              define_domains(size = c(all_pars$L,
                                    all_pars$U,
                                    all_pars$mat_siz )) %>% 
              define_pop_state(
                        n_size_yr = rep( 1 / 200, 200 ) )


ipmr_yr       <- make_ipm( proto_ipm = proto_ipm_yr,
                          iterations = 200 )
lam_mean_ipmr <- lambda( ipmr_yr )

lam_out       <- data.frame( coefficient = 
                 names ( lam_mean_ipmr ), 
                 value = lam_mean_ipmr )
                rownames( lam_out ) <- 1:38

lam_out_wide <- as.list( 
                pivot_wider( lam_out, 
                names_from = "coefficient", 
                values_from = "value" ) )


write.csv( lam_out_wide, "ks_bogr/data/lambdas_yr.csv", row.names = F )


#Populating the PADRINO database template
# SKIP for now as I don't have all the relevant information.


