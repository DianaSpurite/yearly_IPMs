# 1. Setting up workspace and data----------------------------------------------------------------------------------------

#Project ks_ange
#By Diana Spurite
#diana.spurite@posteo.de

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


plot <- and_ger_long %>% ggplot( aes( x = size_value ) ) +
  geom_histogram( binwidth = 1 ) +
  facet_grid( Year_fac ~ size, 
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
surv_bin_yrs   <- lapply( 1:13, df_binned_prop_year, and_ger, 20, 
                          logsize_t0, survives_tplus1, surv_yrs )

surv_yr_pan_df <- bind_rows( surv_bin_yrs ) %>% 
  mutate( transition = paste( paste0( Year ),
                              substr( paste0( Year + 1 ), 3, 4 ),
                              sep = '-' ) ) %>% 
  mutate( year       = as.integer( Year - 31  ) )


plot <- ggplot( data   = surv_yr_pan_df, 
        aes( x = logsize_t0, 
             y = survives_tplus1 ) ) +
  geom_point( alpha = 0.5,
              pch   = 16,
              size  = 1,
              color = 'red' ) +
  scale_y_continuous( breaks = c( 0.1, 0.5, 0.9 ) ) +
  # split in panels
  facet_wrap( .~ transition, nrow = 4 ) +
  theme_bw( ) +
  theme( axis.text = element_text( size = 8 ),
         title     = element_text( size = 10 ),
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


ggsave(filename = "ks_ange/results/survival_binned_yr.png", plot = plot)

#code works, but the graphs are not displayed, but are saved in the "results" folder. Dimenstions bad.


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
         facet_wrap( .~ transition, nrow = 4 ) +
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
        facet_wrap( .~ Year, nrow = 4 ) +
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
              nrow = 4 ) +
        labs( x = expression('log( size )'[t0]),
              y = "Frequency" )

ggsave(filename = "ks_ange/results/recruit_histograms.png", plot = plot)

# the plots so far have not-so-very-nice dimensions, 
# but recycling "width = 10, height = 6, units = "in", res = 150" doesn't do any good for ggsave

# 6. Fitting vital rate models
##### Survival model####-----------------------

surv_df      <- surv %>% 
  mutate( logsize_t0 = log( basalArea_genet ) )

surv_mod_yr <- glmer( survives_tplus1 ~ logsize_t0 + ( logsize_t0 | Year ), data = surv_df, family = binomial )

ranef_su <- data.frame( coef( surv_mod_yr )[1] )
years_v  <- c( 32:69 )

surv_yr_plots <- function( i ){
  surv_temp   <- as.data.frame( surv_bin_yrs[[i]] )
  x_temp      <- seq( min( surv_temp$logsize_t0, na.rm = T ), 
                      max( surv_temp$logsize_t0, na.rm = T ), 
                      length.out = 100)
  pred_temp   <- boot::inv.logit( ranef_su[i,1] + ranef_su[i,2] * x_temp ) 
  pred_temp_df <- data.frame( logsize_t0 = x_temp, survives_tplus1 = pred_temp )
  temp_plot <- surv_temp %>% ggplot( ) +
    geom_point( aes( x = logsize_t0, y = survives_tplus1 ) ) +
    geom_line( data = pred_temp_df, aes( x     = logsize_t0,
                                         y     = survives_tplus1 ),
               color = 'red',
               lwd   = 1  ) +
    labs( title = paste0( years_v[i] ),
          x = expression( 'log( size )'[t0] ),
          y = expression( 'Survival probability  '[ t1] ) )
  if( i %in% c( 2:4, 6:8, 10:12 ) ){
    temp_plot <- temp_plot + theme( axis.title.y = element_blank( ) )
  }
  
  return(temp_plot)
}
length(surv_bin_yrs)
surv_yrs <- lapply( 1:13, surv_yr_plots )
surv_years <- wrap_plots( surv_yrs ) + plot_layout( nrow = 4 )

#only shows 13 plots, need 38, but replacing 13 to 38 gives "out of bounds" error

surv_years


ggsave(filename = "ks_ange/results/survival_pred.png", plot = surv_years)
#code works, but the plot's gone missing, I will come back to this 

# I choose not to open the pandora's box of the 13 warnings


#### Growth model####------------------------------------------------

grow_df      <- grow %>% 
  mutate( logsize_t0 = log( basalArea_genet ),
          logsize_t1 = log( size_tplus1 ) )

grow_mod_yr <- lmer( logsize_t1 ~ logsize_t0 + ( logsize_t0 | Year ), data = grow_df )


ranef_gr <- data.frame( coef( grow_mod_yr )[1] )

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
  if( i %in% c( 1998:2000, 2002:2004, 2006:2008 ) ){
    temp_plot <- temp_plot + theme( axis.title.y = element_blank( ) )
  }
  
  return(temp_plot)
}
grow_yrs <- lapply( 1997:2009, grow_yr_plots )
grow_years <- wrap_plots( grow_yrs ) + plot_layout( nrow = 4 )

# png( 'results/Bou_gra_yr/grow_pred.png', width = 10, height = 6, units = "in", res = 150 )

grow_years

ggsave(filename = "ks_ange/results/survival_pred.png", plot = grow_years)
#code works, but the plot's gone missing, I will come back to this.
# these plots were there, but after a re-run can't be made as the data object is a list.





