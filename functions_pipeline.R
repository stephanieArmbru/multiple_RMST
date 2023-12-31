########################################################################
# PIPELINE FOR RMST CALCULATION FOR MULTIPLE EVENT TYPES

## FUNCTIONS ##

# Prof. LJ Wei, Dr. Zack McCaw, Stephanie Armbruster 
########################################################################



# Formatting --------------------------------------------------------------
construct_composite_endpoints <- function(data, 
                                          t2event_names, 
                                          c4event_names, 
                                          t2death_name = "t2dth", 
                                          c4death_name = "c4dth",
                                          composite_name) {
  K <- length(t2event_names)
  
  # guarantee that all times are above 0
  # loop through all event types 
  for (k in seq(1, K)) {
    t2event_name <- t2event_names[k]
    c4event_name <- c4event_names[k]
    
    data <- data %>% 
      filter(get(t2event_name) > 0 &
               get(t2death_name) > 0)
  } # end loop through all event types 
  

  # loop through all event types 
  for (k in seq(1, K)) {
    t2event_name <- t2event_names[k]
    c4event_name <- c4event_names[k]
    
    # generate composite time to either event or death 
    data[, paste0(t2event_name, composite_name)] <- data %>% 
      dplyr::select(c(all_of(t2event_name), 
                      all_of(t2death_name))) %>% 
      apply(1, FUN = min)
    
    # generate composite event indicator for either event or death 
    data[, paste0(c4event_name, composite_name)] <- data %>% 
      dplyr::select(c(all_of(c4event_name), 
                      all_of(c4death_name))) %>% 
      apply(1, FUN = max)
  } # end loop through all event types 
  return(data)
}

# rename column indicating trial arm 
rename_arm <- function(data, 
                       arm) {
  data$arm <- data %>% 
    pull(all_of(arm))
  return(data)
}


# Truncation time ---------------------------------------------------------
find_truncation_time <- function(data, 
                                 t2event_names,
                                 c4event_names) {
  
  event_tau <- numeric(length = length(t2event_names))
  
  K <- length(t2event_names)
  # loop through all event types 
  for (k in seq(1, K)) {
    
    t2event_name <- t2event_names[k]
    c4event_name <- c4event_names[k]
    
    # max time-to-composite-event or censoring observed for treatment and 
    # placebo group 
    max_observations <- data %>% 
      group_by(arm) %>% 
      filter(get(t2event_name) == max(get(t2event_name)))
    
    t2event_arm1 <- max_observations %>% 
      filter(arm == 1) %>% 
      pull(get(t2event_name)) %>% 
      unique()
    c4event_arm1 <- max_observations %>% 
      filter(arm == 1) %>% 
      pull(get(c4event_name)) %>% 
      unique()
    
    t2event_arm0 <- max_observations %>% 
      filter(arm == 0) %>% 
      pull(get(t2event_name)) %>% 
      unique()
    c4event_arm0 <- max_observations %>% 
      filter(arm == 0) %>% 
      pull(get(c4event_name)) %>% 
      unique()
    
    # 1. case: event happens for both groups
    if (c4event_arm0 == 1 & 
        c4event_arm1 == 1) {
      event_tau[k] <- max(t2event_arm0, t2event_arm1)
      
      # 2. case: event happens in one group, censoring in another, and 
      # censoring happens AFTER event
    } else if (c4event_arm0 == 0 & 
               c4event_arm1 == 1 & 
               t2event_arm1 <= t2event_arm0) {
      event_tau[k] <- max(t2event_arm0, t2event_arm1)
    } else if (c4event_arm0 == 1 & 
               c4event_arm1 == 0 & 
               t2event_arm1 >= t2event_arm0) {
      event_tau[k] <- max(t2event_arm0, t2event_arm1)
      
      # 3. case: event happens in one group, censoring in another, and 
      # censoring happens BEFORE event 
    } else if (c4event_arm0 == 0 & 
               c4event_arm1 == 1 & 
               t2event_arm1 >= t2event_arm0) {
      event_tau[k] <- min(t2event_arm0, t2event_arm1)
    } else if (c4event_arm0 == 1 & 
               c4event_arm1 == 0 & 
               t2event_arm1 <= t2event_arm0) {
      event_tau[k] <- min(t2event_arm0, t2event_arm1)
    } else if (c4event_arm0 == 0 &
               c4event_arm1 == 0) {
      event_tau[k] <- min(t2event_arm0, t2event_arm1)
    }
  } # end loop through all event types 
  
  return(min(event_tau))
}


# Kaplan-Meier survival time ----------------------------------------------
KM_survival <- function(data, 
                        t2event_names, 
                        c4event_names, 
                        
                        alpha = 0.05, 
                        nar_grid = 10, 
                        round_to = 50) {
  # separate data set
  data_arm0 <- data %>% filter(arm == 0)
  data_arm1 <- data %>% filter(arm == 1)
  
  K <- length(t2event_names)
  KM_S <- list()
  curves <- list()
  
  # estimate Kaplan-Meier survival 
  # loop through all event types
  for (k in seq(1, K)) {
    # pull time, censoring time and status for each arm 
    time_arm0_k <- data_arm0 %>% 
      pull(get(t2event_names[k]))
    status_arm0_k <- data_arm0 %>% 
      pull(get(c4event_names[k]))
    
    
    time_arm1_k <- data_arm1 %>% 
      pull(get(t2event_names[k]))
    status_arm1_k <- data_arm1 %>% 
      pull(get(c4event_names[k]))
    
    
    
    # estimate survival probability for each arm 
    arm0_k <- survfit(Surv(time_arm0_k, status_arm0_k) ~ 1, 
                      data = data_arm0, 
                      conf.type = "plain", 
                      conf.int = 1 - alpha)
    arm0_k_df <- data.frame(time = arm0_k$time, 
                            surv = arm0_k$surv, 
                            lower = arm0_k$lower,
                            upper = arm0_k$upper)
    
    arm1_k <- survfit(Surv(time_arm1_k, status_arm1_k) ~ 1, 
                      data = data_arm1, 
                      conf.type = "plain", 
                      conf.int = 1 - alpha)
    arm1_k_df <- data.frame(time = arm1_k$time, 
                            surv = arm1_k$surv, 
                            lower = arm1_k$lower,
                            upper = arm1_k$upper)
    
    KM_S[[t2event_names[k]]] <- list(arm0 = arm0_k_df, 
                                     arm1 = arm1_k_df)
    
    
    # visualize KM survival and save 
    g <- ggplot() +
      geom_line(data = arm1_k_df, 
                aes(x = time, 
                    y = surv, 
                    color = "Treatment (arm 1)")) +
      geom_ribbon(data = arm1_k_df,
                  aes(x = time, 
                      y = surv,
                      ymin = lower,
                      ymax = upper, 
                      fill = "Treatment (arm 1)"),
                  alpha = 0.3) +
      geom_line(data = arm0_k_df, 
                aes(x  = time, 
                    y = surv, 
                    color = "Placebo (arm 0)")) +
      geom_ribbon(data = arm0_k_df,
                  aes(x = time, 
                      y = surv,
                      ymin = lower,
                      ymax = upper, 
                      fill = "Placebo (arm 0)"),
                  alpha = 0.3) +
      theme_bw() +
      theme(legend.position = "bottom") +
      labs(color = "", 
           fill = "",
           y = "KM survival (pointwise CI)") +
      ylim(0, 1)
    
    # add table with numbers at risk 
    # generate reporting times 
    reporting_times <- seq(min(c(arm1_k$time, 
                                 arm0_k$time)), 
                           max(c(arm1_k$time, 
                                 arm0_k$time)), 
                           length = nar_grid) %>% 
      plyr::round_any(round_to)
    # assign numbers at risk to reporting times
    fx_arm1 <- arm1_k$n.risk[findInterval(reporting_times, 
                                               c(-Inf, arm1_k$time))]
    fx_arm0 <- arm0_k$n.risk[findInterval(reporting_times, 
                                               c(-Inf, arm0_k$time))]
    
    fx_df <- data.frame(time = rep(reporting_times, 3), 
                        nar = c(reporting_times, 
                                fx_arm1, 
                                fx_arm0), 
                        arm =  c(rep(1, nar_grid), 
                                 rep(0.5, nar_grid),
                                 rep(0, nar_grid))
                        )
    
    tbl <- ggplot(fx_df,  
                  aes(x = time, 
                      y = factor(arm),
                      label = nar)) +
      geom_text() +
      theme_bw() + 
      theme( 
        legend.position = "none",
        plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, 
                                   face="bold", 
                                   color = 'black'),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_blank()
      ) + 
      scale_y_discrete(breaks = c(0, 0.5, 1), 
                       labels = c("Arm 0", 
                                  "Arm 1", 
                                  "Time")) +
      scale_x_continuous(breaks = reporting_times)
    
    # combine KM survival curve and table for number at risk
    g_combined <- cowplot::plot_grid(
      plotlist = list(g, tbl),
      align = "v",
      axis = "l",
      ncol = 1,
      rel_heights = c(3, 1)
    )
    curves[[t2event_names[k]]] <- g_combined
    
    # save plot 
    ggsave(paste0("Figures/KM_survival_", 
                  t2event_names[k], 
                  ".pdf"), 
           plot = g_combined, 
           height = 13, width = 17, unit = "cm")
    } # end loop through all event types
  
  return(list(estimate = KM_S, 
              curves = curves))
}



# Hazard Ratio ------------------------------------------------------------
HR <- function(data, 
               t2event_names, 
               c4event_names) {
  
  K <- length(t2event_names)
  HR_est <- list()
  
  # estimate Cox proportional hazard model
  # loop through all event types
  for (k in seq(1, K)) {
    # pull time, censoring time and status
    time <- data %>% 
      pull(get(t2event_names[k]))
    status <- data %>% 
      pull(get(c4event_names[k]))
    arm <- data$arm
    
    CPHM <- coxph(Surv(time, status) ~ arm) %>% 
      print_cox_outputs()
    HR_est[[k]] <- data.frame(event = t2event_names[k],
                              HR = CPHM$HR, 
                              p = CPHM$pvalue, 
                              CI = paste0(round(CPHM$CI, 3), collapse = ","))
  } # end loop through all event types
  
  return(do.call(rbind, HR_est))
}


# WLW procedure -----------------------------------------------------------
WLW_TwoSamples <- function(data, 
                           time, 
                           status) {
  
  data_wlw <- data %>% 
    select(arm, 
           all_of(time), 
           all_of(status)) %>% 
    mutate(idx = seq_len(nrow(data)))
  
  # prepare time
  data_wlw_time <- data_wlw %>%
    dplyr::select(idx, 
                  arm, 
                  all_of(time)) %>% 
    tidyr::pivot_longer(
      cols = all_of(time), 
      names_to = "outcome",
      values_to = "time"
    )
  
  data_wlw_time$outcome <- gsub("t2", "", 
                                data_wlw_time$outcome)
  
  # prepare status
  data_wlw_status <- data_wlw %>%
    dplyr::select(idx, 
                  all_of(status)) %>% 
    tidyr::pivot_longer(
      cols = all_of(status), 
  names_to = "outcome",
      values_to = "status"
    )
  data_wlw_status$outcome <- gsub("c4", "", data_wlw_status$outcome)
  
  # merge to final data
  data_wlw <- data_wlw_time %>% 
    dplyr::inner_join(data_wlw_status, 
                      by = c("idx", 
                             "outcome"))
  
  # fit Cox proportional hazard model 
  fit_wlw <- coxph(Surv(time, status) ~ arm + cluster(idx), 
                   data = data_wlw) 
  
  return(fit_wlw)
}




