#########################################################
# In the following code I am giving demonstrations about the Psi value as a function 
# of a compound symmetric off-diagonal,c. Then comparing existing design's correlations
# to the optimal c curve. Then I give the new criterion for an exact C matrix of
# sum(psi(c_ij)) for existing designs. 
# Let me know if you still have any questions 
##########################################################


#### Change paths below
source('/Users/hkyoung/Desktop/SIC_Paper_Code/SIC_Paper_Function_Library.R')
library(abind)


# Set the dimensions of the designs
n = 9L
p = 10L
k=3L
B=3


# Set the log lamdda list to generate curves over

log_lambda_list = seq(-3,2, by = 0.05)
lambda_list = exp(log_lambda_list)

#' The function below "get_c_curve_plot" is a function that calculates the psi
#' (or psi_pm for the all signs case), integrated over a lambda grid for each c 
#' in a given grid. Below are the inputs and the outputs. Currently this only 
#' works for integrating over lambda, but will update to take the max over lambda.
#'
#' @param n: positive integer; number of experimental runs (rows of X)
#' @param p: positive integer; number of factors (cols of X)
#' @param k: positive integer; number of assumed active effects
#' @param lamba_list: vector; grid of lambda to optimize over either by integration or point max
#' @param B1: positive scalar; magnitude of active effects
#' @param signs: string with a value of "all_pos" or "all". If "all_pos" only the 
#' positive sign is used. If "all", the measure is averaged over all posible sign vecs.
#' @param num_c_for_grid: positive integer; number of points on the c-grid. Default is 200
#' @param sigma_sq: positive scalar; variance. Default is 1
#' @param cores: positive integer; number of cores to calculate the c_curve in parallel
#'
#' @return list with items "df_curve" a dataframe that gives the psi or psi_pm 
#' (depending on sign vectors chosen) for each c in the c grid and also a "gg1" 
#' ggplot object that shows the psi_corve as a function of c. Also has "psi_max",
#' the maximum value of psi over the grid, and "c_opt" the optimal value of c. 
#' @export
#'

output <- get_c_curve_plot(n, p, k, lambda_list, B1=B, signs="all_pos", 
                           num_c_for_grid=200, sigma_sq=1, cores = 4)





# Demonstration on specific designs from paper these are the designs from the original paper
# HILS design is with the old criterion. 
# Update path to these designs

d1<-as.matrix(read_in_design("./n9_k10_designs/d1.txt"))
VarSPlus<-as.matrix(read_in_design("./n9_k10_designs/n9_k10_vars_pos_Eff80_100.txt"))
HILS_n9 <- as.matrix(read_csv("./n9_k10_kstar3_allpos_designs_EPCEA_HILS/n9_k10_allpos_HILS.csv"))

design_list <- list(
  "DCD"       = d1,
  "VarSPlus" = VarSPlus,
  "HILS"  = HILS_n9
)


criteria_comparison <- evaluate_designs_new_approach( design_list, 
                                                      df_curve= output$df_curve, 
                                                      c_opt=output$c_opt)
print(criteria_comparison)

# creating the plot that compares the actual correlations to the c regions
# need to wrap this, but this will do for now. 
df_hist = corr_df <- do.call(rbind, lapply(names(design_list), function(design_name) {
  design       <- design_list[[design_name]]
  correlations <- get_correl(design)
  
  data.frame(
    design      = design_name,
    cor = correlations
  )
}))

df_density_scaled <- df_hist |>
  group_by(design) |>
  reframe(
    density_obj = list(density(as.numeric(cor), adjust = 0.5, from = -1, to = 1)),
    x           = density_obj[[1]]$x,
    y           = density_obj[[1]]$y
  )

# Compute scaling factor to match maxes
psi_max     <- output$psi_max
density_max <- max(df_density_scaled$y)
scale_factor <- psi_max / density_max

# Overlay with sec_axis
gg_overlay <- output$c_curve_plot +
  geom_line(data = df_density_scaled,
            aes(x = x, y = y * scale_factor, group = design, color = design),
            lwd = 1, inherit.aes = FALSE) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . / scale_factor, name = "Density")
  ) +
  theme(legend.title = element_blank())

gg_overlay




#Second Demonstration with the same dimensions as above, but with all signs


output2 <- get_c_curve_plot(n, p, k, lambda_list, B1=B, signs="all", 
                           num_c_for_grid=200, sigma_sq=1, cores = 4)

output2$c_curve_plot

# Demo 3
# n=14 p=20, k=5, B=2 all pos signs demo
n = 14L
p = 20L
k=5L
B=2

log_lambda_list = seq(-3,2, by = 0.05)
lambda_list = exp(log_lambda_list)

# At this size this takes a longer time. I used 8 cores here. 
output3 <- get_c_curve_plot(n, p, k, lambda_list, B1=B, signs="all_pos", 
                           num_c_for_grid=150, sigma_sq=1, cores = 8)


output3$c_curve_plot



# Demonstration on specific designs from paper these are the designs from the original paper
# HILS design is with the old criterion. 
# Update path to these designs


UES_sq_n14 <- as.matrix(read_csv("./n14_p20_kstar5_allpos_designs_EPCEA/UES_sq_opt_n14_p20.csv"))
HILS_n14 <- as.matrix(read_csv("./n14_p20_kstar5_allpos_designs_EPCEA/n14_p20_allpos_HILS.csv"))
d1_n14<-as.matrix(read_in_design("./n14_k20_kstar5_allpos_designs_HILS/d1.txt"))
VarSPlus_n14 <- as.matrix(read_csv("./n14_k20_kstar5_allpos_designs_HILS/best_huer_Var_s_design_18.csv"))

design_list3 <- list(
  "DCD"       = d1_n14,
  "VarSPlus" = VarSPlus_n14,
  "HILS"  = HILS_n14,
  "UES_sq" = UES_sq_n14
)


criteria_comparison3 <- evaluate_designs_new_approach( design_list3, 
                                                      df_curve= output3$df_curve, 
                                                      c_opt=output3$c_opt)
print(criteria_comparison3)

# creating the plot that compares the actual correlations to the c regions
# need to wrap this, but this will do for now. 
df_hist3 = do.call(rbind, lapply(names(design_list3), function(design_name) {
  design       <- design_list3[[design_name]]
  correlations <- get_correl(design)
  
  data.frame(
    design      = design_name,
    cor = correlations
  )
}))

df_density_scaled3 <- df_hist3 |>
  group_by(design) |>
  reframe(
    density_obj = list(density(as.numeric(cor), adjust = 0.5, from = -1, to = 1)),
    x           = density_obj[[1]]$x,
    y           = density_obj[[1]]$y
  )

# Compute scaling factor to match maxes
psi_max     <- output3$psi_max
density_max <- max(df_density_scaled$y)
scale_factor <- psi_max / density_max

# Overlay with sec_axis
gg_overlay3 <- output3$c_curve_plot +
  geom_line(data = df_density_scaled3,
            aes(x = x, y = y * scale_factor, group = design, color = design),
            lwd = 1, inherit.aes = FALSE) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . / scale_factor, name = "Density")
  ) +
  theme(legend.title = element_blank())

gg_overlay3