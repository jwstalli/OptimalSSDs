#########################################################
# In the following example, I am generating 2 criteria pareto front
# Optimal designs for p factors each with 3 levels, we are assuming linearity here.
# (-1,0, 1). The 2 criteria are the sum of the interpolated Psi for each off diagonal of an exact C matrix,
# and the Trace of the V matrix
# Let me know if you still have any questions 
##########################################################



#### Change paths below
source('/Users/hkyoung/Desktop/SIC_Paper_Code/SIC_Paper_Function_Library.R')
library(abind)

# Creating a folder for the output
output_folder = "./n9_k10_kstar3_allpos_designs_EPCEA_HILS_3_levels_new_measure"

if(file.exists(output_folder)){
  print("Output folder already exists, overwritting output")
}else{dir.create(output_folder)}


# Set the dimensions of the designs
n = 9L
p = 10L
k=3L
B=3
level = c(-1L,0L, 1L)


# Set the log lamdda list to generate curves over

log_lambda_list = seq(-3,2, by = 0.05)
lambda_list = exp(log_lambda_list)


#' This function takes in user inputs and does the EPCEA for 2 or 3 level designs (all columns have same number of levels).
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
#' @param frob_or_sum: string with a value of "frob" or "sum" to delineate what pareto criteria to use
#' "sum" is the new one, and is default
#'
#' @return list with items "df_curve" a dataframe that gives the psi or psi_pm 
#' (depending on sign vectors chosen) for each c in the c grid and also a "c_curve_plot" 
#' ggplot object that shows the psi_corve as a function of c. Also has "psi_max",
#' the maximum value of psi over the grid, and "c_opt" the optimal value of c. Additionally, 
#' has "pareto_front_designs which is a list of the designs on the pareto front, and 
#' "pareto_front_results" which show the criteria values for each of the PF designs.
#' "pareto_front_designs" are the unlisted pareto front designs for secondary measures.
#' Lastly, "EPCEA_results" is a list that is used to further evaluate the output.

EPCEA_results_new_measure <- EPCEA_wrapper(n, p, k, lambda_list, B1=B,level,ns=30, signs = "all_pos", 
                                                       num_c_for_grid = 150, sigma_sq = 1, cores = 8,
                                                       frob_or_sum = "sum")

# Here are the designs on the pareto front
candidate_designs_list = EPCEA_results_new_measure$pareto_front_designs_list
candidate_designs = EPCEA_results_new_measure$pareto_front_designs
saveRDS(EPCEA_results_new_measure, file = "./n9_k10_kstar3_allpos_designs_EPCEA_HILS_3_levels_new_measure/EPCEA_3L_results_n9_p10_new_measure.RData")

#Next line is to read in results if you have already run it
#EPCEA_results_new_measure = readRDS("./n9_k10_kstar3_allpos_designs_EPCEA_HILS_3_levels_new_measure/EPCEA_3L_results_n9_p10_new_measure.RData")



# Read in other designs to compare with the PF designs

d1<-as.matrix(read_in_design("./n9_k10_designs/d1.txt"))
VarSPlus<-as.matrix(read_in_design("./n9_k10_designs/n9_k10_vars_pos_Eff80_100.txt"))
UES_sq_opt<-as.matrix(read_csv("./n9_k10_designs/UES_sq_opt.csv"))

# Add the named designs in a list to the PF designs
candidate_designs<- abind(candidate_designs, d1, VarSPlus, UES_sq_opt, along=3)
candidate_designs_list[[length(candidate_designs_list)+1]]=d1
candidate_designs_list[[length(candidate_designs_list)+1]]=VarSPlus
candidate_designs_list[[length(candidate_designs_list)+1]]=UES_sq_opt

# name them for tracking
des.name.list = generate_names("pf", length=length(candidate_designs_list))
des.name.list[length(candidate_designs_list)-2]<-"PED"
des.name.list[length(candidate_designs_list)-1]<-"Var(s+)"
des.name.list[length(candidate_designs_list)]<-"UE(s^2)"

df_results <- eval_candidate_designs(candidate_designs_list,candidate_designs,
                                                 des.name.list, 
                                                 n, p, k, 
                                                 lambda_list, B1=B, 
                                     df_curve= EPCEA_results_new_measure$df_curve,
                                     signs = "all_pos",
                                                 sigma_sq = 1, cores = 8)

# Plot of the phi_lambda for all designs
ggjoint = ggplot(data = df_results, aes(x=log(as.numeric(lambda)), y= as.numeric(P_joint), color=design, linetype = design))+
  geom_line(size=1.2)+
  xlab(latex2exp::TeX("$log(\\lambda)"))+
  ylab(latex2exp::TeX("$\\Phi_{\\lambda}"))+
  #scale_color_manual(values =c("blue","red"))+
  theme_bw()+
  scale_x_continuous(limits = c(-2,2))+
  #scale_linetype_manual(values =c("solid", "dashed"))+
  theme( legend.title=element_blank(), text=element_text(size=20))


# get the integrated area of phi_curve
avg_P_joint <- df_results |>
  mutate(P_joint = as.numeric(P_joint)) |>
  group_by(design, k_star, B) |>
  summarise(mean_P_joint = mean(P_joint, na.rm = TRUE), .groups = "drop")

