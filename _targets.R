#targets pipeline setup

library(targets)
library(crew)

tar_source() #configure metadata with matching parameters
tar_option_set(packages = c("readr", "dplyr", "MatchIt", "TBSS", "optmatch"),
               controller = crew_controller_local(workers = 2) #specify workers
)

list(
  tar_target(file, 
             "data/test_data_j.csv",    #use your own data
             format = "file"
  ),
  tar_target(data, 
             get_data(file)
  ),
  tar_target(iteration_id, 
             1:10        #specify number of replicates
  ),
  tar_target(subsample, 
             select_subsample(data, 1000, 1000), #specify sample sizes
             pattern = map(iteration_id)
  ),
  tar_target(sim_data, 
             plasmode(subsample,
                      target_outcome = NULL, #specify either outcome or incidence
                      target_inc = .01, 
                      rr=2),                 #specify relative risk
             pattern = map(subsample)
  ),
  tar_target(ps_data,
             ps_match(sim_data,
                      match_args),
             pattern = map(sim_data)
  ),
  tar_target(tree_file, 
             'data/test_tree.csv',                   #include tree file
             format = "file"
  ),
  tar_target(tree, 
             get_data(tree_file)
  ),
  tar_target(tbss, 
             tbss_mask(ps_data, 
                       model = "bernoulli",
                       tree), 
             pattern = map(ps_data)
  ),
  tar_target(final_vector, 
             unlist(tbss)
  ),
  tar_target(final_power, 
             print(mean(final_vector))
  )
  
)

