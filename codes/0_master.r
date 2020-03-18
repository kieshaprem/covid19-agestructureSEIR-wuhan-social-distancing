# load relevant the data files
source('codes/1_loadData.r')

# source the age-structured SEIcIscR model functions 
source('codes/function_modelSEIcIscR.r')

# source the age-structured SEIcIscR model functions 
source('codes/function_postprocessing.r')

# simulate N oubtreaks
source('codes/2_simOutbreak_ncov_SEIR.r')
source('codes/2_simOutbreak_ncov_SEIcIscR.r')





