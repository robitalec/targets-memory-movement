# === Targets: Memory Movement --------------------------------------------
# Alec L. Robitaille
# Based on https://github.com/kimx3725/Memory_Movement
# Kim, D., Thompson, P., Wolfson, D., Merkle, J., Oliveira-Santos, L.G.R.,
#  Forester, J., Avgar, T., Lewis, M. and Fieberg, J., 2023.
#  Identifying signals of memory from observations of animal movements
#  in Plato's cave. bioRxiv, pp.2023-08.



# Source ------------------------------------------------------------------
targets::tar_source('R')




# Variables ---------------------------------------------------------------
path_mule_deer_445 <- 'data/Mule_Deer-Merkle_etal_2019/GPSdata_muldeer445.rds'




# Adapting Mule_Deer-Merkle_etal_2019/MemorySSF_workflow.R ----------------
c(
  tar_target(
    mule_deer_445,
    readRDS(path_mule_deer_445)
  ),
  tar_target(
    migration_table,
    get_migration_table()
  ),
  tar_target(
    filter_spring_second_year,
    filter_mule_deer(mule_deer_445, migration_table, 2, 'spring')
  )
)

