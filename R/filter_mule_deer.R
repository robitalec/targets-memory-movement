filter_mule_deer <- function(mule_deer, migration_table, year, season) {
  mule_deer

  sub_migration_table <- migration_table[
    migration_table$year_cumm == year &
      migration_table$season == season,]

  sub_mule_deer <- mule_deer[
    mule_deer$date >= sub_migration_table$date_start &
      mule_deer$date <= sub_migration_table$date_end,
  ]

  sub_mule_deer$season <- season

  return(sub_mule_deer)

}
