get_migration_table <- function() {
  data.frame(
    id = 445,
    year = c(2014, 2014, 2015, 2015),
    season = c("spring", "fall", "spring", "fall"),
    date_start = c("2014-05-04", "2014-11-02", "2015-04-22", "2015-11-13"),
    date_end = c("2014-05-17", "2014-11-21", "2015-05-13", "2015-11-24")
  ) %>%
    mutate(
      date_start = ymd(date_start, tz = "America/Denver"),
      date_end = ymd(date_end, tz = "America/Denver")
    ) %>%
    group_by(id) %>%
    arrange(year, desc(season)) %>%
    mutate(year_cumm = year - min(year) + 1)
}