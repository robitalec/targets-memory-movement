cast_line <- function(pts) {
  pts %>%
    st_combine() %>%
    st_cast("LINESTRING") %>%
    data.frame(id=mig.i$id,
               season="spring",
               year=mig.i$year,
               geometry=.) %>%
    st_as_sf(sf_column_name = "geometry")
}