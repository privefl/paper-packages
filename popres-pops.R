list.pops <- list(
  "SW Europe" = c("Portugal", "Spain"),
  "Italy" = "Italy",
  "SE Europe" = c("Greece", "Turkey", "Cyprus", "Macedonia"),
  "Anglo-Irish Isles" = c("United", "Scotland", "Ireland"),
  "FennoScandia" = c("?orway", "Finland", "Sweden", "Denmark"),
  "Western Europe" = c("France", "Swiss-French", "Swiss-German",
                       "Swiss-Italian", "Belgium"),
  "Central Europe" = c("Germany", "Austria", "Poland", "?etherlands"),
  "Eastern Europe" = c("Czech", "Romania", "Bulgaria", "Hungary",
                       "Albania",  "Macedonia,", "Slovakia"),
  "Former Yugoslavia" = c("Croatia", "Serbia", "Bosnia", "Slovenia", "Kosovo"),
  "Former USSR" = c("Russian", "Ukraine", "Latvia")
)
mat.pops <- cbind(
  unlist(list.pops, use.names = FALSE),
  rep(names(list.pops), times = sapply(list.pops, length))
)
ind.match <- match(popres$fam$family.ID, mat.pops[, 1])
Population <- mat.pops[ind.match, 2]