path <- 'data/spectra/'
spectra <- read.csv(path, skip = 11, header = T)

df <- spectra

# Changing all the '-9999.0' or '-9999' (missing data) to NA
for (i in 1:length(df)){
  df[i][df[i] == '-9999' | df[i] == '-9999.0'] <- NA
}
rm(i)

e2 <- df


plot(log10(a2$natural_frequency),
     log10(a2$f_nat.cospec.w_u..cov.w_u.),type='l',
     ylim = c(-4,0.5), xlim = c(-3,1))
lines(log10(b2$natural_frequency),
      log10(b2$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'blue')
lines(log10(c2$natural_frequency),
      log10(c2$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'green')
lines(log10(d2$natural_frequency),
      log10(d2$f_nat.cospec.w_u..cov.w_u.),pch=19, col = 'yellow')
lines(log10(e2$natural_frequency),
      log10(e2$f_nat.cospec.w_u..cov.w_u.),pch=19, col = 'red')
