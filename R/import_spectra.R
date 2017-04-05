path <- 'data/spectra/E6_no_20130923-1600_binned_cospectra_2017-04-05T161542_adv.csv'
spectra <- read.csv(path, skip = 11, header = T)

df <- spectra

# Changing all the '-9999.0' or '-9999' (missing data) to NA
for (i in 1:length(df)){
  df[i][df[i] == '-9999' | df[i] == '-9999.0'] <- NA
}
rm(i)

e6no<- df

# Lateral with VAWT
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
# Lateral without VAWT
plot(log10(a2no$natural_frequency),
     log10(a2no$f_nat.cospec.w_u..cov.w_u.),type='l',
     ylim = c(-4,0.5), xlim = c(-3,1))
lines(log10(b2no$natural_frequency),
      log10(b2no$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'blue')
lines(log10(c2no$natural_frequency),
      log10(c2no$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'green')
lines(log10(d2no$natural_frequency),
      log10(d2no$f_nat.cospec.w_u..cov.w_u.),pch=19, col = 'yellow')
lines(log10(e2no$natural_frequency),
      log10(e2no$f_nat.cospec.w_u..cov.w_u.),pch=19, col = 'red')

# Vertical with VAWT
plot(log10(e2$natural_frequency),
     log10(e2$f_nat.cospec.w_u..cov.w_u.),type='l',
     ylim = c(-4,0.5), xlim = c(-3,1))
lines(log10(e3$natural_frequency),
      log10(e3$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'blue')
lines(log10(e4$natural_frequency),
      log10(e4$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'green')
lines(log10(e5$natural_frequency),
      log10(e5$f_nat.cospec.w_u..cov.w_u.),pch=19, col = 'yellow')
lines(log10(e6$natural_frequency),
      log10(e6$f_nat.cospec.w_u..cov.w_u.),pch=19, col = 'red')



# Vertical without VAWT
plot(log10(e2no$natural_frequency),
     log10(e2no$f_nat.cospec.w_u..cov.w_u.),type='l',
     ylim = c(-4,0.5), xlim = c(-3,1))
lines(log10(e3no$natural_frequency),
      log10(e3no$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'blue')
lines(log10(e4no$natural_frequency),
      log10(e4no$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'green')
lines(log10(e5no$natural_frequency),
      log10(e5no$f_nat.cospec.w_u..cov.w_u.),pch=19, col = 'yellow')
lines(log10(e6no$natural_frequency),
      log10(e6no$f_nat.cospec.w_u..cov.w_u.),pch=19, col = 'red')


# Comparison between with and without VAWT A2
plot(log10(a2$natural_frequency),
     log10(a2$f_nat.cospec.w_u..cov.w_u.),type='l',
     ylim = c(-4,0.5), xlim = c(-3,1))
lines(log10(a2no$natural_frequency),
      log10(a2no$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'blue')

# Comparison between with and without VAWT B2
plot(log10(b2$natural_frequency),
     log10(b2$f_nat.cospec.w_u..cov.w_u.),type='l',
     ylim = c(-4,0.5), xlim = c(-3,1))
lines(log10(b2no$natural_frequency),
      log10(b2no$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'blue')

# Comparison between with and without VAWT C2
plot(log10(c2$natural_frequency),
     log10(c2$f_nat.cospec.w_u..cov.w_u.),type='l',
     ylim = c(-4,0.5), xlim = c(-3,1))
lines(log10(c2no$natural_frequency),
      log10(c2no$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'blue')

# Comparison between with and without VAWT D2
plot(log10(d2$natural_frequency),
     log10(d2$f_nat.cospec.w_u..cov.w_u.),type='l',
     ylim = c(-4,0.5), xlim = c(-3,1))
lines(log10(d2no$natural_frequency),
      log10(d2no$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'blue')

# Comparison between with and without VAWT E2
plot(log10(e2$natural_frequency),
     log10(e2$f_nat.cospec.w_u..cov.w_u.),type='l',
     ylim = c(-4,0.5), xlim = c(-3,1))
lines(log10(e2no$natural_frequency),
      log10(e2no$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'blue')

# Comparison between with and without VAWT E3
plot(log10(e3$natural_frequency),
     log10(e3$f_nat.cospec.w_u..cov.w_u.),type='l',
     ylim = c(-4,0.5), xlim = c(-3,1))
lines(log10(e3no$natural_frequency),
      log10(e3no$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'blue')

# Comparison between with and without VAWT E4
plot(log10(e4$natural_frequency),
     log10(e4$f_nat.cospec.w_u..cov.w_u.),type='l',
     ylim = c(-4,0.5), xlim = c(-3,1))
lines(log10(e4no$natural_frequency),
      log10(e4no$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'blue')

# Comparison between with and without VAWT E5
plot(log10(e5$natural_frequency),
     log10(e5$f_nat.cospec.w_u..cov.w_u.),type='l',
     ylim = c(-4,0.5), xlim = c(-3,1))
lines(log10(e5no$natural_frequency),
      log10(e5no$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'blue')

# Comparison between with and without VAWT E6
plot(log10(e6$natural_frequency),
     log10(e6$f_nat.cospec.w_u..cov.w_u.),type='l',
     ylim = c(-4,0.5), xlim = c(-3,1))
lines(log10(e6no$natural_frequency),
      log10(e6no$f_nat.cospec.w_u..cov.w_u.),pch =19, col = 'blue')

