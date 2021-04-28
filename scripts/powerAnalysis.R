# Power Analysis
library(pwr)
# Hypothesis 1 - territory
# A - ANOSIM (nw = 10, nb = 241)

# A - ANOVA (F17 = 8, F18 = 9, F20 = 16, S20 = 12)
# 5 or 6 PCNM axes depending on season (doing both)
pwr.anova.test(k = 5, f = 0.3, sig.level = 0.05, power = 0.8)
pwr.anova.test(k = 6, f = 0.3, sig.level = 0.05, power = 0.8)
# looking for current power
pwr.anova.test(k = 5, f = 0.3, sig.level = 0.05, n = 8)
pwr.anova.test(k = 6, f = 0.3, sig.level = 0.05, n = 9)
pwr.anova.test(k = 6, f = 0.3, sig.level = 0.05, n = 16)
pwr.anova.test(k = 5, f = 0.3, sig.level = 0.05, n = 12)

# B - ANOVA (n = )
pwr.anova.test(k = 3, f = 0.3, sig.level = 0.05, power = 0.8)

# Hypothesis 2 - host factors
# A - ANOVA (F17 = 15, F18 = 14, F20 = 19, S20 = 34)
pwr.f2.test(u = 4, f2 = 0.15, power = 0.8, sig.level = 0.05)
# 11, 10, 8, 6
pwr.anova.test(k = 11, f = 0.3, sig.level = 0.05, power = 0.8)
pwr.anova.test(k = 10, f = 0.3, sig.level = 0.05, power = 0.8)
pwr.anova.test(k = 8, f = 0.3, sig.level = 0.05, power = 0.8)
pwr.anova.test(k = 6, f = 0.3, sig.level = 0.05, power = 0.8)

pwr.anova.test(k = 11, f = 0.3, sig.level = 0.05, n = 15)
pwr.anova.test(k = 10, f = 0.3, sig.level = 0.05, n = 14)
pwr.anova.test(k = 8, f = 0.3, sig.level = 0.05, n = 19)
pwr.anova.test(k = 6, f = 0.3, sig.level = 0.05, n = 34)

# Hypothesis 3 - diet
# A - Linear model, ANOVA (n = 39)
# Breeder n = 2, 2, 1, 4, 5
# Non-breeder n = 2, 1, 15, 4, 3
pwr.anova.test(k = 5, f = 0.3, sig.level = 0.05, power = 0.8)
pwr.f2.test(u = 1, f2 = 0.15, power = 0.8, sig.level = 0.05)

# B - Two-sample t-test
pwr.t.test(sig.level = 0.05, power = 0.8, d=0.3, type = "two.sample",
           alternative = "greater") # getting an error with less but works with greater

# hypothesis 4 - parent
# A - Paired t-test (n = )
pwr.t.test(sig.level = 0.05, power = 0.8, d=0.3, type = "paired",
           alternative = "greater") # getting an error with less but works with greater

# B - Two sample t-tests (n = )
pwr.t.test(sig.level = 0.05, power = 0.8, d=0.3, type = "two.sample",
           alternative = "greater") # getting an error with less but works with greater

# hypothesis 5 - dispersal
# A - Linear model, ANOVA (n = 5,5,6)
pwr.f2.test(u = 1, f2 = 0.15, power = 0.8, sig.level = 0.05)
