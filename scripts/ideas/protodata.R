### This script creates a data frame of fake data to illustrate what I think
### I will need to perform the clinical and cognitive analyses.
###
### Ellyn Butler
### July 16, 2020

library('ggplot2')

proto_df <- data.frame(bblid=c(rep(82071, 4), rep(87773, 3), rep(81277, 3)),
  age=c(12.5, 13, 14, 18, 16, 16.5, 17.2, 15, 17.2, 20),
  sex=c(rep('F', 4), rep('M', 3), rep('F', 3)),
  diagnosis=c('PS', 'OP', 'PS', 'PS', 'TD', 'OP', 'OP', 'OP', 'OP', 'TD'),
  first_diagnosis=c(rep('PS - First Diagnosis', 4), rep('TD - First Diagnosis', 3),
    rep('OP - First Diagnosis', 3)),
  last_diagnosis=c(rep('PS - Last Diagnosis', 4), rep('OP - Last Diagnosis', 3),
    rep('TD - Last Diagnosis', 3)),
  SIPS=c(3, 2, 3, 4, 1, 2, 2, 2, 2, 1), PMAT=c(-.5, -.2, -.1, -.12,
  .1, .67, 1.2, .5, .6, .72))

proto_df$first_diagnosis <- ordered(proto_df$first_diagnosis,
  c('PS - First Diagnosis', 'OP - First Diagnosis', 'TD - First Diagnosis'))
proto_df$last_diagnosis <- ordered(proto_df$last_diagnosis,
  c('TD - Last Diagnosis', 'OP - Last Diagnosis', 'PS - Last Diagnosis'))

write.csv(proto_df, '~/Documents/pncLongitudinalPsychosis/scripts/ideas/protodata.csv', row.names=FALSE)

############ Plot trajectories ############

sips_plot <- ggplot(proto_df, aes(x=age, y=SIPS, color=last_diagnosis)) +
  theme_linedraw() + geom_point() + geom_line() +
  facet_grid(first_diagnosis ~ last_diagnosis) +
  scale_color_manual(values=c('deepskyblue2', 'darkorchid1', 'firebrick4')) +
  theme(legend.position = 'none')

pdf(file='~/Documents/pncLongitudinalPsychosis/plots/illustration3x3.pdf', width=6, height=6)
sips_plot
dev.off()
