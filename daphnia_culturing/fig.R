mort = read.csv("C:/Users/Mike/git/seven_lakes/FA/daphnia_culturing/mortality_data.csv", skip=7)[,4:6]

pdf(width=7, height=6, file='D:/Dropbox/Grad/Projects/Thesis/Seven Lakes Project 2014/Manuscript/figures/fig6.pdf',
    compress=FALSE)
barplot(t(as.matrix(mort[2:3])), beside=TRUE, col=c('white','black'), width=.5,
        names.arg=c('C. ozolinii','Mont.','Mont. pr.','Alpine','Alpine pr.'),
        space=c(0,.5), ylim=c(0,80), ylab='# individuals', xlab='Treatment', yaxt='n')
axis(2, seq(0,80,10), las=2)
legend('topleft', legend=c('Total Mortality','Total Harvest'), bty='n', fill=c('white','black'))
dev.off()
shell('D:/Dropbox/Grad/Projects/Thesis/Seven^ Lakes^ Project^ 2014/Manuscript/figures/fig6.pdf')
