# size of capital grid
nk = 2^(4:14);

# times
cputimes = c(0.22,0.48,1.14,2.32,5.64,14.17,39.8,132.97,453.71,1732.68,6735.47);
gputimes = c(1.27,1.29,1.33,1.39,1.55,2.01,3.58,10.44,35.32,136.98,550.29);
ratio = gputimes/cputimes;

# regress sqrt(time) on capital grid length
cpureg = lm(sqrt(cputimes)~nk);
gpureg = lm(sqrt(gputimes)~nk);

# fine grid for capital
nkgrid = seq(min(nk),max(nk),length=10000);

# interpolate and extrapolate times
cpufit = (cpureg$coef[1] + cpureg$coef[2]*nkgrid)^2;
gpufit = (gpureg$coef[1] + gpureg$coef[2]*nkgrid)^2;

# plot
pdf(file="plottimes.pdf", height=8, width=10)
plot(nk,cputimes,xlab="N_k",ylab="Computation Times");
points(nk,gputimes);
lines(nkgrid,cpufit,col="blue");
lines(nkgrid,gpufit,col="red");
legend(2500,6000,legend=c("CPU","GPU"),lty=c(1,1),col=c("blue","red"))
dev.off()

# extrapolate
nk_x = 2^(15:20);
cputimes_x = (cpureg$coef[1] + cpureg$coef[2]*nk_x)^2;
gputimes_x = (gpureg$coef[1] + gpureg$coef[2]*nk_x)^2;
ratio_x = gputimes_x/cputimes_x;
