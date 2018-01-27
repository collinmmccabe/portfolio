install.packages("statnet")
library(statnet)

# set working directory
# setwd("<setme>")
epidat <- read.csv("epidemiology_data.csv")
event.times <- c(epidat$infection_begin, epidat$infection_end)
event.times <- unique(event.times)
event.times <- event.times[order(event.times)]
event.times <- c(-12, event.times, event.times[length(event.times)] + 12)
# write.csv(event.times, "eventtimes.csv")
cont <- read.csv("contact_data.csv")
cont.inf <- cont[which(cont$infection == TRUE),]

inf.status <- array(0, dim=c(nrow(epidat), length(event.times)))
for (i in 1:length(event.times))
{
	for (j in 1:nrow(epidat))
	{
		if (epidat$infection_begin[j] == event.times[i])
		{
			inf.status[j,i:length(event.times)] = 2
		}
		else 
		{
			if (epidat$infection_end[j] == event.times[i])
			{
				inf.status[j,i:length(event.times)] = 4
			}
		}
	}
}
inf.status[,length(event.times)] = 4
inf.status <- rbind(inf.status[1,1:length(event.times)],inf.status[10:19,1:length(event.times)],inf.status[2,1:length(event.times)],inf.status[20:29,1:length(event.times)],inf.status[3,1:length(event.times)],inf.status[30:38,1:length(event.times)],inf.status[4:9,1:length(event.times)])
# write.csv(inf.status, "infstatus.csv")

el=read.csv("undirHEBV.csv",header=FALSE) # read a .csv file
el[,1]=as.character(el[,1])
el[,2]=as.character(el[,2])
n = network(el,matrix.type="edgelist",directed=FALSE)
gcoord <- gplot.layout.kamadakawai(n, layout.par=NULL)
gplot(n, gmode="graph", coord=gcoord, vertex.col=1, cex=2, edge.col=1, displaylabels=FALSE)
text(mean(gcoord[,1]), max(gcoord[,2])+.5, "HEBV Transmission Network", cex=2, font=4)
text(mean(gcoord[,1]), min(gcoord[,2])-.5, paste("t =",round(event.times[1], digits=2),"hours from outbreak", sep=" "), cex=1.5)
legend(x=-4, y=-2, legend=c("Susceptible", "Infectious", "Recovered"), fill=c(0,2,4), border=1, cex=1)

#library(animation)
#ani.options(convert=shQuote("C:/Program Files/ImageMagick-6.7.9-Q16/convert.exe")) # tell R where to look for ImageMagick's convert.exe
#saveGIF(
  for (i in 1:length(event.times)){
	gplot(n, gmode="graph", coord=gcoord, vertex.col=inf.status[,i], cex=2, edge.col=1, displaylabels=FALSE)
	text(mean(gcoord[,1]), max(gcoord[,2])+.5, "HEBV Transmission Network", cex=2, font=4)
	text(mean(gcoord[,1]), min(gcoord[,2])-.5, paste("t =",round(event.times[i], digits=2),"hours from outbreak", sep=" "), cex=1.5)
	legend(x=-4, y=-2, legend=c("Susceptible", "Infectious", "Recovered"), fill=c(0,2,4), border=1, cex=1)
}
#, movie.name = "animation.gif", img.name = "Rplot", convert = "convert", outdir="c:/",cmd.fun = system, clean = TRUE)

