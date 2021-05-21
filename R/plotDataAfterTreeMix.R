args <- commandArgs(trailingOnly = TRUE)

if (1 == 0) {

    args <- c(
        "/well/myers/rwdavies/primates/treemix",
        "~/personal/proj/primates/R/treemix_plotting_funcs.R",
        "felidae.csbjacessltj.GATKug.treemix.migrants.0.out",
        "felidae.csbjacessltj.GATKug.treemix.frq.gz"
    )
    
}

print(args)
dir <- args[1]
treemix_script <- args[2]
stem <- args[3]
input <- args[4]

    

setwd(dir)
source(treemix_script)


T <- TRUE

check_file <- function(f) {
    if (file.exists(f) == FALSE) {
        stop(paste0("Cannot find file:", f))
    }
}
check_file(paste0(stem, ".vertices.gz"))
check_file(paste0(stem, ".edges.gz"))
check_file(paste0(stem, ".covse.gz"))
check_file(input)

remove_felis6.2 <- function(x) return(gsub("_FelisCatus6.2", "", x))

### new version of this function - fixed loading of o - turn into matrix
plot_resid_mod <- function(stem, input,
         min = -0.009,
         max = 0.009,
         cex = 1,
         usemax = T,
         wcols = "r"
         )
{
	c = read.table(gzfile(paste(stem, ".cov.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
	m = read.table(gzfile(paste(stem, ".modelcov.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
        rownames(c) <- remove_felis6.2(rownames(c))
        colnames(c) <- remove_felis6.2(colnames(c))
        rownames(m) <- remove_felis6.2(rownames(m))
        colnames(m) <- remove_felis6.2(colnames(m))    
	names(c) = rownames(c)
	names(m) = rownames(m)
        h <- read.table(input, as.is = T, comment.char = "", quote = "", nrow = 1)
	o = matrix(remove_felis6.2(h),ncol=1)        
	##o = matrix(read.table(pop_order, as.is = T, comment.char = "", quote = ""),ncol=1)
	se = read.table(gzfile(paste(stem, ".covse.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
        rownames(se) <- remove_felis6.2(rownames(se))
        colnames(se) <- remove_felis6.2(colnames(se))
	mse = apply(se, 1, mean)
	mse = mean(mse)
	print(mse)	
	c = c[order(names(c)), order(names(c))]
	m = m[order(names(m)), order(names(m))]
	tmp = c -m 
	#tmp = m - c
	#tmp = (m-c)/m
	#print(tmp)
	toplot = data.frame(matrix(nrow = nrow(tmp), ncol = ncol(tmp)))
	for(i in 1:nrow(o)){
	        for( j in 1:nrow(o)){
			#print(paste(o[i,1], o[j,1]))
			if (o[i,1] %in% names(tmp) ==F){
				print(paste("not found", o[i,1]))
			}
			if (o[j,1] %in% names(tmp) ==F){
				print(paste("not found", o[j,1]))
			}
        	        toplot[i, j] = tmp[which(names(tmp)==o[i,1]), which(names(tmp)==o[j,1])]
        	}
	}
	#print(toplot)
	if (usemax){
		m1 = max(abs(toplot), na.rm = T)
		max = m1*1.02
		min = -(m1*1.02)	
	}
	print("here")
	names(toplot) = o[,1]
	toreturn = plot_resid_internal(toplot, max = max, min = min, wcols = wcols, mse = mse, o = o, cex = cex)
	return(toreturn)
}

plot_tree=function(stem, o = NA, cex = 1, disp = 0.003, plus = 0.01, flip = vector(), arrow = 0.05, scale = T, ybar = 0.1, mbar = T, plotmig = T, plotnames = T, xmin = 0, lwd = 1, font = 1){
	d = paste(stem, ".vertices.gz", sep = "")
	e = paste(stem, ".edges.gz", sep = "")
	se = paste(stem, ".covse.gz", sep = "")
	d = read.table(gzfile(d), as.is = T, comment.char = "", quote = "")
	e = read.table(gzfile(e), as.is  = T, comment.char = "", quote = "")
	if (!is.na(o)){
		o = read.table(o, as.is = T, comment.char = "", quote = "")
	}
	e[,3] = e[,3]*e[,4]
	e[,3] = e[,3]*e[,4]
	se = read.table(gzfile(se), as.is = T, comment.char = "", quote = "")
	m1 = apply(se, 1, mean)
	m = mean(m1)
	#m = 0
	for(i in 1:length(flip)){
		d = flip_node(d, flip[i])
	}
	d$x = "NA"
	d$y = "NA"
	d$ymin = "NA"
	d$ymax = "NA"
	d$x = as.numeric(d$x)
	d$y = as.numeric(d$y)
	d$ymin = as.numeric(d$ymin)
	d$ymax = as.numeric(d$ymax)
	d = set_y_coords(d)
	d = set_x_coords(d, e)
	print(d)
	d = set_mig_coords(d, e)
	plot_tree_internal(d, e, o = o, cex = cex, xmin = xmin, disp = disp, plus = plus, arrow = arrow, ybar = ybar, mbar = mbar, mse = m, scale = scale, plotmig = plotmig, plotnames = plotnames, lwd = lwd, font = font)
	return(list( d= d, e = e))
}



plot_tree_internal=function(d, e, o = NA, cex = 1, disp = 0.005, plus = 0.005, arrow = 0.05, ybar = 0.01, scale = T, mbar = T, mse = 0.01, plotmig = T, plotnames = T, xmin = 0, lwd = 1, font = 1){
	plot(d$x, d$y, axes = F, ylab = "", xlab = "Drift parameter", xlim = c(xmin, max(d$x)+plus), pch = "")
	axis(1)
	mw = max(e[e[,5]=="MIG",4])
	mcols = rev(heat.colors(150))
	for(i in 1:nrow(e)){
		col = "black"
		if (e[i,5] == "MIG"){
			w = floor(e[i,4]*200)+50
			if (mw > 0.5){
				w = floor(e[i,4]*100)+50
			}
			col = mcols[w]
			if (is.na(col)){
				col = "blue"
			}
		}
		v1 = d[d[,1] == e[i,1],]
		v2 = d[d[,1] == e[i,2],]
		if (e[i,5] == "MIG"){
			if (plotmig){
			arrows( v1[1,]$x, v1[1,]$y, v2[1,]$x, v2[1,]$y, col = col, length = arrow)
			}
		}
		else{
			lines( c(v1[1,]$x, v2[1,]$x), c(v1[1,]$y, v2[1,]$y), col = col, lwd = lwd)
		}
	}
        ### HACK
        ##print("APPLY NAME HACK")
        ##d[d[,5]=="TIP",2]=newNames[match(d[d[,5] == "TIP",2],oriNames)]
        d[d[, 5] == "TIP", 2] <- gsub("_FelisCatus6.2", "", d[d[, 5] == "TIP", 2])
        ##print(head(d)        )
        ##print("DONE HACK")
        ### HACK
	tmp = d[d[,5] == "TIP",]
        ##print("DONE HACK")        
	print(tmp[,2])
	print(disp)
	if ( !is.na(o)){
		for(i in 1:nrow(tmp)){
			tcol = o[o[,1] == tmp[i,2],2]
			if(plotnames){
				print(tmp[i,2])
				text(tmp[i,]$x+disp, tmp[i,]$y, labels = tmp[i,2], adj = 0, cex = cex, col  = tcol, font = font)
			}
		}
	}
	else{
		if (plotnames){
                  print("printing names")
                  print(tmp[,2])
		text(tmp$x+disp, tmp$y, labels = tmp[,2], adj = 0, cex = cex, font = font)
		}
	}
	if (scale){
	#print (paste("mse", mse))
        lines(c(0, mse*10), c(ybar, ybar))
	text( 0, ybar - 0.04, lab = "10 s.e.", adj = 0, cex  = 0.8)
	lines( c(0, 0), c( ybar - 0.01, ybar+0.01))
	lines( c(mse*10, mse*10), c(ybar- 0.01, ybar+ 0.01))
	}
        if (mbar){
                mcols = rev( heat.colors(150) )
                mcols = mcols[50:length(mcols)]
                ymi = ybar+0.15
                yma = ybar+0.35
                l = 0.2
                w = l/100
                xma = max(d$x/20)
                rect( rep(0, 100), ymi+(0:99)*w, rep(xma, 100), ymi+(1:100)*w, col = mcols, border = mcols)
                text(xma+disp, ymi, lab = "0", adj = 0, cex = 0.7)
		if ( mw >0.5){ text(xma+disp, yma, lab = "1", adj = 0, cex = 0.7)}
                else{
			text(xma+disp, yma, lab = "0.5", adj = 0, cex =0.7)
		}
		text(0, yma+0.06, lab = "Migration", adj = 0 , cex = 0.6)
		text(0, yma+0.03, lab = "weight", adj = 0 , cex = 0.6)
        }	
}


## get names
file <- paste0(stem, ".pdf")
pdf(file,height=6,width=6)
plot_tree(stem,plus=0.15)
print(file)
dev.off()

file <- paste0(stem, ".both.pdf")
pdf(file,height=4,width=8)
par(mfrow=c(1,2))
par(oma=c(3,3,0,0))
par(mar=c(0,0,0,0))
plot_tree(stem,plus=0.15)
par(mar=c(3,3,0,0))
plot_resid_mod(stem, input = input)
print(file)
dev.off()




