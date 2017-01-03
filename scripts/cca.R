library(mixOmics)

otu_table <- read.table('../data/deblur.txt', skip=1, header=TRUE, 
                        comment.char='', sep='\t')
ms_table <- read.table('../data/MS1_features.txt', skip=1, header=TRUE, 
                       comment.char='', sep='\t')

otu_table <- otu_table[, order(names(otu_table))]
ms_table <- ms_table[, order(names(ms_table))] 

ids <- intersect(colnames(ms_table), colnames(otu_table))
otu_table <- otu_table[ids]
ms_table <- ms_table[ids]

ms_table.T <- t(ms_table[,2:ncol(ms_table)])
# Set the column headings
colnames(ms_table.T) <- ms_table[,1]

otu_table.T <- t(otu_table[,2:ncol(otu_table)])
# Set the column headings
colnames(otu_table.T) <- otu_table[,1]

X <- ms_table.T
Y <- otu_table.T

X <- X / rowSums(X)
Y <- Y / rowSums(Y)

grid1 <- seq(0.000001, 1, length = 51)
grid2 <- seq(0.000001, 1, length = 51)
## Depending on how fine your grid is, tune.rcc may take some time
cv.score <- tune.rcc(Y, X, grid1 = grid1, grid2 = grid2)

pdf("../results/score.pdf",width=6,height=4,paper='special') 
image(cv.score)
dev.off()

levels <- c(0.66, 0.7, 0.8, 0.88, 0.92, 0.94)
pdf("../results/contour.pdf",width=6,height=4,paper='special') 
contour(cv.score$grid1, cv.score$grid2, cv.score$mat, add = T, levels = levels,
        col = "blue")
dev.off()

result = rcc(Y, X, lambda1=0.04, lambda2=0.008)

pdf("../results/scree.pdf",width=6,height=4,paper='special') 
plot(result, scree.type = "barplot")
dev.off()


pdf("../results/correlation_plot.pdf",width=6,height=4,paper='special') 
plotVar(result, comp = 1:2, cutoff = 0.01, cex = c(5, 5))
dev.off()

pdf("../results/dendrogram.pdf",width=6,height=4,paper='special') 
cim(result, comp = 1:2, xlab = "microbes", ylab = "metabolites", margins = c(5,6))
dev.off()

