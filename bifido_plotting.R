# Bifido plot 

df1 <- "taxatable-filtered-relative.tsv"
df_table <- read.delim(df, sep = "\t", 
                       skip = 1, check.names = FALSE)
df_split <- df_table %>% separate(tax,into = c('K','P','C','O','F','G','S','T'),sep = ";")

df_split_cleaned <- subset(df_split,select = -c(K,P,C,O,F,G))
df_split_cleaned$Taxonomy <-paste(df_split_cleaned$S,"_",df_split_cleaned$T)
df_split_cleaned <- subset(df_split_cleaned,select = -c(S,T))
df_split_cleaned<-df_split_cleaned[,c(9, 1:8)]

##TODO add other for Genus


## plot on the x put samples , Y put numbers
2 plots one with the constant bifido and the other with the constant depth

want to tell them we can detect strain if at least this % of the total depth if have a depth of this 
-----------
  make heat map -> where bifido is its self, and all bifidos are seen in heatmap, but everything else is placed in other



sample <- df_try_cut["k__Bacteria;p__Actinobacteria;c__Actinobacteria;o__Bifidobacteriales;f__Bifidobacteriaceae;g__Bifidobacterium;s__Bifidobacterium_longum;t__Bifidobacterium_longum_ATCC_15697",]
t_sample <- t(sample)
# bad code I know but test is the name of the adta frame that I used
test <- as.data.frame(t_sample)
test$names <- rownames(test)
colnames(test)[1]  <- "Bifido_15697"
row.names(test) <- NULL
test<-test[,c(2, 1:1)]
# THis plot the labels are not showing
plot(test$Bifido_15697, pch=20, ylim=c(0, 600),xaxt="n", xlab="Sample", ylab="Reads Mapped")
 # Plot the axis separately
axis(1, at=1:8, labels=test$names)

test$x <- c(100,1000,10000,100000,500000,1,1,1)
test$x2 <- log(test$x)
test$y2 <- log(test$Bifido_15697)

plot_1 <- ggplot(test[1:5,],aes(x = names, y = Bifido_15697)) + geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
png("bifido_seqs.png")
print(plot_1)
dev.off()

plot_2 <- ggplot(test[6:8,],aes(x = names, y = Bifido_15697)) + geom_point() +
theme(axis.text.x = element_text(angle = 90, hjust = 1))
png("bifido_depths_cont.png")
print(plot_2)
dev.off()

q <- ggplot(test[1:5,],aes(x = x2, y = y2)) + geom_point() +
 theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
xlab("Number of Bifido Sequences") + ylab("Log reads mapped")

png("bifido_seqs_logged.png")
print(q)
dev.off()

# plot where x values is numeric 


# divide into two groups 


# maing my heat map
library("stringr")   

df_bifido <- df_table[str_detect(df_table$tax, "Bifidobacterium"), ] 
df_other <- anti_join(df_table,df_bifido)
df_other <- rbind(df_other, data.frame(tax = "other", t(colSums(df_other[, -1]))))
# make column names cleaner
df_bifido_c <- df_bifido %>% separate(tax,into = c('K','P','C','O','F','G','S','T'),sep = ";")
df_bifido_c <- subset(df_bifido_c,select = -c(K,P,C,O,F,G))
df_bifido_c$tax <-paste(df_bifido_c$S,"_",df_bifido_c$T)
df_bifido_c <- subset(df_bifido_c,select = -c(S,T))
df_bifido_c<-df_bifido_c[,c(9, 1:8)]


df_bifido_other <- rbind(df_bifido,df_other[1059,])
df_bifido_other_row <- df_bifido_other
rownames(df_bifido_other_row) <- df_bifido_other_row[,1]
df_bifido_other_row[,1] <- NULL
df_bifido_other_row_trim <- df_bifido_other_row[,1:5]
df_bifido_other_row_trim_matrix <- data.matrix(df_bifido_other_row_trim)

png("rel_heatmap_bifido_other.png")
heatmap<-heatmap.2(df_bifido_other_row_trim_matrix,col= brewer.pal(9,"RdYlBu"),scale="row",Rowv=NA,Colv=NA,trace="none",density.info="none")
dev.off()


df_bifido_c_matrix <- df_bifido_c
rownames(df_bifido_c_matrix) <- df_bifido_c_matrix[,1]
df_bifido_c_matrix[,1] <- NULL
df_bifido_c_matrix <- df_bifido_c_matrix[,1:5]
df_bifido_c_matrix <- data.matrix(df_bifido_c_matrix)
---------------------
single_column 
rownames(single_column) <- single_column[,1]
single_column[,1] <- NULL
single_column <- single_column[,1]



df1 <- "taxatable-filtered-relative.tsv"
df_rel <- read.delim(df1, sep = "\t", 
                     +                         check.names = FALSE)
df_rel_bifido <- df_rel[str_detect(df_rel$tax, "Bifidobacterium"), ] 
df_rel_bifido <- df_rel_bifido %>% separate(tax,into = c('K','P','C','O','F','G','S','T'),sep = ";")
df_rel_bifido <- subset(df_rel_bifido,select = -c(K,P,C,O,F,G))
df_rel_bifido$tax <-paste(df_rel_bifido$S,"_",df_rel_bifido$T)
df_rel_bifido <- subset(df_rel_bifido,select = -c(S,T))
df_rel_bifido<-df_rel_bifido[,c(9, 1:8)]


#df_rel_bifido <- df_rel_bifido
rownames(df_rel_bifido) <- df_rel_bifido[,1]
df_rel_bifido[,1] <- NULL
df_rel_bifido_m <- df_rel_bifido[,1:5]
df_rel_bifido_m <- data.matrix(df_rel_bifido_m)
png("rel_heatmap_bifido.png")
heatmap_2<-heatmap.2(df_rel_bifido_m,col= brewer.pal(9,"RdYlBu"),scale="row",Rowv=NA,Colv=NA,trace="none",density.info="none",margins=c(12,16))
dev.off()
