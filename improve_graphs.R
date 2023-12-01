# Improve cluster line-graphs
library(ggplot2)
raw <- read.csv("/home/harry/github/asp_lncrna/v11.5.1/test1_full_run/itra_k20.longdata.csv", row.names = 1)
# Make a condition zero variable with value set to zero
cond_zero <- raw[raw$condition == 0.25,]
cond_zero$value <- rep(0, nrow(cond_zero))
cond_zero$condition <- rep(0, nrow(cond_zero))
data <- rbind(cond_zero, raw)




lnc_data <- data[which(data$type == "lncrna"),]
pcg_data <- data[which(data$type == "pcg"),]

# Get the stat values per cluster

in_data <- data

all_clust_names <- c()
all_clust_cond <- c()
all_clust_stat <- c()
all_clust_sd <- c()
for(clust in unique(in_data$cluster)){
  # Extract data per cluster
  clust_data <- in_data[in_data$cluster == clust,]
  cond_res <- sort(unique(clust_data$condition))
  all_cond_stat <- c()
  all_cond_sd <- c()
  for(cond in cond_res){
    # Get stat
    cond_stat <- median(clust_data[clust_data$condition == cond,]$value)
    # Get X * STD
    cond_sd <- sd(clust_data[clust_data$condition == cond,]$value) * 3
    # Store information, per condition
    all_cond_sd <- c(all_cond_sd, cond_sd)
    all_cond_stat <- c(all_cond_stat, cond_stat)
  }
  # Store all of the information, per cluster
  clust_names <- rep(clust, length(all_cond_sd))
  all_clust_names <- c(all_clust_names, clust_names)
  all_clust_cond <- c(all_clust_cond, cond_res)
  all_clust_stat <- c( all_clust_stat, all_cond_stat)
  all_clust_sd <- c(all_clust_sd, all_cond_sd)
}
# Convert S.D NA values to zero
all_clust_sd[is.na(all_clust_sd)] <- 0
# Generate negative SD
sd_neg <- all_clust_sd * -1
# Generate appropriate dataframes
stat_df <- data.frame("cluster"=all_clust_names,
                        "condition"=all_clust_cond,
                        "value"=all_clust_stat,
                        "lower"=sd_neg,
                        "upper"=all_clust_sd)

stat_df$lower <- stat_df$lower + stat_df$value
stat_df$upper <- stat_df$upper + stat_df$value
ggplot_stat <- cbind(rep("stat", nrow(stat_df)), rep("stat", nrow(stat_df)), stat_df)
colnames(ggplot_stat)[1:2] <- c("transcript","type")

ggplot_data <- data

# Order based on type
ggplot_data <-  ggplot_data[order(ggplot_data$transcript, decreasing = F),]
# Convert type to factor
ggplot_data$transcript <- factor(ggplot_data$transcript, levels=(sort(unique(ggplot_data$transcript), decreasing = T)))

p <- ggplot(data = ggplot_data) +
  geom_line(aes(x=condition, y=value, group=transcript, colour=type), size =1.1) +
  scale_colour_manual(values=c("#5E5FF6","#F77373"), breaks=levels(ggplot_data$type)) +
  geom_hline(yintercept=0) +
  geom_line(data=ggplot_stat, aes(x=condition, y=value, group=transcript), color="black", size = 1.2) +
  geom_ribbon(data=ggplot_stat, aes(x=condition, ymin=lower, ymax=upper), alpha = 0.2, color = "#6E6E6E") +
  facet_wrap(~ cluster, ncol = 5) +
  xlab("xMIC") + ylab(expression(Standardised~log["2"]-fold~change)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 15),
        panel.spacing.x = unit(2, "lines"))

print(p)


## Barplots for number of lncRNA per cluster
# Get frequencies

cluster_list <- c()
transcript_list <- c()
total_list <- c()
percent_list <- c()
raw_single <- subset(raw, condition==0.25)
for(c in unique(raw_single[,"cluster"])){
  c_tot <- nrow(subset(raw_single, cluster==c))
  for(t in unique(raw_single[,"type"])){
    tot <- nrow(subset(raw_single, type==t & cluster==c))
    perc <- (tot/c_tot) * 100
    cluster_list <- c(cluster_list, c)
    transcript_list <- c(transcript_list, t)
    total_list <- c(total_list, tot)
    percent_list <- c(percent_list, perc)
  }
}

df_bar <- data.frame(transcript=transcript_list, cluster=cluster_list, 
                     total=total_list, percentage=percent_list)

 ggplot(data=df_bar, aes(x=cluster, y=percentage, fill=transcript)) +
                      geom_bar(stat="identity") +
                      scale_fill_manual(values=c("#5E5FF6","#F77373")) +
                      scale_x_continuous(breaks=1:max(raw_single$cluster)) +
                      xlab("Cluster") + ylab("Percentage") +
                      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                            panel.background = element_blank(), axis.line = element_line(colour = "black"),
                            axis.text.x = element_text(size = 15, color="black"),
                            axis.text.y = element_text(size = 15, color="black"),
                            axis.title = element_text(size = 20),
                            strip.text = element_text(size = 15),
                            legend.position = "none",
                            panel.spacing.x = unit(2, "lines"))

ggplot(data=df_bar, aes(x=cluster, y=total, fill=transcript)) +
                  geom_bar(stat="identity", position=position_dodge()) +
                  scale_fill_manual(values=c("#5E5FF6","#F77373")) +
                  scale_x_continuous(breaks=1:max(raw_single$cluster)) +
                  xlab("Cluster") + ylab("Total") +
                  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        panel.background = element_blank(), axis.line = element_line(colour = "black"),
                        axis.text.x = element_text(size = 15, color="black"),
                        axis.text.y = element_text(size = 15, color="black"),
                        axis.title = element_text(size = 20),
                        strip.text = element_text(size = 15),
                        legend.position = "none",
                        panel.spacing.x = unit(2, "lines"))


# ggplot(dat) + 
#   geom_line(aes(x = x, y = y, color = id2), size = 2) + 
#   # use your dummy column for the line order
#   scale_colour_manual('id',values = cols, breaks = levels(dat$id))  
# 
# 
# as.character(data$cluster)
# ggplot(data=stat_df, aes(x=condition, y=values)) +
#   geom_line(color="#00BFC4") + 
#   geom_hline(yintercept=0) +
#   geom_ribbon(aes(ymin =lower, ymax = upper), alpha = 0.2) +
#   facet_wrap(~ cluster) +
#   xlab("xMIC") + ylab(expression(Standardised~log["2"]-fold~change)) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# 
# p <- ggplot(data=MULTIlong_data, aes(x=condition, y=value, group=transcript)) +
#   geom_line(color="#00BFC4") + 
#   geom_line(data=MULTIstat_long, aes(x=condition, y=value, group=transcript), color="red") +
#   geom_hline(yintercept=0) +
#   facet_wrap(~ cluster) +
#   xlab("MIC") + ylab("Standardized Log2Fold Change")

