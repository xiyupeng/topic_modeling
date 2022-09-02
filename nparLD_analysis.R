### 
# Author: Xiyu Peng
# stat analysis and visualization for single clusters.
#


library(tidyverse)
library(ggpubr)
library(patchwork)
library(ggplot2)

merge.data <- read.csv("~/flow_cytometry/lda/merge_data_X50_1.5_new.csv")
rm_samples<-c("17-162-32C","17-162-34B","17-162-07B")  ## not contains enough cells
merge.data<-merge.data[!merge.data$sample %in% rm_samples,]
merge.data %>% select(starts_with("cluster")) %>% colnames()->colnames_data

merge.data$res<-factor(merge.data$res)
levels(merge.data$res)
levels(merge.data$res)<-c("PD/SD","CR/PR")

merge.data.reformat<-merge.data %>% gather(.,key = "cluster",value = "freq",colnames_data)

library(nparLD)
merge.data.sub.reformat<-merge.data.reformat %>% filter(pt %in% names(which(colSums(table(merge.data$time,merge.data$pt))==3)))
plot_clusters<-paste("cluster",0:19,sep = "")


### across imtypes
for (c in plot_clusters){
  re<-merge.data.sub.reformat %>% filter(cluster == c) %>% nparLD(freq ~ imtype*time, data = .,subject = "pt",description = FALSE)
  print(c)
  print(re)
}

for (c in plot_clusters){
  re<-merge.data.sub.reformat %>% filter(cluster == c) %>% nparLD(freq ~ res*time, data = .,subject = "pt",description = FALSE)
  print(c)
  print(re)
}

for (c in plot_clusters){
  re<-merge.data.sub.reformat %>% filter(cluster == c) %>% nparLD(freq ~ tox*time, data = .,subject = "pt",description = FALSE)
  print(c)
  print(re)
}


## organize re$ANOVA.test and re$ANOVA.test.mod.Box into tables

result<-list()

for (c in plot_clusters){
  
  ### immunotypes
  re<-merge.data.sub.reformat %>% filter(cluster == c) %>% nparLD(freq ~ imtype*time, data = .,subject = "pt",description = FALSE)
  imm_array<-c(array(t(re$ANOVA.test)),re$ANOVA.test.mod.Box)
  imm_header<-outer(rownames(t(re$ANOVA.test)),colnames(t(re$ANOVA.test)),paste, sep = ".")
  imm_header<-c(imm_header, outer(rownames(t(re$ANOVA.test.mod.Box)),colnames(t(re$ANOVA.test.mod.Box)),paste, sep = "."))
  imm_header[13]<-"p-value.mod.imtype"
  
  ### response
  re<-merge.data.sub.reformat %>% filter(cluster == c) %>% nparLD(freq ~ res*time, data = .,subject = "pt",description = FALSE)
  res_array<-c(array(t(re$ANOVA.test)),re$ANOVA.test.mod.Box)
  res_header<-outer(rownames(t(re$ANOVA.test)),colnames(t(re$ANOVA.test)),paste, sep = ".")
  res_header<-c(res_header, outer(rownames(t(re$ANOVA.test.mod.Box)),colnames(t(re$ANOVA.test.mod.Box)),paste, sep = "."))
  res_header[13]<-"p-value.mod.res"
  
  ### toxicity
  re<-merge.data.sub.reformat %>% filter(cluster == c) %>% nparLD(freq ~ tox*time, data = .,subject = "pt",description = FALSE)
  tox_array<-c(array(t(re$ANOVA.test)),re$ANOVA.test.mod.Box)
  tox_header<-outer(rownames(t(re$ANOVA.test)),colnames(t(re$ANOVA.test)),paste, sep = ".")
  tox_header<-c(tox_header, outer(rownames(t(re$ANOVA.test.mod.Box)),colnames(t(re$ANOVA.test.mod.Box)),paste, sep = "."))
  tox_header[13]<-"p-value.mod.tox"
  
  result[[c]]<-c(imm_array,res_array,tox_array)
  
}

resultdf<-t(as.data.frame(result))
colnames(resultdf)<-c(imm_header,res_header,tox_header)
resultdf<-as.data.frame(resultdf)

## adjust pvalue with selected columns
select_col<-c("p-value.mod.imtype","p-value.mod.res","p-value.mod.tox","p-value.imtype:time","p-value.res:time","p-value.tox:time")
adj_pvalues<-resultdf %>% 
  dplyr::select(select_col) %>% 
  mutate(across(select_col,p.adjust,method = "BH"))
colnames(adj_pvalues)<-paste0("adj.",colnames(adj_pvalues))
resultdf<-cbind(resultdf,adj_pvalues)

## write the result
write.csv(resultdf,file = "~/flow_cytometry/lda/single_cluster_stat_analysis.csv")


#### make boxplots
colors<-c("#D55E00", "#009E73", "#0072B2")
merge.data.reformat$week<-as.factor(merge.data.reformat$time)
levels(merge.data.reformat$week)<-c("0","3","6")
## immunotypes
plot_clusters<-paste("cluster",0:19,sep = "")
ps_imtype<-list()
for (c in plot_clusters){
  ## boxplot with data points connected with lines
  p<-merge.data.reformat %>% 
    mutate (freq = freq/100) %>%
    filter(cluster == c) %>%
    ggboxplot(.,x = "week", y = "freq", add = "point", fill = "imtype", palette = colors, facet.by = "imtype")+geom_line(aes(group = pt),color = 'gray')+theme(legend.position = "none")+ggtitle(c)+ylab("Percentage")
  ps_imtype[[c]]<-p
}

wrap_plots(ps_imtype[1:20] ,guides = 'collect')
ggsave("imtype_clusters.png",width = 20, height = 20,units = "in")



## response 
ps_res<-list()
for (c in plot_clusters){
  ## boxplot with data points connected with lines
  p<-merge.data.reformat %>% 
    mutate (freq = freq/100) %>%
    filter(cluster == c) %>%
    ggboxplot(.,x = "week", y = "freq", add = "point", fill = "res", palette = colors, facet.by = "res")+geom_line(aes(group = pt),color = 'gray')+theme(legend.position = "none")+ggtitle(c)+ylab("Percentage")
  ps_res[[c]]<-p
}

wrap_plots(ps_res[1:20] ,guides = 'collect')
ggsave("res_clusters.png",width = 16, height = 20,units = "in")


## tox 
ps_tox<-list()
for (c in plot_clusters){
  p<-merge.data.reformat %>% 
    mutate (freq = freq/100) %>%
    filter(cluster == c) %>%
    ggboxplot(.,x = "week", y = "freq", add = "point", fill = "tox", palette = colors, facet.by = "tox")+geom_line(aes(group = pt),color = 'gray')+theme(legend.position = "none")+ggtitle(c)+ylab("Percentage")
  ps_tox[[c]]<-p
}

wrap_plots(ps_tox[1:20] ,guides = 'collect')
ggsave("tox_clusters.png",width = 16, height = 20,units = "in")


