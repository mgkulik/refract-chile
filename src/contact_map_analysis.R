library(pacman)
p_load(dplyr, tidyr, tibble, ggplot2, ggpubr, seqinr, stringr, xlsx, mgsub)

setwd("../")
comp_path <- getwd()
data_path <- paste0(comp_path, '/refract-data/')
docs_path <- paste0(data_path, 'drSasa/structure.DNA_vs_PROTEIN.by_res.tsv')
plots_folder <- paste0(comp_path, "plots/")
cov_folder <- paste0(data_path, "covariance/")

SAVE_OUTPUT <- FALSE

start_posC <- 427
end_posB <- 424
start_pos <- 427

sigma <- 1/50

df_cont <- read.csv(paste0(comp_path, docs_path), sep="\t", stringsAsFactors=FALSE)

target_helix3 <- c(40,43,44)
target_helix6 <- c(90,91,94)
nt_ori_ord <- c("A", "T", "C", "G")
aa_order <- c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')

######### GETTING TARGET CONTACTS #########

df_contacts0 <- df_cont%>%
  separate(DNA..PROTEIN, into=c("aa", NA, "pos"))%>%
  mutate(X=rowSums(across(where(is.numeric))),
         pos=as.numeric(pos),
         aa=a(str_to_title(aa)))%>%
  filter(X>0)%>%
  select(-X)%>%
  gather(dna, vals, -aa, -pos)%>%
  separate(dna, into=c("nt", "chain", "pos_nt"))%>%
  mutate(nt=gsub("D", "", nt),
         pos_nt=as.numeric(pos_nt))%>%
  rename(pos_aa=pos)%>%
  #mutate(lg_vals=log10(vals))%>%
  filter(pos_aa%in%target_helix3|pos_aa%in%target_helix6)

df_contacts <- df_contacts0%>%
  mutate(helix=ifelse(pos_aa%in%target_helix3, "h3", "h6"))%>%
  filter(vals>0)%>%
  arrange(pos_aa, desc(vals))

df_contacts_top <- df_contacts%>%
  group_by(pos_aa)%>%
  slice_head(n=4)%>%
  group_by(pos_aa)%>%
  mutate(order = dplyr::row_number())

#chainB <- rev(427:448-24)
#chainC <- 427:448
# Selected h3: h3_40_420, h3_40_418, h3_43_417, h3_43_418, h3_44_419, h3_44_422
# Selected h6: h6_90_407, h6_90_408, h6_91_412, h6_91_411, h6_94_412, h6_94_408

df_contacts_AllB <- df_contacts%>%
  group_by(pos_aa)%>%
  mutate(order = dplyr::row_number(),
         pos_nt=ifelse(chain=="C", pos_nt-((pos_nt-start_posC+1)*2+1), pos_nt),
         nt=ifelse(chain=="C", ifelse(nt=="A", "T",
                                      ifelse(nt=="T", "A",
                                             ifelse(nt=="C", "G", "C"))), nt))%>%
  select(-chain)

df_contacts_AllC <- df_contacts%>%
  group_by(pos_aa)%>%
  mutate(order = dplyr::row_number(),
         pos_nt=ifelse(chain=="B", pos_nt+(end_posB-pos_nt+1)*2+1, pos_nt),
         nt=ifelse(chain=="B", ifelse(nt=="A", "T",
                                      ifelse(nt=="T", "A",
                                             ifelse(nt=="C", "G", "C"))), nt))%>%
  select(-chain)

if (SAVE_OUTPUT) {

  #------------------- HELIX 3 Plots -------------------
  
  #pall <- c("#FFFFFF", rev(heat.colors(3)))
  pall <- c("#FFFFFF", "#DDE8CB", "#00A933")
  
  df_contacts_topC_h3 <- df_contacts0%>%
    filter(chain=="C")%>%
    filter(pos_aa%in%target_helix3)%>%
    select(-chain)%>%
    unite(aa, pos_aa:aa, sep="/")%>%
    unite(nt, pos_nt:nt, sep="/")%>%
    mutate(lbls=ifelse(vals==0, "", round(vals, 2)))
  
  ord <- rev(unique(df_contacts_topC_h3$aa))
  
  plt_h3_C <- ggplot(df_contacts_topC_h3, aes(x = nt, y = ordered(aa, levels=!!ord), fill = vals, label=lbls)) +
    scale_fill_gradientn(colours=pall)+
    geom_tile()+
    geom_text(size=5)+
    theme_light()+
    scale_x_discrete(position = "top", guide=guide_axis(n.dodge=2))+
    theme(legend.position="none",
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  df_contacts_topB_h3 <- df_contacts0%>%
    filter(chain=="B")%>%
    filter(pos_aa%in%target_helix3)%>%
    select(-chain)%>%
    unite(aa, pos_aa:aa, sep="/")%>%
    unite(nt, pos_nt:nt, sep="/")%>%
    mutate(lbls=ifelse(vals==0, "", round(vals, 2)))%>%
    arrange(desc(nt))
  df_contacts_topB_h3 <- df_contacts_topB_h3[, c('aa', rev(colnames(df_contacts_topB_h3)[2:ncol(df_contacts_topB_h3)]))]
  df_contacts_topB_h3 <- df_contacts_topB_h3%>%mutate(nt=factor(nt, levels=!!unique(df_contacts_topB_h3$nt)))
  
  plt_h3_B <- ggplot(df_contacts_topB_h3, aes(x = nt, y = ordered(aa, levels=!!ord), fill = vals, label=lbls)) +
    scale_fill_gradientn(colours=pall)+
    geom_tile()+
    geom_text(size=5)+
    theme_light()+
    scale_x_discrete(guide=guide_axis(n.dodge=2))+
    theme(legend.position="none",
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  plt_h3 <- ggarrange(plt_h3_C, plt_h3_B, nrow=2)+bgcolor("#ffffff")
  #ggplot2::ggsave("heat_pos_h3.png", plt_h3, path = plots_folder, units = "in", width = 16, height = 3, dpi = 300, bg="white")
  
  #------------------- HELIX 6 Plots -------------------
  df_contacts_topC_h6 <- df_contacts0%>%
    filter(chain=="C")%>%
    filter(pos_aa%in%target_helix6)%>%
    select(-chain)%>%
    unite(aa, pos_aa:aa, sep="/")%>%
    unite(nt, pos_nt:nt, sep="/")%>%
    mutate(lbls=ifelse(vals==0, "", round(vals, 2)))
  
  ord <- rev(unique(df_contacts_topC_h6$aa))
  
  plt_h6_C <- ggplot(df_contacts_topC_h6, aes(x = nt, y = ordered(aa, levels=!!ord), fill = vals, label=lbls)) +
    scale_fill_gradientn(colours=pall)+
    geom_tile()+
    geom_text(size=5)+
    theme_light()+
    scale_x_discrete(position = "top", guide=guide_axis(n.dodge=2))+
    theme(legend.position="none",
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  df_contacts_topB_h6 <- df_contacts0%>%
    filter(chain=="B")%>%
    filter(pos_aa%in%target_helix6)%>%
    select(-chain)%>%
    unite(aa, pos_aa:aa, sep="/")%>%
    unite(nt, pos_nt:nt, sep="/")%>%
    mutate(lbls=ifelse(vals==0, "", round(vals, 2)))%>%
    arrange(desc(nt))
  df_contacts_topB_h6 <- df_contacts_topB_h6[, c('aa', rev(colnames(df_contacts_topB_h6)[2:ncol(df_contacts_topB_h6)]))]
  df_contacts_topB_h6 <- df_contacts_topB_h6%>%mutate(nt=factor(nt, levels=!!unique(df_contacts_topB_h6$nt)))
  
  
  plt_h6_B <- ggplot(df_contacts_topB_h6, aes(x = nt, y = ordered(aa, levels=!!ord), fill = vals, label=lbls)) +
    scale_fill_gradientn(colours=pall)+
    geom_tile()+
    geom_text(size=5)+
    theme_light()+
    scale_x_discrete(guide=guide_axis(n.dodge=2))+
    theme(legend.position="none",
          axis.text.x = element_text(size=14),
          axis.text.y = element_text(size=14),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  plt_h6 <- ggarrange(plt_h6_C, plt_h6_B, nrow=2)+bgcolor("#ffffff")
  #ggplot2::ggsave("heat_pos_h6.png", plt_h6, path = plots_folder, units = "in", width = 16, height = 3, dpi = 300, bg="white")
}

######### SUMMARIZATION FOR MIXED PROT_MARBOX #########

df_co_occur_all <- read.csv(paste0(data_path, 'all_coocurrence_counts.csv'), stringsAsFactors=FALSE)

df_co_occur_invNT <- df_co_occur_all%>%
  filter(nt=="N")%>%
  summarize(n=sum(cnts))
print(paste0("Invalid NTs: ", df_co_occur_invNT$n))

df_co_occur_invAAs <- df_co_occur_all%>%
  filter(aa=="X")%>%
  summarize(n=sum(cnts))
print(paste0("Invalid AAs: ", df_co_occur_invAAs$n))

df_co_occur_summ <- df_co_occur_all%>%
  summarize(n=sum(cnts))
print(paste0("Global counts: ", df_co_occur_summ$n))

# Removing invalid
df_co_occur <- df_co_occur_all%>%
  filter(nt!="N"&aa!="X")

df_co_occur_summ <- df_co_occur%>%
  summarize(n=sum(cnts))
print(paste0("NEW counts: ", df_co_occur_summ$n))

df_co_occur <- df_co_occur%>%
  complete(aa, nt, mat)%>%
  replace(is.na(.), 0)%>%
  separate(mat, into=c('marb', 'exp_time', 'mat'), sep=' ')%>%
  separate(mat, into=c('helix', 'pos_aa', 'pos_nt'), remove=FALSE)%>%
  mutate(aa=factor(aa, levels=!!aa_order))

marb_order <- c("marRAB", "yba0", "rob", "acnA", "acrAB", "fldB", "fldA", "fpr", "hdeA", "mdtG", "poxB", "purA", "ribA", "slp")
exp_time <- c("INI", "96H")

######### COVARIANCE NT X AA #########

gen_log_ratio <- function(marb_order, exp_time, df_co_occur, aa_order, sigma, save=FALSE) {
  prev_marb <-""
  for (m in marb_order){
    prev_helix<-""
    h3_sz <- sum(grepl("h3", unique(df_co_occur$mat)))
    h6_sz <- sum(grepl("h6", unique(df_co_occur$mat)))
    if (h3_sz!=h6_sz) {
      print("Number of NTs for H3 and H6 must be the same!\nRe-run the python code.")
      break
    }
    
    for (e in unique(df_co_occur$mat)){
      df_exp <- df_co_occur[which(df_co_occur$marb==m & df_co_occur$exp_time==exp_time[2] & df_co_occur$mat==e), ]
      # Calculate relative counts
      df_exp_summ <- df_exp%>%
        group_by(aa)%>%
        summarise(tots=sum(cnts))
      df_exp_sel <- df_exp%>%
        left_join(df_exp_summ, by="aa")%>%
        mutate(rel_cnts=cnts/tots)%>%
        select(aa, nt, rel_cnts)%>%
        spread(nt, rel_cnts)%>%
        column_to_rownames(var="aa")%>%
        select(A, `T`, C, G)
      # To save
      df_exp_ori <- df_exp%>%
        select(aa, nt, cnts)%>%
        spread(nt, cnts)%>%
        left_join(df_exp_summ, by="aa")%>%
        column_to_rownames(var="aa")%>%
        select(A, `T`, C, G, tots)
      
      df_obs <- df_co_occur[which(df_co_occur$marb==m & df_co_occur$exp_time==exp_time[1] & df_co_occur$mat==e), ]
      # Calculate relative counts
      df_obs_summ <- df_obs%>%
        group_by(aa)%>%
        summarise(tots=sum(cnts))
      df_obs_sel <- df_obs%>%
        left_join(df_obs_summ, by="aa")%>%
        mutate(rel_cnts=cnts/tots)%>%
        select(aa, nt, rel_cnts)%>%
        spread(nt, rel_cnts)%>%
        column_to_rownames(var="aa")%>%
        select(A, `T`, C, G)
      
      df_obs_ori <- df_obs%>%
        select(aa, nt, cnts)%>%
        spread(nt, cnts)%>%
        left_join(df_obs_summ, by="aa")%>%
        column_to_rownames(var="aa")%>%
        select(A, `T`, C, G, tots)
      
      df_ratio <- as.data.frame(as.matrix(df_exp_sel)/as.matrix(df_obs_sel))
      df_ratio <- df_ratio%>%
        mutate_if(is.numeric, ~ replace_na(., 0) %>% 
                    replace(., is.infinite(.), 0))%>%
        rownames_to_column(var="aa")%>%
        mutate(aa=factor(aa, levels=!!aa_order))%>%
        arrange(aa)
      
      df_log_ratio <- df_ratio%>%
        left_join(df_exp_summ, by="aa")%>%
        mutate(part1=tots*!!sigma)%>%
        gather(nt, ratio, -aa, -tots, -part1)%>%
        mutate(part2=part1*ratio,
               covars=log(1+part1)-log(1+part2))%>%
        select(-tots, -ratio, -part1, -part2)%>%
        spread(nt, covars)%>%
        select(aa, A, `T`, C, G)
      
      #df_log_ratio1 <- df_ratio%>%
      #  mutate_if(is.numeric, ~ -log(. +1))
      
      df_log_ratio_save <- df_log_ratio
      df_log_ratio$marb <- df_exp$marb[1]
      df_log_ratio$helix <- df_exp$helix[1]
      df_log_ratio$pos_aa <- as.integer(df_exp$pos_aa[1])
      df_log_ratio$pos_nt <- as.integer(df_exp$pos_nt[1])
      
      df_obs_ori$mat <- e
      df_obs_ori$marbox <- m
      df_exp_ori$mat <- e
      df_exp_ori$marbox <- m
      
      df_obs_ori <- df_obs_ori%>%rownames_to_column(var="aa")%>%relocate("aa", .before = 1)
      df_exp_ori <- df_exp_ori%>%rownames_to_column(var="aa")%>%relocate("aa", .before = 1)
      
      if (!exists("df_log_ratio_all")) {
        df_log_ratio_all <- df_log_ratio
        df_cnt_obs <- df_obs_ori
        df_cnt_exp <- df_exp_ori
      } else {
        df_cnt_obs <- rbind(df_cnt_obs, df_obs_ori)
        df_cnt_exp <- rbind(df_cnt_exp, df_exp_ori)
        df_log_ratio_all <- rbind(df_log_ratio_all, df_log_ratio)
      }
      
      if ((prev_marb!=m)&prev_marb!="") {
        if (!exists("df_log_ratio_save_all")) {
          df_log_ratio_save_all <- rbind(df_log_ratio_save_h3, df_log_ratio_save_h6)
          rm(df_log_ratio_save_h3, df_log_ratio_save_h6)
        } else {
          df_log_ratio_save_all <- rbind(df_log_ratio_save_all, rbind(df_log_ratio_save_h3, df_log_ratio_save_h6))
          rm(df_log_ratio_save_h3, df_log_ratio_save_h6)
        }
      }
      
      if (df_exp$helix[1]!=prev_helix) {
        # To export the excel file
        title <- as.data.frame(t(c(paste0(m, " ", df_exp$mat[1]), "A", "T", "C", "G")))
        colnames(title) <- colnames(df_log_ratio_save)
        if (df_exp$helix[1]=="h3") {
          df_log_ratio_save_h3 <- rbind(title, df_log_ratio_save)
        } else {
          df_log_ratio_save_h6 <- rbind(title, df_log_ratio_save)
        }
      } else {
        title <- as.data.frame(t(c(paste0(m, " ", df_exp$mat[1]), "A", "T", "C", "G")))
        colnames(title) <- colnames(df_log_ratio_save)
        df_log_ratio_save <- rbind(title, df_log_ratio_save)
        df_log_ratio_save <- df_log_ratio_save%>%add_column(empty = "", .before = 1)
        if (df_exp$helix[1]=="h3") {
          df_log_ratio_save_h3 <- cbind(df_log_ratio_save_h3, df_log_ratio_save)
        } else {
          df_log_ratio_save_h6 <- cbind(df_log_ratio_save_h6, df_log_ratio_save)
        }
      }
      
      prev_helix <- df_exp$helix[1]
      prev_marb <- m
    }
    
  }
  df_log_ratio_save_all <- rbind(df_log_ratio_save_all, rbind(df_log_ratio_save_h3, df_log_ratio_save_h6))
  #rm(df_log_ratio_save_all, df_cnt_obs, df_cnt_exp, df_log_ratio_all)
  
  if (save) {
    #write.xlsx(df_log_ratio_save_all, file = "covariance_counts.xlsx", sheetName = "Ratios", append = FALSE)
    write.table(df_log_ratio_save_all, paste0(data_path, 'covariance_logs_ratio.csv'), sep=",", row.names = FALSE, col.names = FALSE)
    write.csv(df_cnt_obs, paste0(data_path, 'covariance_counts_ini.csv'), row.names = FALSE)
    write.csv(df_cnt_exp, paste0(data_path, 'covariance_counts_96H.csv'), row.names = FALSE)
  }
  return (df_log_ratio_all)
}

get_ori_nts <- function(df_contacts_AllC, data_path){
  df_marb_pos <- df_contacts_AllC%>%
    select(pos_aa, pos_nt)%>%
    distinct()%>%
    mutate(pos_nt2 = pos_nt-start_pos)
  df_marb_default <- read.csv(paste0(data_path, 'marboxes.txt'), stringsAsFactors=FALSE, sep="\t")
  
  for (n in 1:nrow(df_marb_pos)) {
    col_name <- paste0("pos_", df_marb_pos$pos_aa[n] ,"_",df_marb_pos$pos_nt[n])
    idx <- df_marb_pos$pos_nt2[n]
    df_marb_default[[col_name]] <- substr(df_marb_default[["C"]], idx, idx)
  }
  
  df_marb_default <- df_marb_default%>%
    select(-B, -C)%>%
    gather(pos, nt_ori, -marb)%>%
    separate(pos, into=c(NA, "pos_aa", "pos_nt"))%>%
    mutate(pos_aa=as.integer(pos_aa),
           pos_nt=as.integer(pos_nt))
  return (df_marb_default)
}

get_title_text <- function(df_top) {
  
  df_contacts_pairs <- df_top%>%
    select(pos_aa, pos_nt)%>%distinct()%>%
    group_by(pos_aa)%>%
    mutate(ord=dplyr::row_number())%>%
    spread(pos_aa, pos_nt)%>%
    select(-ord)
  
  df_contacts_pairs <- as.data.frame(t(df_contacts_pairs))
  
  df_contacts_pairs <- df_contacts_pairs%>%
    rownames_to_column(var="pos_aa")%>%
    unite(nt_name, c(V1, V2), sep=" and ")%>%
    mutate(pos_aa=as.integer(pos_aa))
  
  return(df_contacts_pairs)
}

df_log_ratio0 <- gen_log_ratio(marb_order, exp_time, df_co_occur, aa_order, sigma, FALSE)
df_marb_default <- get_ori_nts(df_contacts_AllC, data_path)
df_log_ratio <- df_log_ratio0%>%
  left_join(df_marb_default, by=c("marb","pos_aa", "pos_nt"))

gen_cov_heatmaps <- function(df_contacts_AllC, df_log_ratio, marb_order, cov_folder) {
  
  # Aesthetics
  nt_ord <- c("A_1","T_1","C_1","G_1","A_2","T_2","C_2","G_2")
  pall <- c("#00A933", "#DDE8CB", "#FFFFFF")
  #div_lines <- data.frame(x=c(4.5,4.5), y=c(.5,20.5), 
  #                        xend=c(4.5,4.5), yend=c(20.5,20.5))
  div_lines <- data.frame(x=c(.5,8.5,.5,8.5,4.5,4.5), y=c(20.5,20.5,40.5,40.5,.5,60.5),
                          xend=c(8.5,8.5,8.5,8.5,4.5,4.5), yend=c(20.5,20.5,40.5,40.5,60.5,60.5))

  df_contacts_top2 <- df_contacts_AllC%>%
    filter(order<=2)%>%
    arrange(pos_aa, desc(vals))

  df_contacts_pairs <- get_title_text(df_contacts_top2)

  aa_labs_order <- paste(rev(rep(unique(df_contacts_top2$pos_aa), each=20)), rev(unique(df_log_ratio$aa)), sep="_")

  df_contacts_top2_mats <- df_contacts_top2%>%
    left_join(df_log_ratio, by=c('helix', 'pos_aa', 'pos_nt'), suffix=c("_ori", ""))%>%
    left_join(df_contacts_pairs, by="pos_aa")%>%
    select(-aa_ori)%>%
    gather(nt, log_ratio, -pos_aa, -nt, -nt_ori, -pos_nt, -vals, -aa, -order, -nt_name, -marb, -helix)%>%
    unite(nt_comp, c("nt", "order"), remove = FALSE)%>%
    unite(aa_comp, c("pos_aa", "aa"), remove = FALSE)%>%
    mutate(nt_comp=factor(nt_comp, levels=!!nt_ord),
           lbls=ifelse(log_ratio>=0, "", round(log_ratio, 2)),
           lbls=ifelse(nt_ori==nt, paste0(lbls, "*"), lbls))
  
  for (marb in marb_order) { # 14 marboxes
    for (s2D in unique(df_contacts_top2$helix)) { # 2 helices locations
      
      df_contacts_top2_mat <- df_contacts_top2_mats%>%
        filter(marb==!!marb&helix==!!s2D)
      
      pos_aas <- unique(df_contacts_top2_mat$pos_aa)
      
      #for (id_pos in 1:length(pos_aas)) { # 3 amino acid residues per helix
        
        #df_contacts_top2_mat <- df_contacts_top2_mat%>%
        #  filter(pos_aa==pos_aas[id_pos])
        
      x_axis <- paste0(unique(df_contacts_top2_mat$nt_name), collapse="  |  ")
      y_axis <- paste0(rev(unique(df_contacts_top2_mat$pos_aa)), collapse=" | ")

      plot_covariance <- ggplot(df_contacts_top2_mat, aes(x = nt_comp, y = ordered(aa_comp, levels=!!aa_labs_order), fill = log_ratio)) +
        scale_fill_gradientn(colours=pall)+
        geom_tile()+
        theme_light()+
        scale_x_discrete(position = "top", labels=rep(nt_ori_ord, 2))+
        scale_y_discrete(labels=rep(rev(unique(df_log_ratio$aa)), 3))+
        geom_segment(data=div_lines, aes(x,y,xend=xend, yend=yend), linewidth=.3, inherit.aes=F)+
        theme(legend.position="top",
              axis.text.x = element_text(size=14),
              axis.text.y = element_text(size=14),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
      
      ggplot2::ggsave(paste0(marb, "_", s2D, "_simple.png"), plot_covariance, path = cov_folder, units = "in", width = 16, height = 20, dpi = 300, bg="white")

      plot_covariance <- plot_covariance+
        geom_text(aes(label=lbls), size=5)+
        labs(y=paste0("AA Position ", y_axis), x=paste0("NTs positions: ", x_axis))+
        theme(legend.position="none",
              axis.text.x = element_text(size=14),
              axis.text.y = element_text(size=14),
              axis.title.x = element_text(size=18, face="bold"),
              axis.title.y = element_text(size=18, face="bold", angle=90))
      
      ggplot2::ggsave(paste0(marb, "_", s2D, "_det.png"), plot_covariance, path = cov_folder, units = "in", width = 16, height = 20, dpi = 300, bg="white")
      #}
    }
  }
}

gen_cov_heatmaps(df_contacts_AllC, df_log_ratio, marb_order, cov_folder)

######### NOTATION FOR PYMOL #########

# Automatic annotate positions with pymol
pos_nt1 <- t(df_contacts_top[which(df_contacts_top$order==1), 'pos_nt'])
pos_nt2 <- t(df_contacts_top[which(df_contacts_top$order==2), 'pos_nt'])
pos_nt3 <- t(df_contacts_top[which(df_contacts_top$order==3), 'pos_nt'])
paste0("color blue, resi ", unique(df_contacts_top$pos_aa), 
       "; color red, resi ", pos_nt1, 
       "; color orange, resi ", pos_nt2,
       "; color yellow, resi ", pos_nt3)

paste0("color green, resi ", unique(df_contacts_top$pos_aa), 
       "; color red, resi ", pos_nt1, 
       "; color orange, resi ", pos_nt2,
       "; color yellow, resi ", pos_nt3)


