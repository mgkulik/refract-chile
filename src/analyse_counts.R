library(pacman)
p_load(dplyr, tidyr, gplots, corrplot, ggalluvial, data.table, ggpubr)

source(here::here("data_management.R"))

sigma <- 1/50

target_helix3 <- c(40,43,44)
target_helix6 <- c(90,91,94)

aa_order = c('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V')

SAVE_OUTPUT <- FALSE
REMOVE_SINGLETON <- FALSE

######### TESTING SINGLETON IMPACT #########

if (REMOVE_SINGLETON) {
  
  # TAKE OFF the singletons of the data (can help on the analysis of metrics)
  df_ini <- df_ini%>%filter(cnt>1)
  df_96h <- df_96h%>%filter(cnt>1)
  
  df_ini_sel <- df_ini%>%select(marbox, peptide, cnt)
  df_96h_sel <- df_96h%>%select(marbox, peptide, cnt, order_marbox)
  
  df_intersect_noSingleton <- df_ini_sel%>%
    inner_join(df_96h_sel, by=c("marbox", "peptide"), suffix=c("_ini", "_96h"))%>%
    group_by(marbox)%>%tally()
  
} else {
  # When the singletons are there, this complete session makes sense
  logo_lowcounts <- function(df, min_cnt, tp, exec_name="", sv=FALSE) {
    
    if (exec_name!="") {
      # Separating couts to join later
      df_small <- df%>%select(marbox, peptide, cnt)
      
      if (tp==1) {
        # Gets all values but the one bellow n
        df_1 <- df%>%
          filter(cnt>!!min_cnt)%>%
          select(marbox, peptide, res40, res43, res44, res90, res91, res94)
      } else if (tp==2) {
        # Gets only the values bellow or equal to n (reverse logo)
        df_1 <- df%>%
          filter(cnt<=!!min_cnt)%>%
          select(marbox, peptide, res40, res43, res44, res90, res91, res94)
      } else if (tp==3) {
        # Gets only the top n values (top logos)
        df_1 <- df%>%
          slice_head(n=min_cnt)%>%
          select(marbox, peptide, res40, res43, res44, res90, res91, res94)
      }
      
      df_1 <- df_1%>%
        # Gathering the positions
        gather(pos, aa, -marbox, -peptide)%>%
        # Joining the counts (they'll repeat for each AA)
        left_join(df_small, by=c("marbox", "peptide"))%>%
        arrange(marbox, peptide, pos)%>%
        group_by(marbox, pos, aa)%>%
        # Sum the counts and generate the fractions
        summarize(cnt=sum(cnt))%>%
        mutate(aa=factor(aa, levels=!!aa_order))
      
      # Generate and save the fractions
      df_2 <- df_1%>%
        mutate(cnt=round(cnt/sum(cnt, rm.na=T), 4))%>%
        # Spread the fractions
        spread(pos, cnt)%>%
        group_by(marbox)%>%
        complete(aa)%>%
        replace(is.na(.), 0)
      colnames(df_2) <-  c("marbox", "P0", 0:5)
      
      df_1 <- df_1%>%
        spread(pos, cnt)%>%
        group_by(marbox)%>%
        complete(aa)%>%
        replace(is.na(.), 0)
      colnames(df_1) <-  c("marbox", "P0", 0:5)
    }

    if (exec_name!=""&sv) {
      # Save txt files with all 
      if (tp==2){
        exec_name <- paste0(exec_name, "_rev")
      } else if (tp==3) {
        exec_name <- paste0(exec_name, "_top")
      }
      file_name <- paste0(data_out, exec_name, '_fractions.txt')
      unlink(file_name)
      for (marbox in unique(df_2$marbox)){
        df <- df_2%>%filter(marbox==!!marbox)%>%ungroup()%>%select(-marbox)
        write.table(df,
                    file = file_name,
                    sep = "\t",
                    row.names=F,
                    append=T, quote=F)
      }
      file_name <- paste0(data_out, exec_name, '_counts.txt')
      unlink(file_name)
      for (marbox in unique(df_1$marbox)){
        df <- df_1%>%filter(marbox==!!marbox)%>%ungroup()%>%select(-marbox)
        write.table(df,
                    file = file_name,
                    sep = "\t",
                    row.names=F,
                    append=T, quote=F)
      }
    }
    return(df_1)
  }
  
  get_lowpep <- function(df, min_cnt, tp) {
    # Summarize the peptide/counts
    df_tot <- df%>%
      filter(cnt>!!min_cnt)
    if (tp==1){
      # Sum the counts
      df_tot <- df_tot%>%group_by(marbox)%>%summarize(tot=sum(cnt))
    } else {
      # Count peptides
      df_tot <- df_tot%>%group_by(marbox)%>%tally()
    }
    colnames(df_tot) <- c("marbox", min_cnt)
    return(df_tot)
  }
  
  
#' Main function to generate the counts and fractions in the right format
#' to generate the logos using the automation created for enoLOGOS.
#' @param df source dataframe.
#' @param max_range an integer as maximum value for the loop.
#' @param tp type of filter: 
#'    1) returns the values of all counts above the value of the current loop index;
#'    2) returns only the values bellow or equal the current loop index (reverse logo);
#'    3) returns the top n cases based on the max_range informed (only 1 execution is performed)
#' @param exec_name prefix of the name to be saved to disk
#' @param vals2save a numeric list defining what executions of the loop must 
#' be saved. At list 1 item must be provided when tp=1. Ignored for tp=3 (default=NA).
#' @param max_range an integer as minimum value for the loop (default=0).
#' @return Save specified counts and fractions to disk. For tp=1 also returns
#'         a list with 2 dataframes with the cumulative peptide/counts for the remaining
#'         data when each cout is removed from the set in each step of the loop.
#' @examples
#' gen_logocounts(df1, 1000, 1, "ini", c(0,1,500,1000))  # For all accumulated values when 0,1...1000 is removed
#' gen_logocounts(df1, 1000, 2, "ini", c(0,1,500,1000))  # For all accumulated values for cases with 1,2...1000 counts 
#' gen_logocounts(df1, 1000, 3, "ini")                   # Filters the top 1000 rows
gen_logocounts <- function(df, max_range, tp, exec_name, vals2save=NA, min_range=0, sv=FALSE) {
  
  lst_dfs <- list()
  if (tp==2) {
    # Avoid an empty counts df when only the removed are returned
    if (min_range==0) {
      min_range = 1
    }
  } else if (tp==3){
    # No loop for the top and the files are always saved to disk
    min_range <- max_range
    vals2save <- c(max_range)
  }

  for (n in min_range:max_range){
    exec_name_full = ""
    if ((n%in%vals2save)) {
      exec_name_full = paste0(exec_name, n)
      temp <- logo_lowcounts(df, n, tp, exec_name_full, sv)
    }
    if (tp==1) {
      df_cnt <- get_lowpep(df, n, 1)
      df_pep <- get_lowpep(df, n, 2)
      if (n==min_range) {
        df_cnt_diff <- df_cnt
        df_pep_diff <- df_pep
      } else {
        df_cnt_diff <- df_cnt_diff%>%left_join(df_cnt, by=c("marbox"))
        df_pep_diff <- df_pep_diff%>%left_join(df_pep, by=c("marbox"))
      }
    }
  }
  if (tp==1) {
    lst_dfs[[1]] <- df_cnt_diff
    lst_dfs[[2]] <- df_pep_diff
  }
  return(lst_dfs)
}

  # Removing the singletons
  x <- gen_logocounts(df_ini, 1, 1, "ini", c(0,1), sv=TRUE)
  x <- gen_logocounts(df_96h, 1, 1, "96h", c(0,1), sv=TRUE)
  rm(x)
  
  # Just the singletons to compare
  gen_logocounts(df_ini, 1, 2, "ini", c(0, 1))
  gen_logocounts(df_96h, 1, 2, "96h", c(0, 1))
  
  # Removing all up to 1000
  lst_diff_ini <- gen_logocounts(df_ini, 1000, 1, "ini", c(1000), sv=TRUE)
  lst_diff_96h <- gen_logocounts(df_96h, 1000, 1, "96h", c(1000), sv=TRUE)
  
  df_ini_diff_cnt <- lst_diff_ini[[1]]
  df_ini_diff_pep <- lst_diff_ini[[2]]
  
  # Taking just the cases up to 1000 counts to compare
  gen_logocounts(df_ini, 1000, 2, "ini", c(1000), sv=TRUE)
  gen_logocounts(df_96h, 1000, 2, "96h", c(1000), sv=TRUE)
  
  # Taking the top 10 to get the fractions
  gen_logocounts(df_ini, 10, 3, "ini")
  gen_logocounts(df_96h, 10, 3, "96h")
  
  # Taking the first 1000 to get the fractions
  gen_logocounts(df_ini, 1000, 3, "ini", sv=TRUE)
  gen_logocounts(df_96h, 1000, 3, "96h", sv=TRUE)
  
  
  # Tried to check the impact of removing the cases with 50 less counts, doesn't worth it
  get_diffNOlow <- function(df_diff) {
    df_diff1 <- transpose(df_diff)
    colnames(df_diff1) <- df_diff1[1,]
    df_diff1 <- df_diff1[2:nrow(df_diff1), ]
    rownames(df_diff1) <- 1:nrow(df_diff1)
    df_diff1$min <- 0:(nrow(df_diff1)-1)
    
    df_diff1 <- df_diff1%>%
      mutate_if(is.character, as.integer)%>%
      gather(marbox, cnt, -min)%>%
      group_by(marbox)%>%
      mutate(diff=abs(cnt-lag(cnt, default=first(cnt))))
    
    df_diff_sel <- df_diff%>%select(marbox, `0`)
    
    df_diff1_summ <- df_diff1%>%
      group_by(marbox)%>%
      summarize(diff=sum(diff))%>%
      left_join(df_diff_sel, by="marbox")%>%
      mutate(tot_diff=`0`-diff,
             frac_diff=tot_diff/`0`,
             marbox=factor(marbox, levels=!!marbox_order))%>%
      arrange(marbox)
    return(df_diff1_summ)
  }
  
  df_ini_low_summ <- get_diffNOlow(df_ini_cnt_diff)
  df_96h_low_summ <- get_diffNOlow(df_96h_cnt_diff)
}

######### NORMALIZATION STRATEGIES #########
df_merged <- df_ini_sel%>%
  full_join(df_96h_sel, by=c("marbox", "peptide"), suffix=c("_ini", "_96h"))%>%
  mutate_at(vars(c(3,4)), ~replace(., is.na(.), 0))%>%
  filter(cnt_ini>5&cnt_96h>5)%>%
  arrange(marbox, desc(cnt_96h))%>%
  mutate(cum_sum=cumsum(cnt_96h),
         diff_96h=abs(cnt_96h-lag(cnt_96h)),
         log_diff_96h=log(1+diff_96h),
         frac_ini=cnt_ini/sum(cnt_ini),
         frac_96h=cnt_96h/sum(cnt_96h),
         ratio_frac=replace_na(frac_96h/frac_ini, 0),
         ratio_frac=ifelse(is.infinite(ratio_frac), 0, ratio_frac),
         log_ratio=log(1+ratio_frac))%>%
  arrange(marbox, desc(cnt_96h))

if (SAVE_OUTPUT) {
  metric_name <- "log_ratio"
  marb <- "marRAB"
  #marb <- "fpr"
  
  metric_marb1 <- df_merged%>%filter(marbox==!!marb)%>%
    ungroup()%>%
    arrange(desc(cnt_96h))%>%
    mutate(new_order=row_number())%>%
    #slice(2:nrow(metric_marb))%>%
    slice_head(n=1000)%>%
    filter(cnt_96h>1)%>%
    select(new_order, peptide, cnt_ini, frac_ini, cnt_96h, frac_96h, ratio_frac, !!metric_name)
  
  plot_metric <- ggplot(metric_marb1, aes(x=new_order, y=!!sym(metric_name)))+
    geom_line()+
    theme(axis.text.x = element_text(size=12), 
          axis.text.y = element_text(size=12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=14, face="bold"))
  
  metric_marb1 <- df_merged%>%filter(marbox==!!marb)%>%
    ungroup()%>%
    arrange(desc(cnt_96h))%>%
    mutate(new_order=row_number())%>%
    #slice(2:nrow(metric_marb))%>%
    slice_head(n=1000)%>%
    #filter(!!sym(metric_name)<1)%>%
    select(new_order, !!metric_name)
  
  plot_metric <- ggplot(metric_marb1, aes(x=new_order, y=!!sym(metric_name)))+
    geom_line()+
    theme(axis.text.x = element_text(size=12), 
          axis.text.y = element_text(size=12),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=14, face="bold"))
  
  ggplot2::ggsave(paste0("metrics/metric_", metric_name, "_", marb,".png"), plot_metric, units = "in", width = 3, height = 2, dpi = 300, bg="white")
}

gen_cumsum <- function(df_merged, n=NA) {
  
  lst_plt <- list()
  i<-1
  
  for (marb in marbox_order) {
    
    merge_marb <- df_merged%>%filter(marbox==!!marb)%>%
      ungroup()
    
    if (is.na(n)) {
      n_top <- nrow(merge_marb)
    } else {
      n_top <- n
    }
    
    merge_marb1 <- merge_marb%>%
      arrange(desc(cnt_96h))%>%
      mutate(new_order=row_number(),
             cumsum=cumsum(cnt_96h))%>%
      slice(1:!!n_top)%>%
      select(new_order,cumsum)%>%
      mutate(tp="Descent counts")
    
    merge_marb2 <- merge_marb%>%
      arrange(desc(diff_exp))%>%
      mutate(new_order=row_number(),
             cumsum=cumsum(cnt_96h))%>%
      slice(1:!!n_top)%>%
      select(new_order,cumsum)%>%
      mutate(tp="Difference of expression")
    
    merge_marb3 <- merge_marb%>%
      arrange(desc(diff_frac))%>%
      mutate(new_order=row_number(),
             cumsum=cumsum(cnt_96h))%>%
      slice(1:!!n_top)%>%
      select(new_order,cumsum)%>%
      mutate(tp="Difference in fractions")
    
    merge_marb4 <- merge_marb%>%
      arrange(desc(ratio_frac))%>%
      mutate(new_order=row_number(),
             cumsum=cumsum(cnt_96h))%>%
      slice(1:!!n_top)%>%
      select(new_order, cumsum)%>%
      mutate(tp="Ratio of fractions")
    
    merge_marb5 <- merge_marb%>%
      arrange(desc(diff_log))%>%
      mutate(new_order=row_number(),
             cumsum=cumsum(cnt_96h))%>%
      slice(1:!!n_top)%>%
      select(new_order, cumsum)%>%
      mutate(tp="Difference of logs")
    
    merge_marb6 <- merge_marb%>%
      arrange(desc(ratio_log))%>%
      mutate(new_order=row_number(),
             cumsum=cumsum(cnt_96h))%>%
      slice(1:!!n_top)%>%
      select(new_order, cumsum)%>%
      mutate(tp="Ratio of logs")
    
    merge_marb_all <- rbind(merge_marb1, merge_marb2)
    merge_marb_all <- rbind(merge_marb_all, merge_marb3)
    merge_marb_all <- rbind(merge_marb_all, merge_marb4)
    merge_marb_all <- rbind(merge_marb_all, merge_marb5)
    merge_marb_all <- rbind(merge_marb_all, merge_marb6)
    
    order_metrics <- c("Descent counts", "Difference of expression", 
                       "Difference in fractions", "Ratio of fractions",
                       "Difference of logs", "Ratio of logs")
    merge_marb_all <- merge_marb_all%>%
      mutate(tp=factor(tp, levels=!!order_metrics))
    
    
    plot_cumsum <- ggplot(merge_marb_all, aes(x=new_order, y=cumsum, color=tp))+
      geom_line()+
      scale_fill_brewer(palette="Set2")+
      labs(x=marb)+
      guides(fill=guide_legend(ncol=3, nrow=2, byrow=TRUE))+
      theme(axis.text.x = element_text(size=14),
            axis.title.x = element_text(size=18, face="bold"),
            legend.text = element_text(size=18),
            legend.title = element_text(size=20),
            legend.spacing.x = unit(1, 'cm'))
    
    lst_plt[[i]] <- plot_cumsum
    i=i+1
  }
  
  marb_cumsum <- ggarrange(plotlist=lst_plt, ncol=5, nrow=3, common.legend = TRUE, legend="top", labels=NULL)+bgcolor("#ffffff")
  return (marb_cumsum)
}

plt_all <- gen_cumsum(df_merged, 1000)
ggplot2::ggsave(paste0("metrics/metrics_noSingletons_1k.png"), plt_all, units = "in", width = 26, height = 12, dpi = 300, bg="white")

######### TESTING AA COVARIANCE (Later implemented in python) #########
cols_h3 <- paste0("res", target_helix3)
all_combs_h3 <- as.data.frame(gtools::permutations(n = length(cols_h3), r = 2, v = cols_h3, repeats.allowed = FALSE))

filter1 <- "W"

#i=1

for (i in 1:nrow(all_combs_h3)){

  if ((i==4|i==6)&filter1!=""){
    df_obs_AA <- df_ini%>%filter(!!sym(all_combs_h3[1,1])==!!filter1)
    df_exp_AA <- df_96h%>%filter(!!sym(all_combs_h3[1,1])==!!filter1)
  } else {
    df_obs_AA <- df_ini
    df_exp_AA <- df_96h
  }
  
  # INI
  df_obs_AA <- tibble(df_obs_AA%>%
    select(marbox, !!all_combs_h3[i,1], !!all_combs_h3[i,2], cnt))%>%
    complete(marbox, !!sym(all_combs_h3[i,1]), !!sym(all_combs_h3[i,2]))%>%
    group_by(marbox, !!sym(all_combs_h3[i,1]), !!sym(all_combs_h3[i,2]))%>%
    summarise(cnt=sum(cnt))%>%
    replace(is.na(.), 0)%>%
    mutate(aa1=factor(!!sym(all_combs_h3[i,1]), levels=!!aa_order))%>%
    ungroup()%>%
    select(-!!all_combs_h3[i,1])%>%
    arrange(marbox, aa1)
  
  # 96H
  df_exp_AA <- tibble(df_96h%>%
                        select(marbox, !!all_combs_h3[i,1], !!all_combs_h3[i,2], cnt))%>%
    complete(marbox, !!sym(all_combs_h3[i,1]), !!sym(all_combs_h3[i,2]))%>%
    group_by(marbox, !!sym(all_combs_h3[i,1]), !!sym(all_combs_h3[i,2]))%>%
    summarise(cnt=sum(cnt))%>%
    replace(is.na(.), 0)%>%
    mutate(aa1=factor(!!sym(all_combs_h3[i,1]), levels=!!aa_order))%>%
    ungroup()%>%
    select(-!!all_combs_h3[i,1])%>%
    arrange(marbox, aa1)
  
  #marb<-"yba0"
  
  df_obs_AA_marb <- df_obs_AA%>%
    ungroup()%>%filter(marbox==!!marb)%>%
    select(-marbox)
  
  df_obs_AA_marb_summ <- df_obs_AA_marb%>%
    group_by(aa1)%>%
    summarise(tots=sum(cnt))
  
  #df_obs_AA2_marb_summ <- df_obs_AA_marb%>%
  #  group_by(!!sym(all_combs_h3[i,2]))%>%
  #  summarise(tots=sum(cnt))
    
  df_obs_AA_marb_rel <- df_obs_AA_marb%>%
    left_join(df_obs_AA_marb_summ, by="aa1")%>%
    mutate(rel_cnt=cnt/tots, .keep="unused")%>%
    pivot_wider(names_from = !!sym(all_combs_h3[i,2]), values_from = rel_cnt)%>%
    tibble::column_to_rownames(var="aa1")%>%
    select(!!aa_order)
  
  df_obs_AA_marb <- df_obs_AA_marb%>%
    pivot_wider(names_from = !!sym(all_combs_h3[i,2]), values_from = cnt)%>%
    tibble::column_to_rownames(var="aa1")%>%
    select(!!aa_order)
  
  write.csv(df_obs_AA_marb, paste0(files_path, paste0(marb, '_test_ini_', paste(all_combs_h3[i,1], all_combs_h3[i,2], sep="_"), '.csv')), row.names = FALSE)
  
  df_exp_AA_marb <- df_exp_AA%>%
    ungroup()%>%filter(marbox==!!marb)%>%
    select(-marbox)
  
  df_exp_AA_marb_summ <- df_exp_AA_marb%>%
    group_by(aa1)%>%
    summarise(tots=sum(cnt))
  
  df_exp_AA_marb_rel <- df_exp_AA_marb%>%
    left_join(df_exp_AA_marb_summ, by="aa1")%>%
    mutate(rel_cnt=cnt/tots, .keep="unused")%>%
    pivot_wider(names_from = !!sym(all_combs_h3[i,2]), values_from = rel_cnt)%>%
    tibble::column_to_rownames(var="aa1")%>%
    select(!!aa_order)
  
  df_exp_AA_marb <- df_exp_AA_marb%>%
    pivot_wider(names_from = !!sym(all_combs_h3[i,2]), values_from = cnt)%>%
    tibble::column_to_rownames(var="aa1")%>%
    select(!!aa_order)
  
  write.csv(df_exp_AA_marb, paste0(files_path, paste0(marb, '_test_96h_', paste(all_combs_h3[i,1], all_combs_h3[i,2], sep="_"), '.csv')), row.names = FALSE)
  
  rank_red <- df_exp_AA_marb-df_obs_AA_marb
  rank_red <- rank_red%>%
    tibble::rownames_to_column(var="aa1")%>%
    gather("aa2","value", -aa1)
  
  df_ratio <- as.data.frame(as.matrix(df_exp_AA_marb_rel)/as.matrix(df_obs_AA_marb_rel))
  df_ratio <- df_ratio%>%
    mutate_if(is.numeric, ~ replace_na(., 0) %>% 
                replace(., is.infinite(.), 0))%>%
    tibble::rownames_to_column(var="aa1")
  
  df_log_ratio <- df_ratio%>%
    left_join(df_exp_AA_marb_summ, by="aa1")%>%
    mutate(part1=tots*!!sigma)%>%
    gather(aa2, ratio, -aa1, -tots, -part1)%>%
    mutate(part2=part1*ratio,
           covars=log(1+part1)-log(1+part2))%>%
    select(-tots, -ratio, -part1, -part2)%>%
    pivot_wider(names_from = aa2, values_from = covars)
  colnames(df_log_ratio)[1] <- paste(all_combs_h3[i,1], all_combs_h3[i,2], sep="_")
  
  write.csv(df_log_ratio, paste0(files_path, paste0(marb, '_test_covariance_', paste(all_combs_h3[i,1], all_combs_h3[i,2], sep="_"), '.csv')), row.names = FALSE)

}
