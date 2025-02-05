# Current Version: 1.3.0
# Version history
# ...
# 1.1.1 - added conditions to bypass percentile table generation if tumor or stroma segments are absent
# 1.1.2 - changed 'Exported dataset' to 1 when reading Excel sheets
# 1.1.3 - sanitize datadir paths via wrapper
# 1.1.5 - fixed outlier detection to include all antibodies instead of just pos controls
# 1.1.6 - expand antibody plots to include failed antibodies from QC step
# 1.2.0 - added new input file of control genes as RUV-III feature selection
#       - applied RUV-III to TMA data for outlier detection and plotting
#       - created variable k for changing default PCs used for RUV-III
# 1.2.1 -
# 1.2.2 - create output for melt.tma for pre-RUV data for Westgard QC
# 1.3.0 - introduce QC check input functions to screen data entry errors

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gtable))
suppressPackageStartupMessages(library(tidyr))

source('dsp_inputs.R')
source('evaluation.R')
source('helpers.R')
source('normalization.R')
source('plots.R')
source('summarization.R')
source('qc_check.R')

options(datatable.rbindlist.check="warning")
options(datatable.optimize=1)
options(error=traceback)

### START ARGS
args <- commandArgs(trailingOnly=TRUE)
my_samp <- args[1]
runid <- args[2]
coh <- args[3]

# Metadata filepaths
datadir <- args[4]
ab_info <- args[5]
low_probes <- args[6]
control_type <- args[7]
dsp_meta <- args[8]

# Output options
exp.out <- args[9]
seg.proc.out <- args[10]
melt.tma.out.pre <- args[11]
melt.tma.out.post <- args[12]
report.out <- args[13]
excel.out.tum <- args[14]
excel.out.str <- args[15]

# positive cntrl file
pos_cntrls <- args[16]
cntrl_genes <- args[17]

# Set constants
exp.regex <- "[0-9]{8}-[0-9]{2}"
#good_tma <- c("12212022-01", "01122023-01", "01192023-01", "01202023-01", "01252023-01", "01262023-01", "02032023-01", "02082023-01", "02102023-01", "02152023-01", "02282023-01", "03012023-01", "03082023-01", "03152023-01", "04052023-01")

# Load metadata
paths <- data.table(read.xlsx(ab_info, sheet="parsed"))
exp.low <- data.table(openxlsx::read.xlsx(low_probes, colNames=F))
control.type <- data.table(openxlsx::read.xlsx(control_type, startRow=2))
#stopifnot(control.type[,.N,by=.(lower_secondary, name, type)][,all(N==1)])
if(!(control.type[,.N,by=.(lower_secondary, name, type)][,all(N==1)])){
  message("Error: check control type input file")
  quit(status=5)
}
dsp.meta <- data.table(read.xlsx(dsp_meta))
good_tma <- unique(as.vector(dsp.meta %>% filter(use_tma=='y') %>% select(Date_run))[[1]])
pos.cntrls <- data.table(read.csv(pos_cntrls, sep = '\t', header=F))

# Check to see if my_samp is in dsp.meta, stop if it's not.
#stopifnot(nrow(dsp.meta[`Specimen.ID` %in% my_samp]) > 0)
if(nrow(dsp.meta[`Specimen.ID` %in% my_samp]) == 0){
  message("Sample ID invalid.")
  quit(status=5)
}

batch.dt <- data.table(files=list.files(datadir, pattern="[-0-9A-Za-z]_[0-9]{8}-[0-9]{2}.xlsx", recursive = F, full.names=T))
batch.dt[,batch:=str_extract(files, exp.regex)]

res.list <- process_batches(batch.dt, sheet=1)

qc.meta <- res.list$meta

#clean up sample labeling if necessary
qc.meta[`Segment (Name/ Label)` == "Segment 1",.(.N, max_roi=length(unique(croi))),by=.(sample_id, batch)][(N != 3) | (N != max_roi)]
qc.meta[sample_id=="01142", `:=`(sample_id="2001142")]
#Pull out the QC corrected data

# Perform data entry QC
input_qc()
qc.meta <- remove_dummy()

all.abund <- Reduce(function(x,y){

  merge(x,y, by=c("ProbeName"))

}, lapply(res.list$data, function(x) x$qc[,-c(1:3)]))
#stopifnot(nrow(all.abund) == res.list$data[[1]]$qc[,.N])
#stopifnot((ncol(all.abund)-1) == sum(sapply(res.list$data, function(x) ncol(x$qc)-4 )))
if((nrow(all.abund) != res.list$data[[1]]$qc[,.N])|| (ncol(all.abund)-1) != sum(sapply(res.list$data, function(x) ncol(x$qc)-4 ))){
    message("Dimension error with processed batch data")
}
abund.mat <- log2(as.matrix(all.abund[,-1,with=F]))
rownames(abund.mat) <- all.abund[,ProbeName]

# Get TMA submatrix
tma.meta <- qc.meta[`Segment (Name/ Label)`=="Full ROI" | `Segment (Name/ Label)`=="Geometric Segment"]
tma.meta[,lower_sample:=tolower(sample_id)]

#stopifnot(]control.type[,.N,by=.(lower_secondary, name, type)[,all(N==1)])
if(!(control.type[,.N,by=.(lower_secondary, name, type)][,all(N==1)])){
    message("Ambiguous and/or duplicate alternate TMA names -- check control type input file.")
    quit(status=5)
}

tma.meta <- merge(tma.meta, control.type[,.(lower_sample=lower_secondary, name, type)], by="lower_sample", all=F)
tma.meta <- tma.meta[`batch` %in% good_tma | `batch` %in% runid]
#stopifnot(tma.meta[,.N,by=batch][,all(N==19)])

tma.abund <- abund.mat[,tma.meta$barcode]

## QC of experimental samples with respect to ROI
exp.meta <- qc.meta[`sample_id` %in% dsp.meta$Specimen.ID]
exp.meta$`Segment (Name/ Label)` <- ifelse(exp.meta$`Segment (Name/ Label)` == 'Geometric Segment', "Segment 1", exp.meta$`Segment (Name/ Label)`)
exp.meta$`Segment (Name/ Label)` <- ifelse(exp.meta$`Segment (Name/ Label)` == 'Full ROI', "Segment 1", exp.meta$`Segment (Name/ Label)`)

exp.abund <- abund.mat[, exp.meta$barcode]
#note we do this with the raw data as otherwise will constantly have to monitor due to normalization fluxes
exp.meta <- flag_roi(exp.abund, exp.meta, roi.thresh=.85)
#note the handful that didn't pass this QC, note that 432 would likely be considered marginal and kept
exp.meta[rm_croi == T,.(`Segment (Name/ Label)`, batch, sample_id, croi, max_cor, rm_croi)]
#for now, removing all
exp.meta <- exp.meta[rm_croi == F]
exp.abund <- abund.mat[, exp.meta$barcode]
#at this point save the metadata and raw abundance
save(exp.meta, exp.abund, tma.abund, tma.meta, file=paste0(exp.out))

### Normalization and summarization of cohort
# Determine relevant samples / batches and compute normalization factors
relevant.meta <- exp.meta[sample_id %in% dsp.meta$Specimen.ID]
# Figure out how we want to number batches.
relevant.meta[,num_batch:=batch]
relevant.abund <- exp.abund[,relevant.meta$barcode]

tma.abund <- abund.mat[,tma.meta$barcode]
tma.meta <- tma.meta[batch %in% relevant.meta$batch]
# Figure out how we want to number batches.
tma.meta[,num_batch:=batch]

# PRE-RUV melt tma for westgard rules
melt.tma <- data.table(reshape2::melt(tma.abund, as.is=T))
names(melt.tma) <- c("ProbeName", "barcode", "abundance")
melt.tma <- merge(melt.tma, tma.meta[,.(barcode, name, batch)], by="barcode")
melt.tma <- merge(paths[,.(ProbeName, igg)], melt.tma, by="ProbeName", all=T)
melt.tma[,fac_batch:=factor(batch)]
#add in values for corresponding igg
melt.tma[ProbeName %in% c("Ms IgG1",  "Ms IgG2a", "Rb IgG"), igg:=ProbeName]
melt.tma <- merge(melt.tma, melt.tma[ProbeName %in% c("Ms IgG1",  "Ms IgG2a", "Rb IgG"),.(igg=ProbeName, barcode, igg_abund=abundance)], by=c("igg", "barcode"), all.x=T, all.y=F)
melt.tma[,perc_igg:=(abundance/igg_abund)*100]
melt.tma[,ceil_igg:=pmin(perc_igg, 100)]
# Add columns splitting batch date into month-day and year for sorting
melt.tma$monthday<- substr(melt.tma$batch, 1, 4)
melt.tma$monthday<- str_remove(melt.tma$monthday, "^0+")
melt.tma$year<- substr(melt.tma$batch, 5, 8)
# Create combined name_ProbeName column
melt.tma$name_ProbeName<- str_c(melt.tma$name,'_',melt.tma$ProbeName)
write.csv(melt.tma, file=paste0(melt.tma.out.pre), row.names=F)

# Enter custom set of control genes for RUV feature selection
#cntrl.abs <- setdiff(rownames(tma.abund), exp.low[[1]])
cntrl.abs <- read.csv(cntrl_genes, header=F)[[1]]
use_k <- 4

# Remove targets in list
rm_ab <- c('p-4EBP1')
tma.abund <- tma.abund[!(row.names(tma.abund) %in% rm_ab),]
relevant.abund <- relevant.abund[!(row.names(relevant.abund) %in% rm_ab),]

segment.proc <- preprocess_dsp_tma(tma.meta, tma.abund, relevant.meta, relevant.abund, igg.map=paths, bg.method=c('none'), controls=cntrl.abs, use.type='quant', k=use_k, num.roi.avg=3)
# Re-arrange segment.proc output
tma.abund <- segment.proc[[3]]
segment.proc <- segment.proc[1:2]

# Use the summarized metadata to deal with replicates and form groups
#create per segment summarized metadata
avg.meta <- rbindlist(lapply(segment.proc, function(x){
  x$meta
}), idcol="Segment (Name/ Label)")
## Batch 19 samples lost here
my.meta <- merge(dsp.meta[,.(num_batch=Date_run, sample_id=Specimen.ID, Specimen.ID, cohort, code)],
                 avg.meta[,.(`Segment (Name/ Label)`, num_batch, sample_id, avg_barcode)], by=c("num_batch", "sample_id"))
#minus two from above

#choose a technical replicate
## First determine how correlated they are to each other in segment abundance
sapply(names(segment.proc), function(x){

  tmp.meta <- my.meta[`Segment (Name/ Label)`==x]
  tmp.abund <- segment.proc[[x]]$avg_abund[,tmp.meta$avg_barcode]
  pair.cors <- sapply(split(tmp.meta, by=c("Specimen.ID", "code")), function(y){
    min(cor(tmp.abund[,y$avg_barcode]))
  })
  pair.cors
})

my.meta <- my.meta[!duplicated(cbind(`Segment (Name/ Label)`, Specimen.ID, num_batch, code)),]
#Here define reference vs experimental
my.meta[cohort==coh,Best_Response:="Ref"]
save(my.meta, segment.proc, file=paste0(seg.proc.out))

# Form antibody scores as the quantiles relevant to reference cohort
ref_samps <- my.meta[Best_Response == "Ref",unique(sample_id)]
ref_samps <- ref_samps[!ref_samps %in% my_samp]
# If there is no ref data for segment 3 (sarcomas) then we'll just look at segment 1.
if (all(unique(my.meta[`sample_id` %in% ref_samps]$`Segment (Name/ Label)`) == "Segment 1")) {
  quant.list <- score_abs(segment.proc$`Segment 1`, ref.samples=ref_samps, stroma=F,score.type="quant")
  pat.quants <- quant.list$scores[,"Segment (Name/ Label)":="Segment 1"]
}else{
  quant.list <- score_abs(segment.proc, ref.samples=ref_samps, stroma=T,score.type="quant")
  #combine
  pat.quants <- rbindlist(lapply(quant.list, "[[", "scores"), idcol="Segment (Name/ Label)")
}
my.scores <- merge(my.meta, pat.quants, by=c("Segment (Name/ Label)", "avg_barcode"))
# JHL: In the case of identical sample id's, get rid of the one that we are not currently analyzing.
my.scores <- my.scores[!(`sample_id` == my_samp & `num_batch` != runid)]

#getting pathways in order
#use.paths <- paths[analysis_pathway %in% c("Expression Controls", "N/A")==F]
#use.paths[analysis_pathway == "Tumor Markers", analysis_pathway:="Other Markers"]
use.paths <- paths
#path.ord <- c("Cell Cycle", "PI3K/AKT pathway", "RAS/MAPK pathway", "Tumor Markers", "Cell Death", "Immune Markers")
path.ord <- unique(paths$analysis_pathway)
use.paths[,`:=`(path_ord=factor(analysis_pathway, levels=path.ord, ordered=T),
                ab_ord=factor(ab, levels=ab, ordered=T))]

my.scores <- merge(use.paths[,.(ab_ord, ProbeName, path_ord)], my.scores, by="ProbeName", all=F)
my.scores[,comb_id:=Specimen.ID]
my.scores[,segment_label:=ifelse(`Segment (Name/ Label)` == "Segment 1", "tumor", "stroma")]
my.scores[,sample_ord:=sample_id]

pt.ord <- my.scores[,.N,by=.(comb_id, Best_Response)][order(Best_Response)]
my.scores[,patient_ord:=factor(comb_id, levels=pt.ord$comb_id, ordered=T)]
#output abundance for reference samples
if (all(unique(my.meta[`sample_id` %in% ref_samps]$`Segment (Name/ Label)`) == "Segment 1")) {
  ref.abund <- quant.list$ref_abund[,"Segment (Name/ Label)":="Segment 1"]
}else{
  ref.abund <- rbindlist(lapply(quant.list, "[[", "ref_abund"), idcol="Segment (Name/ Label)")
}
ref.abund[,segment_label:=ifelse(`Segment (Name/ Label)` == "Segment 1", "tumor", "stroma")]
ref.abund <- merge(use.paths[,.(ab_ord, ProbeName, path_ord)], ref.abund, by="ProbeName", all=F)

# Steps used to provide data for antibody plots.
melt.tma <- data.table(reshape2::melt(tma.abund, as.is=T))
names(melt.tma) <- c("ProbeName", "barcode", "abundance")
melt.tma <- merge(melt.tma, tma.meta[,.(barcode, name, batch)], by="barcode")
melt.tma <- merge(paths[,.(ProbeName, igg)], melt.tma, by="ProbeName", all=T)
melt.tma[,fac_batch:=factor(batch)]
#add in values for corresponding igg
melt.tma[ProbeName %in% c("Ms IgG1",  "Ms IgG2a", "Rb IgG"), igg:=ProbeName]
melt.tma <- merge(melt.tma, melt.tma[ProbeName %in% c("Ms IgG1",  "Ms IgG2a", "Rb IgG"),.(igg=ProbeName, barcode, igg_abund=abundance)], by=c("igg", "barcode"), all.x=T, all.y=F)
melt.tma[,perc_igg:=(abundance/igg_abund)*100]
melt.tma[,ceil_igg:=pmin(perc_igg, 100)]

# Add columns splitting batch date into month-day and year for sorting
melt.tma$monthday<- substr(melt.tma$batch, 1, 4)
melt.tma$monthday<- str_remove(melt.tma$monthday, "^0+")
melt.tma$year<- substr(melt.tma$batch, 5, 8)

# Create combined name_ProbeName column
melt.tma$name_ProbeName<- str_c(melt.tma$name,'_',melt.tma$ProbeName)
# Write out csv for Westgard rules script in Galaxy wf
write.csv(melt.tma, file=paste0(melt.tma.out.post), row.names=F)

ref.batches <- melt.tma[batch != runid]
cur.batch <- melt.tma[batch == runid]
# Get the number of runs from the metadata sheet
run_no <- length(unique(dsp.meta$Date_run))
use.pal <- scales::hue_pal()(run_no)

# Remove unused factors
samp.scores <- my.scores[`sample_id` == my_samp]
samp.scores$patient_ord <- droplevels(samp.scores$patient_ord, except=my_samp)

### WRITE TO PDF

# Check if current sample is in score table
if(!(nrow(samp.scores)>0)){
  message("No data found for current sample ID.")
  quit(status=5)
}

# Overall Plotting
plot.list <- loli_plot(score.dt=samp.scores, ref.dt=ref.abund, coh)
pdf(file=paste0(report.out), width=16, height=16)

# Restrict ab plots and outlier detection to those ab/cell line combos included in pos.cntrls
colnames(pos.cntrls) <- c('ProbeName','name')
ab.batches <- ref.batches %>% inner_join(pos.cntrls, by=c('ProbeName','name'))
ab.batches.cur <- cur.batch %>% inner_join(pos.cntrls, by=c('ProbeName','name'))
clia_abs <- unique(ab.batches$ProbeName)

# Outlier detection -- Zscore method
datalist = list()
k <- 1

for (i in seq(1, length(clia_abs))){
  for (j in seq(1, length(unique(ab.batches.cur$name)))){
    ref <- ab.batches %>% filter(ProbeName==clia_abs[i] & name == unique(ab.batches.cur$name)[j]) %>%
      select(batch, ProbeName, name, abundance) %>%
      group_by(name)
    mu <- mean(ref$abundance)
    sigma <- sd(ref$abundance)

    cur <- ab.batches.cur %>% filter(ProbeName==clia_abs[i] & name == unique(ab.batches.cur$name)[j]) %>%
      select(batch, ProbeName, name, abundance) %>%
      group_by(name) %>%
      mutate(mean = mu,
             sd = sigma,
             zscore = (abundance - mu)/sigma,
             outlier = case_when(abs(zscore) >= 3 ~ 'rare'))
    cur <- na.omit(cur)
    datalist[[k]] <- cur
    k <- k+1
  }
}

outlier_df = do.call(rbind, datalist)

# Get failed antibodies (Ab combos with 5+ flagged outliers)
failed_ab <- as.data.frame(table(outlier_df$ProbeName))
if (nrow(failed_ab)>0){
  failed_ab <- failed_ab %>% filter(Freq >= 5) %>% select(Var1)
  colnames(failed_ab) <- c('Failed Antibodies')
}

# First, produce Cover Sheet
tt <- ttheme_default(base_size = 16)
tt1 <- ttheme_minimal(core=list(fg_params=list(fontface=3, fontsize=12)))
tt_cover <- ttheme_minimal(core=list(bg_params = list(fill = blues9[1:4], col=NA),
                                     fg_params=list(fontface=3, fontsize=23)),
                           colhead=list(fg_params=list(col="darkblue", fontface=4L, fontsize=30)))

summ_df <- as.data.frame(c(my_samp, runid, as.character(Sys.Date())), header=FALSE)
rownames(summ_df) <- c('SAMPLE ID: ', 'RUN ID: ', 'RUN DATE: ')
colnames(summ_df) <- c('Nanostring_DSP')

reference_batches <- ref.batches %>% arrange(year, monthday) %>% select(batch) %>% distinct(batch)
colnames(reference_batches) <- c('Reference Batches')

gc1 <- tableGrob(summ_df, theme=tt_cover)

if (nrow(reference_batches)>20){
  rb1 <- tableGrob(reference_batches[1:20,1], theme=tt1)
  rb2 <- tableGrob(reference_batches[21:nrow(reference_batches),1], theme=tt1)
  gc2 <- gtable_combine(rb1, rb2)
} else{
  gc2 <- tableGrob(reference_batches[,1], theme=tt1)
}


if (nrow(outlier_df)>0){
  outlier_message <- ''
} else{
  outlier_message <- '[No outliers detected]'
}

g.outlier <- textGrob(outlier_message, gp = gpar(col = "blue", fontsize = 20))
if (nrow(failed_ab)>0){
  gc3 <- tableGrob(failed_ab, theme=tt1)
  haligned <- gtable_combine(gc2, gc3)
  cover_sheet<- grid.arrange(gc1, haligned, g.outlier, ncol=1)
} else{
  cover_sheet<- grid.arrange(gc1, gc2, g.outlier, ncol=1)
}
grid.draw(cover_sheet)

# Draw outlier table
if (nrow(outlier_df>0)){
  for (j in seq(1,nrow(outlier_df), by=35)){
    g2 <- tableGrob(na.omit(outlier_df[j:(j+34), 2:7]), rows = NULL, theme = tt)
    g2 <- gtable_add_grob(g2, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                          t = 2, b = nrow(g2), l = 1, r = ncol(g2))
    g2 <- gtable_add_grob(g2, grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                          t = 1, l = 1, r = ncol(g2))

    grid.newpage()
    grid.draw(g2)
  }
}

for (tums in names(plot.list)){
  show(plot.list[[tums]])
}

# Create and write table of normalized counts.
score_out <- samp.scores %>% select(ProbeName,segment_label,norm)

# Separate scores by tumor/stroma segments
score_tum <- score_out[`segment_label` == 'tumor']
score_str <- score_out[`segment_label` == 'stroma']

# Insert percentile in counts table

if (nrow(score_tum)>0){
  percentiles <- sapply(unique(score_tum$ProbeName), function(x){
    v <- ref.abund[`segment_label` == 'tumor'] %>% filter(ProbeName==x) %>% select(norm)
    batch <- score_tum %>% filter(ProbeName == x) %>% select(norm)
    v <- rbind(v,batch)
    z <- (batch$norm[1]-mean(v$norm)) /sd(v$norm)
    percentile <- round(pnorm(z)*100)
    rank <- round(match(batch, sort(as.vector(v)[[1]]))/length(v$norm)*100, digits=2)
    percentile <- cbind(percentile, rank,
                        as.vector(round(quantile(v$norm, probs=0.25), digits=2)),
                        as.vector(round(quantile(v$norm, probs=0.5), digits=2)),
                        as.vector(round(quantile(v$norm, probs=0.75), digits=2)),
                        as.vector(round(quantile(v$norm, probs=1), digits=2))
                        )
  })
  percentiles <- t(percentiles)
  colnames(percentiles) <- c('z-score perc', 'rank perc', 'Q1', 'Q2', 'Q3', 'Q4')
  score_tum <- cbind(score_tum, percentiles)

  g <- tableGrob(score_tum[1:34,1:9], rows = NULL, theme = tt)
  g <- gtable_add_grob(g,
                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 2, b = nrow(g), l = 1, r = ncol(g))
  g <- gtable_add_grob(g,
                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 1, l = 1, r = ncol(g))

  grid.newpage()
  grid.draw(g)

  g1 <- tableGrob(score_tum[35:68,1:9], rows = NULL, theme = tt)
  g1 <- gtable_add_grob(g1,
                        grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                        t = 2, b = nrow(g1), l = 1, r = ncol(g1))
  g1 <- gtable_add_grob(g1,
                        grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                        t = 1, l = 1, r = ncol(g1))

  grid.newpage()
  grid.draw(g1)
}

if(nrow(score_str>0)){
  percentiles <- sapply(unique(score_str$ProbeName), function(x){
    v <- ref.abund[`segment_label` == 'stroma'] %>% filter(ProbeName==x) %>% select(norm)
    batch <- score_str %>% filter(ProbeName == x) %>% select(norm)
    v <- rbind(v,batch)
    z <- (batch$norm[1]-mean(v$norm)) /sd(v$norm)
    percentile <- round(pnorm(z)*100)
    rank <- round(match(batch, sort(as.vector(v)[[1]]))/length(v$norm)*100, digits=2)
    percentile <- cbind(percentile, rank,
                        as.vector(round(quantile(v$norm, probs=0.25), digits=2)),
                        as.vector(round(quantile(v$norm, probs=0.5), digits=2)),
                        as.vector(round(quantile(v$norm, probs=0.75), digits=2)),
                        as.vector(round(quantile(v$norm, probs=1), digits=2))
                        )
  })
  percentiles <- t(percentiles)
  colnames(percentiles) <- c('z-score perc', 'rank perc', 'Q1', 'Q2', 'Q3', 'Q4')
  score_str <- cbind(score_str, percentiles)

  # Draw to PDF
  g <- tableGrob(score_str[1:34,1:9], rows = NULL, theme = tt)
  g <- gtable_add_grob(g,
                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 2, b = nrow(g), l = 1, r = ncol(g))
  g <- gtable_add_grob(g,
                       grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                       t = 1, l = 1, r = ncol(g))

  grid.newpage()
  grid.draw(g)

  g1 <- tableGrob(score_str[35:68,1:9], rows = NULL, theme = tt)
  g1 <- gtable_add_grob(g1,
                        grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                        t = 2, b = nrow(g1), l = 1, r = ncol(g1))
  g1 <- gtable_add_grob(g1,
                        grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                        t = 1, l = 1, r = ncol(g1))

  grid.newpage()
  grid.draw(g1)
}

# Write out to Excel file
write.xlsx(score_tum, excel.out.tum)
write.xlsx(score_str, excel.out.str)

# Produce boxplots
failed_ab <- as.vector(failed_ab$'Failed Antibodies')
clia_abs <- unique(c(failed_ab,clia_abs))

for (i in seq(1,length(clia_abs), by=4)){

  print(clia_abs[i])
  q1.plot <- ggplot(data=ab.batches[ProbeName == clia_abs[i]], mapping=aes(x=name, y=abundance)) +
    geom_boxplot(outlier.shape=NA) + geom_jitter(mapping=aes(color=fac_batch),size=3, height=0, width=.15, show.legend = F) +
    geom_jitter(data=ab.batches.cur[ProbeName == clia_abs[i]], size=3, width=.15, height=0) +
    theme_bw() + xlab("") + ylab("log2 Abundance") + ggtitle(paste("Antibody: ", clia_abs[i])) +
    scale_x_discrete(guide = guide_axis(n.dodge = 3))

  if ((i+1)<=length(clia_abs)){
    print(clia_abs[i+1])
    q2.plot <- ggplot(data=ab.batches[ProbeName == clia_abs[i+1]], mapping=aes(x=name, y=abundance)) +
      geom_boxplot(outlier.shape=NA) + geom_jitter(mapping=aes(color=fac_batch),size=3, height=0, width=.15, show.legend = T) +
      geom_jitter(data=ab.batches.cur[ProbeName == clia_abs[i+1]], size=3, width=.15, height=0) +
      theme_bw() + xlab("") + ylab("log2 Abundance") + ggtitle(paste("Antibody: ", clia_abs[i+1])) +
      scale_x_discrete(guide = guide_axis(n.dodge = 3))
  }

  if ((i+2)<=length(clia_abs)){
    print(clia_abs[i+2])
    q3.plot <- ggplot(data=ab.batches[ProbeName == clia_abs[i+2]], mapping=aes(x=name, y=abundance)) +
      geom_boxplot(outlier.shape=NA) + geom_jitter(mapping=aes(color=fac_batch),size=3, height=0, width=.15, show.legend = F) +
      geom_jitter(data=ab.batches.cur[ProbeName == clia_abs[i+2]], size=3, width=.15, height=0) +
      theme_bw() + xlab("") + ylab("log2 Abundance") + ggtitle(paste("Antibody: ", clia_abs[i+2])) +
      scale_x_discrete(guide = guide_axis(n.dodge = 3))
  }

  if ((i+3)<=length(clia_abs)){
    print(clia_abs[i+3])
    q4.plot <- ggplot(data=ab.batches[ProbeName == clia_abs[i+3]], mapping=aes(x=name, y=abundance)) +
      geom_boxplot(outlier.shape=NA) + geom_jitter(mapping=aes(color=fac_batch),size=3, height=0, width=.15, show.legend = T) +
      geom_jitter(data=ab.batches.cur[ProbeName == clia_abs[i+3]], size=3, width=.15, height=0) +
      theme_bw() + xlab("") + ylab("log2 Abundance") + ggtitle(paste("Antibody: ", clia_abs[i+3])) +
      scale_x_discrete(guide = guide_axis(n.dodge = 3))
  }
  q1.grob <- ggplotGrob(q1.plot)
  q2.grob <- ggplotGrob(q2.plot)
  q3.grob <- ggplotGrob(q3.plot)
  q4.grob <- ggplotGrob(q4.plot)
  grid.arrange(q1.grob, q2.grob, q3.grob, q4.grob)
}

dev.off()