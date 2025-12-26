#=====================#
#### ngs functions ####
#=====================#

library(tidyverse)
library(magrittr)

#### row2colname ####
# a function when pivoted matrix needs a new colname from a row (number `n`)
# and convert columns (excluding .excol) to a certain class by using .fun function
row2colname <- function(df, n = 1, .fun = NULL, .excol){
  tmp <- 
    df %>% 
    magrittr::set_colnames(.[n,]) %>% 
    slice(-n)
  if(length(.fun) > 0){
    tmp %<>% 
      mutate_at(vars(-all_of(.excol)), ~.fun(.x))
  }
  tmp
}
#### end ####

#### functions:: adapter.detect ####
# A function to manually detect the full length of selected sequences in the raw fastq #
adapter.detect <- function(x, .type, .input){
  sub <- filter(x, type == .type)
  library <- str_c(sub$seq, collapse = "|")
  out <- str_locate_all(.input, library)[[1]] %>% as_tibble
  
  out <- 
    str_which(.input, sub$seq) %>% 
    dplyr::slice(sub, .) %>% 
    select(name, type, direction) %>% 
    bind_cols(out,.)
  
  out
}
#### end ####

#### functions:: adapter.check ####
# A function to manually check the full length of primer/adapter sequences in the raw fastq #
adapter.check <- function(path, size, ref, typelist, option = "all"){
  
  df <- 
    str_interp("zcat ${path} | head -n${size*4}") %>% 
    system(intern = TRUE) %>% 
    matrix(ncol = 4, byrow = TRUE) %>% 
    as_tibble() %>% 
    magrittr::set_colnames(c("name", "seq", "+", "qualityscore")) %>%
    mutate(id = 1:nrow(.)) %>% 
    select(!c(name, `+`, qualityscore)) 
  
  for(i in 1:length(typelist)){
    df %<>% 
      mutate(tmpcat = pmap(., function(seq, ...)adapter.detect(x=ref, .type = typelist[i], .input = seq)))
    
    colnames(df)[colnames(df) == "tmpcat"] <- typelist[i]
  }
  
    df %<>%
      unnest_longer(all_of(typelist), keep_empty = TRUE) %>% 
      unnest(cols = all_of(typelist), names_sep = ".")
  
  .colnames <- c("seq", "id", "start", "end", "name", "type", "direction")
  
  tmp <- vector("list", length(typelist))
  
  for(i in 1:length(typelist)){
    tmp[[i]] <- df %>% select(seq, id, starts_with(typelist[i])) %>% set_colnames(.colnames)
  }
  # tmp1 <- df %>% select(seq, id, starts_with("ssr")) %>% set_colnames(.colnames)
  # tmp2 <- df %>% select(seq, id, starts_with("pcrtail")) %>% set_colnames(.colnames)
  # tmp3 <- df %>% select(seq, id, starts_with("adapter")) %>% set_colnames(.colnames)
  
  if(option == "contaminantonly"){
    bind_rows(tmp) %>% 
      arrange(id) %>% 
      filter(!is.na(start)) %>% 
      distinct()
  }else{
    bind_rows(tmp) %>% 
      arrange(id) %>% 
      # filter(!is.na(start)) %>% 
      distinct()
  }
  
}
#### end ####

#### function:: md5sumcheck ####
# A function to conduct md5sumcheck #
md5check <- function(x){
  system(str_interp("md5sum ${x}"), intern = TRUE) %>% str_replace(" +.*","")
}
#### end ####

#### function:: dist2nexus ####
# A function to convert distance matrix to nexus format #

dist2nexus <- function(df,
                       labels = "no", 
                       triangle = "both",
                       diagonal = TRUE,
                       dir,
                       dryrun = TRUE){
  
  ind <- df %>% pull(ind) %>% length()
  taxlabels <- df %>% pull(ind) %>% str_c(collapse = "\n")
  
  if(labels == "no"){
    df <- df %>% select(!ind)  
  }
  m <- df %>% apply(1, str_c, collapse = " ") %>% str_c(collapse = "\n")
  
  header <- "#nexus"
  taxablock <- str_c("BEGIN Taxa;", 
                     str_interp("DIMENSIONS ntax=${ind};"), 
                     "TAXLABELS",
                     taxlabels,
                     ";","END; [Taxa]", 
                     sep = "\n")
  
  .diag <- ifelse(diagonal, "DIAGONAL", "NO DIAGONAL")
  
  distanceblock <- str_c("BEGIN Distances;",
                         str_interp("DIMENSIONS ntax=${ind};"),
                         str_interp("FORMAT labels=${str_to_upper(labels)} ${.diag} triangle=${str_to_upper(triangle)};"),
                         "MATRIX",
                         m,
                         ";",
                         "END; [Distances]",
                         sep = "\n")
  
  nexus <- str_c(header, taxablock, distanceblock, sep = "\n\n")
  
  if(dryrun){
    cat(nexus)
  }else{
    write_lines(nexus, dir)
  }
}
#### end ####

#### dna2geno ####
# previously called multiconvert
dna2geno <- function(x){
  g <- case_when(
    x %in% c("AA", "A") ~ "A",
    x %in% c("CC", "C") ~ "C",
    x %in% c("GG", "G") ~ "G",
    x %in% c("TT", "T") ~ "T",  
    x %in% c("AC", "CA") ~ "M",
    x %in% c("AG", "GA") ~ "R",
    x %in% c("AT", "TA") ~ "W",
    x %in% c("CG", "GC") ~ "S",
    x %in% c("CT", "TC") ~ "Y",
    x %in% c("GT", "TG") ~ "K",
    x %in% c("NN")       ~ "N"
  )
  
  g
}

#### end ####

#### countmap ####
# A function to view the proportion of paired mapped reads # 

countmap <- function(sample, directory){
  n.mappaired <- str_interp("samtools view -f 0x1 -F 0x904 -c ${directory}/${sample}.bam") %>%
    system(intern = TRUE) %>% 
    as.numeric
  n.mapunpaired <- str_interp("samtools view -f 0x8 -F 0x904 -c ${directory}/${sample}.bam") %>% 
    system(intern = TRUE) %>% 
    as.numeric()
  
  data.frame(n.mappaired, n.mapunpaired)
}
#### end #### 

#### readcounter #### 
# A function to count the number of reads #
readcounter <- function(fullpath.fastq){
  str_interp("zcat ${fullpath.fastq}|wc -l") %>% 
    system(intern = TRUE) %>% 
    as.numeric(.)/4
}
#### end ####

#### group_consecutive ####
group_consecutive <- function(vec){
  dt <- data.table(value = vec)
  dt[, group := cumsum(c(1,diff(value) != 1))]
  dt[, .(range = paste0(min(value), "-", max(value))), by = group][, range]
}
#### end ####

#### pcaplot ####
pcaplot <- function(xaxis, yaxis, .group, .pc, .e, .colpal, .labels, .title = NULL){
  
  # combine into one dataframe
  plot.pc <- 
    .pc %>% 
    dplyr::select(all_of(c(xaxis, yaxis))) %>% 
    mutate(group = .group) %>% 
    set_colnames(c("xaxis","yaxis","group"))
  
  # fix the xaxis range always to have positives as larger absolute values
  xaxis.range <- range(plot.pc$xaxis)
  
  if(abs(xaxis.range[1]) > abs(xaxis.range[2])){
    plot.pc$xaxis <- -plot.pc$xaxis
  }
  
  # fix the yaxis range always to have positives as larger absolute values
  yaxis.range <- range(plot.pc$yaxis)
  
  if(abs(yaxis.range[1]) > abs(yaxis.range[2])){
    plot.pc$yaxis <- -plot.pc$yaxis
  }
  
  # percentage of the total variance
  plot.pct <-
    .e %>% 
    filter(pc %in% c(xaxis, yaxis)) %>% 
    pull(pct)
  
  # plot label
  plot.axis <- str_to_upper(c(xaxis, yaxis))
  
  # plot
  p <-
    ggplot(data = plot.pc, aes(x = xaxis, y = yaxis, color = group))+
      geom_hline(aes(yintercept = 0.00)) + 
      geom_vline(aes(xintercept = 0.00)) + 
      geom_point() + 
      # ggrepel::geom_label_repel(aes(label = sample, color = group_reg), max.overlaps = 15) +
      xlab(str_interp("${plot.axis[1]} (${plot.pct[1]}%)")) + 
      ylab(str_interp("${plot.axis[2]} (${plot.pct[2]}%)")) + 
      stat_ellipse(level=0.5) + 
      scale_color_manual(values = .colpal, labels = .labels) +
      theme_bw() 
  
  if(!is.null(.title)){
    p <- p + ggtitle(.title)
  }
  
  p
} 

#### chisquare test ####
qqchi<-function(x,...){
  lambda<-round(median(x)/qchisq(0.5,1),2)
  qqplot(qchisq((1:length(x)-0.5)/(length(x)),1),x,ylab="Observed",xlab="Expected",...);abline(0,1,col=2,lwd=2)
  legend("topleft",paste("lambda=",lambda))
}
#### function end ####

#### CATGformat ####
# A function to create a catg input file for raxml-ng
CATGformat <- function(.glf, .geno, .taxon,
                       directory, out = "input"){
  longornot <- (nrow(.geno)> 20000)
  
  if(longornot){
    sepnum <- floor(nrow(geno)/4)
    .geno <- split(.geno, (seq(nrow(.geno))-1) %/% sepnum)
    .glf <- split(.glf, (seq(nrow(.glf))-1) %/% sepnum)
  }else{
    .glf <- list(.glf)
    .geno <- list(.geno)
  }
  
  catg <- vector("list", length(.geno))
  
  for(k in 1:length(.geno)){
    # concatenate glf and genotype
    geno.concat <- 
      .geno[[k]] %>% 
      mutate_at(vars(starts_with("V")), dna2geno) %>% 
      unite("concat", starts_with("V"), sep="")
    
    # genotype order for 
    # glf
    geno.abbr <- c("A", "M", "R", "W", "C", "S", "Y", "G", "K", "T")
    # catg format
    geno.abbr.level <- c("A", "C", "G", "T", "M", "R", "W", "S", "Y", "K")
    
    # creating CATG
    colnam <- 
      expand.grid(str_c("geno", geno.abbr), .taxon) %>% 
      unite("colnam", sep = ".") %>% 
      pull(colnam)
    
    tmp.catg <-
      left_join(geno.concat, .glf[[k]], by = c("scaf", "pos")) %>% 
      select(!c(scaf, pos)) %>% 
      set_colnames(c("concat", colnam)) 
    
    # normalization of likelihoods
    # normalization required to make the sum of all the possible genotypes to 1
    # or otherwise the raxml-ng will crash
    
    tmp.catg.long <- 
      tmp.catg %>% 
      mutate(sitename = 1:nrow(.)) %>% 
      pivot_longer(starts_with("geno"), names_to = c("genotype","taxon"), names_pattern = "geno(.).(.*)") 
    
    sum.long <- 
      tmp.catg.long %>% 
      select(!concat) %>% 
      group_by(sitename, taxon) %>% 
      summarize(sum = sum(value)) %>% 
      select(sitename, taxon, sum)
    
    catg[[k]] <- 
      tmp.catg.long %>% 
      left_join(., sum.long, by = c("sitename", "taxon")) %>% 
      mutate(pct = if_else(sum != 0, value/sum, 0)) %>% 
      select(concat, sitename, genotype, taxon, pct) %>% 
      mutate(genotype = factor(genotype, levels = geno.abbr.level)) %>% 
      arrange(sitename, taxon, genotype) %>% 
      mutate(genotype = as.character(genotype)) %>% 
      pivot_wider(names_from = c(genotype),
                  names_glue = "geno{genotype}",
                  values_from = pct) %>% 
      unite("glued.like", str_c("geno", geno.abbr.level), sep = ",") %>% 
      pivot_wider(names_from = c(taxon),
                  values_from = glued.like) %>% 
      select(!sitename) 
    
    rm(list = c("tmp.catg", "tmp.catg.long", "sum.long"));gc(reset=TRUE)
  }
  
  if(longornot){
    catg <- bind_rows(catg)
    }else{
      catg <- catg[[1]]
    }
  
  file.catg <-
    c(str_c(ncol(catg)-1, nrow(catg), sep = " "),
      str_c(.taxon, collapse = "\t"),
      catg %>% 
        unite("concat_sample1", concat, .taxon[1], sep = "\t") %>% 
        apply(., 1, str_c, collapse = " ") %>% 
        str_c(collapse = "\n")
    )%>% 
    str_c(collapse = "\n")
  
  write_lines(file.catg, str_interp("${directory}/${out}.txt"))
  
  catg
}

#### end ####

#### CATGcombiner ####
catgcombiner <- function(list, out){
  # collect header information from multiple files
  # number of samples and sites
  df <- 
    lapply(list, function(x){
      read_lines(x, n_max = 1) %>% 
        str_split(" ", simplify = TRUE) %>% 
        as_tibble %>% 
        mutate_all(as.integer) %>% 
        mutate(input = x)
    }) %>% bind_rows
  num <- df$V2 %>% sum
  ind <- df$V1 %>% unique
  
  # sample name header
  head <- lapply(list, function(x){
    read_lines(x, skip = 1, n_max = 1)
  }) %>% unlist %>% unique
  
  # tmp files
  tmp.list <- str_replace(list, ".txt", "_tmp.txt")
  tmp.header <- str_replace(out, ".txt", "_header.txt")
  
  # while checking the header information correct, get rid of the header from each file
  # and concatenate them into one
  if(length(ind) > 1){
    cat("ERROR:: Number of samples are not equal among the listed files")
  }else if(length(head) > 1){
    cat("ERROR:: Sample names are not the same among the listed files")
  }else{
    # remove headers
    rmh <- mapply(function(x,y){str_interp("tail -n +3 ${x} > ${y}")}, list, tmp.list) %>% as.character()
    
    for(k in 1:length(rmh)){
      call <- rmh[k]
      system(call)
    }
    
    # save the header
    new.head <- str_interp('${ind} ${num}\n${head}')
    write_lines(new.head, tmp.header)
    # concatenate multiple of the listed files and header
    tmp.list <- c(tmp.header, tmp.list)
    str_c(tmp.list, collapse = " ") %>% str_c("cat ", ., " > ", out) %>% system()
    # remove files
    str_c(tmp.list, collapse = " ") %>% str_c("rm ", .) %>% system()
  }
}
#### end ####

#### gp2ambig.dt :: catg formatting by using data.table ####

pos.lib <- tibble(alleles      = c("A/C", "A/G", "A/T", "C/G", "C/T", "G/T",
                                   "C/A", "G/A", "T/A", "G/C", "T/C", "T/G"),
                  position_ref = c("1,5,2", "1,6,3", "1,7,4", "2,8,3", "2,9,4", "3,10,4",
                                   "2,5,1", "3,6,1", "4,7,1", "3,8,2", "4,9,2", "4,10,3"))

ambig.lib <- tibble(
  alleles  = c("A/A", "C/C", "G/G", "T/T", 
               "A/C", "A/G", "A/T", "C/G", "C/T", "G/T",
               "C/A", "G/A", "T/A", "G/C", "T/C", "T/G", 
               "."),
  genotype = c("A", "C", "G", "T", 
               "M", "R", "W", "S", "Y", "K",
               "M", "R", "W", "S", "Y", "K", 
               "N")
)

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

gp2ambig.dt <- function(gt_GP, gt_GT, position_ref, .decimal = 15, na_val = 0.00, ...){
  
  if(is.na(gt_GT)){
    p <- str_c(rep(as.character(na_val),10),collapse =",")
  }else{
    v <- data.table(pos = 1:10, val = "0.0")
    val <- 
      c(position_ref, gt_GP) %>% 
      str_split(",", simplify = TRUE) %>% 
      t %>% 
      as.data.table() %>% 
      .[, .(pos = as.integer(V1), GP = as.numeric(V2))] %>% 
      .[, GP := as.numeric(specify_decimal(GP, .decimal))] 
    
    val[GP == max(GP), GP := 1 - sum(GP[GP != max(GP)])]
    
    val[, GP := format(GP, scientific = FALSE)]
    
    v <- merge(v, val, by = "pos", all.x = TRUE)
    v[, val := fifelse(is.na(GP), val, GP)]  # valの置換
    v[, val.num := as.numeric(val)]
    v[, val := fifelse(val.num == 0, "0.0", val)]
    
    p <- paste(v$val, collapse = ",")
  }
  
  p
}
#### end ####

#### TreeMixformat ####
# Treemix input file creater
TreeMixformat <- function(df, target, outgroup, prefix = "JAKCOX", out.path, out.prefix){
  
  colnam <- df %>% select(starts_with(prefix)) %>% colnames()
  freq <- vector("list", length(target))
  
  # calculate per-population frequency of bases for each position
  for(i in 1:length(target)){
    
    freq[[i]] <- 
      df %>% 
      filter(str_detect(pop, target[i])) %>% 
      t %>% 
      as.data.frame %>% 
      format %>% 
      apply(., 1, paste, collapse = "") %>% 
      lapply(., str_count, c("A", "T", "C", "G")) %>% 
      bind_rows %>% 
      select(!c(sample, pop)) %>% 
      as_tibble %>% 
      mutate(base = c("A", "T", "C", "G"), pop = target[i]) %>% 
      select(base, pop, starts_with(prefix))
  }
  
  # create major/minor library and count each in each population
  majmin <- data.frame(position = colnam, maj = NA_character_, min = NA_character_, alt = NA_character_)
  freq <- bind_rows(freq)
  treemix <- target
  
  for(i in 1:(ncol(freq)-2)){
    
    tmp <- 
      freq[,c(1,2,i+2)] %>% 
      filter(.[3] != 0)
    
    pos <- colnames(tmp)[3]
    
    base <- tmp %>% pull(base) %>% unique
    
    majmin[majmin$position == pos,][,2:4] <- if(length(base) > 2){base}else{c(base, rep(NA_character_, 3-length(base)))}
  }
  
  # remove triallelic
  majmin <- filter(majmin, is.na(alt))
  
  for(i in 1:nrow(majmin)){
    
    pos <- majmin$position[i]
    
    base.tmp <- 
      df %>% 
      select(pop, matches(pos)) %>%
      magrittr::set_colnames(c("pop", "pos")) %>% 
      group_by(pop) %>%
      summarise(bases = str_c(pos, collapse ="")) %>% 
      ungroup() %>% 
      mutate(count1 = str_count(bases, as.character(majmin[i,2])),
             count2 = str_count(bases, as.character(majmin[i,3]))) %>% 
      mutate(base.freq = str_c(count1, count2, sep = ",")) 
    
    treemix <- rbind(treemix, base.tmp %>% pull(base.freq))
    rownames(treemix)[i+1] <- pos
  }
  
  out.filt <- str_c("is.na(", outgroup, ")") %>% str_c(collapse = " | ") 
  
  treemix %<>% 
    as.data.frame %>% 
    magrittr::set_colnames(.[1,]) %>% dplyr::slice(-1) %>% 
    filter(!(eval(parse(text = out.filt)))) %>% 
    filter_all(., all_vars(. != "0,0"))
  
  if(!str_detect(out.prefix, "^_") && nchar(out.prefix) > 0){out.prefix <- str_c("_", out.prefix)}
  
  write.table(treemix, file = str_interp("${out.path}/input${out.prefix}.freq"), 
              row.names = FALSE, col.name=TRUE, quote=FALSE)
  
  treemix
  
}

#### end ####

#### site2bed ####
site2bed <- function(.sitefile, dir = NULL, out = NULL){
  
  tmp <- 
    .sitefile %>%
    arrange(scaf, pos) %>% 
    group_by(scaf) %>% 
    mutate(diff = pos - lag(pos)) %>% 
    mutate(diff = if_else(diff != 1 | is.na(diff), 1, 0)) %>% 
    ungroup %>% 
    mutate(group = cumsum(diff)) %>% 
    group_by(group) %>% 
    mutate(start = min(pos)-1,
           end = max(pos)) %>% 
    ungroup %>% 
    filter(diff == 1) 
  
  bed <- 
    tmp %>% 
    select(scaf, start, end)
  
  if(length(dir) > 0){
    write_tsv(bed, str_interp("${dir}/${out}.bed"), col_names = FALSE)
    bgzip <- str_interp("bgzip -f ${dir}/${out}.bed")
    idx <- str_interp("tabix -s 1 -b 2 -e 3 ${dir}/${out}.bed.gz")
    system(str_c(bgzip, idx, sep = ";"))
  }
  
  bed
}
#### end ####

#### fscgen.tpl ####
fscgen.tpl <-
  function(pop = c(0,1,2), popsize = NULL, 
           growth = NULL, migmat = NULL, event.mat = NULL,
           sfs = TRUE, numindloc = NULL, numlkblok = NULL,
           rateparam = NULL,
           dir, out, dryrun = TRUE
           ){
    
    # headers
    headall <- "// Parameters for the coalescence simulation program : fastsimcoal2 in linux"
    n <- length(pop)
    HEAD <- str_interp("${n} samples to simulate")
    
    head1 <- "//Population effective sizes (number of genes)"
    # popname
    POPNAME <- str_c("NPOP", pop) %>% str_c(collapse = "\n")
    
    head2 <- "//Samples sizes and samples age"
    # popsize
    POPSIZE <- popsize %>% str_c(collapse = "\n")
    
    head3 <- "//Growth rates: negative growth implies population expansion"
    # Growth rates
    GROWTH <- growth %>% str_c(collapse = "\n")
    
    head4 <- "//Number of migration matrices : 0 implies no migration between demes"
    # number of migration matrices
    NMIGMAT <- length(migmat)
    
    MIGMAT <- vector("list", NMIGMAT)
    for(m in 1:NMIGMAT){
      MIGMAT[[m]] <- 
        str_interp("//Migration matrix ${m-1}") %>% 
        str_c(., "\n", migmat[[m]] %>% apply(., 1, str_c, collapse = " ") %>% str_c(collapse = "\n"))
    }
    
    MIGMAT <- MIGMAT %>% str_c(collapse = "\n")
    
    head5 <- "//historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index"
    nev <- nrow(event.mat)
    NEVENT <- str_interp("${nev} historical event")
    maxchr <- event.mat %>% apply(., 2,function(x){str_width(x) %>% max()})
    EVENTMAT <-
      apply(event.mat, 
            1, 
            function(x){map2(.x = x, .y = maxchr, .f = function(chr = .x, mxchr = .y){str_pad(chr, mxchr, "right")}) %>% unlist}) %>% 
      t() %>% 
      apply(1, str_c, collapse = " ") %>% 
      str_c(collapse = "\n")
    
    head6 <- "//Number of independent loci [chromosome]"
    head7 <- "//Per chromosome: Number of linkage blocks"
    if(sfs){
      NUMINDLOC <- "1 0"
      NUMLKBLOK <- "1" 
    }else{
      NUMINDLOC <- numindloc
      NUMLKBLOK <- numlkblok
    }
    
    head8 <- "//per Block: data type, num loci, rec. rate and mut rate + optional parameters"
    RATE <- rateparam
    
    tpl <- 
      str_c(headall, HEAD, head1, POPNAME, head2, POPSIZE, head3, GROWTH, 
            head4, NMIGMAT, MIGMAT, head5, NEVENT, EVENTMAT, head6, NUMINDLOC, 
            head7, NUMLKBLOK, head8, RATE, sep = "\n")
    
    if(!dryrun){
      write_lines(tpl, str_interp("${dir}/${out}.tpl"))
    }else{
      cat(tpl)
    }
    
  }

#### end ####

#### fscgen.mod ####
fscgen.mod <- function(tmpl, .mig,
                       merge = TRUE,
                       denom.nancall = NULL, denom.nancsub = NULL,
                       bounded.reference = TRUE,
                       dir, out, dryrun = TRUE
                       ){
  
  # migration parameters to RETAIN (TRUE) or REMOVE (FALSE)
  if(.mig == "mig"){
    mig02 = TRUE; mig01 = TRUE; mig12 = TRUE; mig10 = TRUE; mig20 = TRUE; mig21 = TRUE; migBC = FALSE; migCB = FALSE
  }
  if(.mig == "nonmig"){
    mig02 = FALSE; mig01 = FALSE; mig12 = FALSE; mig10 = FALSE; mig20 = FALSE; mig21 = FALSE; migBC = FALSE; migCB = FALSE
  }
  if(.mig == "diffmig"){
    mig02 = TRUE; mig01 = TRUE; mig12 = TRUE; mig10 = TRUE; mig20 = TRUE; mig21 = TRUE; migBC = TRUE; migCB = TRUE
  }
  if(.mig == "pastmig"){
    mig02 = FALSE; mig01 = FALSE; mig12 = FALSE; mig10 = FALSE; mig20 = FALSE; mig21 = FALSE; migBC = TRUE; migCB = TRUE
  }
  if(.mig == "constmig"){
    mig02 = TRUE; mig01 = TRUE; mig12 = TRUE; mig10 = TRUE; mig20 = TRUE; mig21 = TRUE; migBC = FALSE; migCB = FALSE
  }
  
  # FALSE to be removed, TRUE to be retained
  # those sections are tagged by [BEGIN/END]
  start.int <- str_which(tmpl, "\\[BEGIN")
  start.chr <- str_extract(tmpl, "\\[BEGIN.*") %>% .[!is.na(.)]
  end.int   <- str_which(tmpl, "\\[END")
  
  df.param <-
    tibble(start.int = start.int,
           start.chr = start.chr,
           end.int = end.int) %>% 
      mutate(defexp = str_extract(start.chr, "(DEF|EXP).*(?=_)"),
             param  = str_extract(start.chr, "(?<=_)[^_]*(?=\\]$)")) %>% 
      select(start.int, end.int, defexp, param) %>% 
    mutate(retain = rep(c(merge, mig02, mig01, mig12, mig10, mig20, mig21, migBC, migCB),2)) %>% 
    mutate(range = pmap(., function(start.int, end.int, retain, ...){
             if(retain){
               c(start.int, end.int)
             }else{
               start.int:end.int
             }
           })
    ) 
  
  # line numbers to remove  
  rm.line <- 
    df.param %>% 
    pull(range) %>% 
    unlist
  
  # remove  
  tmpl.edit <- tmpl[-rm.line]
  
  # parameter overwriting - expansion ratio and migration sink population
  tmpl.edit %<>% str_replace_all(., "\\[basepop1\\]", str_c("NPOP", denom.nancall))
  tmpl.edit %<>% str_replace_all(., "\\[basepop2\\]", str_c("NPOP", denom.nancsub))
  
  # parameter overwriting - bounded reference
  if(bounded.reference){
    tmpl.edit %<>% str_replace_all(., "\\[lTDIV2\\]", "100000")
    tmpl.edit %<>% str_replace_all(., "\\[hTDIV2\\]", "1000000")
    tmpl.edit %<>% str_replace_all(., "\\[bounded reference\\]", "bounded reference")
  }else{
    tmpl.edit %<>% str_replace_all(., "\\[lTDIV2\\]", "100000")
    tmpl.edit %<>% str_replace_all(., "\\[hTDIV2\\]", "1000000")
    tmpl.edit %<>% str_replace_all(., "\\[bounded reference\\]", "")
  }
  
  # save/print
  if(!dryrun){
    write_lines(tmpl.edit, str_interp("${dir}/${out}.est"))
  }else{
    cat(tmpl.edit %>% str_c(collapse = "\n"))
  }
 
}
#### end ####

#### fscgen.hist ####
fscgen.hist <- function(model, .mig, .popA, .popB, .popC){
  
  if(.mig == "mig"){
    migindex1 <- 1
    }else if(.mig == "nonmig"){
      migindex1 <- 0
    }else if(.mig == "diffmig"){
      migindex1 <- 1
    }else if(.mig == "pastmig"){
      migindex1 <- 1
    }else if(.mig == "constmig"){
      migindex1 <- 1
    }

  if(.mig == "mig"){
    migindex2 <- 1
    }else if(.mig == "nonmig"){
      migindex2 <- 0
    }else if(.mig == "diffmig"){
      migindex2 <- 2
    }else if(.mig == "pastmig"){
      migindex2 <- 0
    }else if(.mig == "constmig"){
      migindex2 <- 2
    }
  
  # migindex1 <- ifelse(.mig == "mig", 1, ifelse(.mig == "nonmig", 0, 1))
  # migindex2 <- ifelse(.mig == "mig", 1, ifelse(.mig == "nonmig", 0, 2))
  
  # //historical event: time, source, sink, migrants, new deme size, growth rate, migr mat index
  event.base1 <- 
    matrix(c("TDIV1", .popA, .popB, 1, 1,      0, migindex1,
             "TDIV2", .popB, .popC, 1, 1,      0, migindex2), nrow = 2, byrow = TRUE)
  # this is same as
  # event.base1 <- 
  #   matrix(c("TDIV1", .popB, .popB, 1, "RES1", 0, migindex1,
  #            "TDIV1", .popA, .popB, 1, 1,      0, migindex1,
  #            "TDIV2", .popC, .popC, 1, "RES0", 0, migindex2,
  #            "TDIV2", .popB, .popC, 1, 1,      0, migindex2,), nrow = 4, byrow = TRUE)
  
  event.base2 <- 
    matrix(c("TDIV1", .popA, .popB, 1, "RES1", 0, migindex1,
             "TDIV2", .popB, .popC, 1, "RES0", 0, migindex2), nrow = 2, byrow = TRUE)
  
  event.exp <-
    matrix(c("TEXP", .popB, .popB, 0, "RES2", 0, ifelse(.mig == "mig", 1, 0),
             "TEXP", .popA, .popA, 0, "RES3", 0, ifelse(.mig == "mig", 1, 0),
             "TEXP", .popC, .popC, 0, "RES4", 0, ifelse(.mig == "mig", 1, 0)), nrow = 3, byrow = TRUE)
  
  if(model == "base1"){
    rt.model <- event.base1
  }
  
  if(model == "base2"){
    rt.model <- event.base2
  }
  
  if(str_detect(model, "exp")){
    rt.model <- list(event.exp, rt.model)
  }
  
  rt.model
}

#### end ####

#### fscgen.migmat ####
fscgen.migmat <- function(.mig, .popA, .popB, .popC, mig.mat0, mig.mat1){
  
  mig.matA <- mig.mat0
  mig.matA[mig.mat0 == str_c("MIG", .popB, .popC)] <- "MIGBC"
  mig.matA[mig.mat0 == str_c("MIG", .popC, .popB)] <- "MIGCB"
  mig.matA[!(mig.mat0 == str_c("MIG", .popB, .popC)|mig.mat0 == str_c("MIG", .popC, .popB))] <- "0"
  
  mig.matB <- mig.mat0
  mig.matB[!((mig.mat0 == str_c("MIG", .popB, .popC))|mig.mat0 == str_c("MIG", .popC, .popB))] <- "0"
  # migmat <- case_when(.mig == "mig"     ~ list(mig.mat0, mig.mat1),
  #                     .mig == "nonmig"  ~ list(mig.mat1), 
  #                     .mig == "diffmig" ~ list(mig.mat0, mig.matA, mig.mat1),
  #                     .mig == "pastmig" ~ list(mig.mat1, mig.matA))
  
  if(.mig == "mig"){
    migmat <- list(mig.mat0, mig.mat1)
  }else if(.mig == "nonmig"){
    migmat <- list(mig.mat1)
  }else if(.mig == "diffmig"){
    migmat <- list(mig.mat0, mig.matA, mig.mat1)
  }else if(.mig == "pastmig"){
    migmat <- list(mig.mat1, mig.matA)
  }else if(.mig == "constmig"){
    migmat <- list(mig.mat0, mig.matB, mig.mat1)
  }
  migmat
}
#### end ####

#### fscrun ####
fscrun <- function(split = TRUE, start, end, splitn = 0, nsim = "100000", optlECM = 40, opthECM = 100, typesfs, sgltn, 
                   memsz = 3, cpunum = 72, runname, dir, subdir, subdir.clust){
  
  if(split){
    ############ TEST SHORT SIMULATION OR SPLIT LOOP SIMULATION
    # parameters
    
    if(splitn > 0){
      list.split <- split(start:end, cut(seq_along(start:end), splitn, labels=FALSE))
    }else{
      list.split <- list(start:end)
    }
    

    for(i in 1:length(list.split)){
      tmp.start <- list.split[[i]][1]
      tmp.end <- list.split[[i]][length(list.split[[i]])]
      
      call <- str_c(
        str_interp("for chain in {${tmp.start}..${tmp.end}}"),
        str_interp("do  mkdir ${dir}/run$chain"),
        str_interp("cp $(ls ${dir}/*${runname}*) ${dir}/run$chain\"/\""),
        str_interp("cd ${dir}/run$chain"),
        str_interp("/lfs/aokid/fsc27093 -t ${runname}.tpl -e ${runname}.est ${typesfs} -C${sgltn} -n${nsim} -s0 -j -M0.01 -l${optlECM} -L${opthECM} -q --seed1234 -c72 -B72"),
        "cd ..",
        "done", sep = "\n")
      
      cluster.scripter(filename = str_interp("fsc_${runname}_run${tmp.start}-${tmp.end}"),
                       local.directory = subdir,
                       argument = list(call),
                       name = str_interp("fsc_${runname}_r${tmp.start}-${tmp.end}"),
                       memsz = memsz, cpunum = cpunum, envn = "bioinfo",
                       cluster.directory = subdir.clust)
    }
      

    ############ TEST SHORT SIMULATION OR SPLIT LOOP SIMULATION END 
  }else{
    
    for(chain in start:end){ # 1:100

      call.mkdir <- str_interp("mkdir ${dir}/run${chain}")
      call.cp   <-  str_interp("cp $(ls ${dir}/*${runname}*) ${dir}/run${chain}\"/\"")
      call.cd   <- str_interp("cd ${dir}/run${chain}")
      call.run <- str_interp("/lfs/aokid/fsc27093 -t ${runname}.tpl -e ${runname}.est ${typesfs} -C${sgltn} -n${nsim} -s0 -j -M0.01 -l${optlECM} -L${opthECM} -q --seed1234 -c72 -B72")
      
      cluster.scripter(filename = str_interp("fsc_${runname}_run${chain}"),
                       local.directory = subdir,
                       argument = list(call.mkdir, call.cp, call.cd, call.run),
                       name = str_interp("fsc_${runname}_r${chain}"),
                       memsz = memsz, cpunum = cpunum, envn = "bioinfo",
                       cluster.directory = subdir.clust)
    }
  }
  
  
  
}
#### end ####

#### LDdecay ####
# This is a function to fit a LD-decay curve
LDdecay <- function(df, Cstart, rho, fit){
  if(fit){
    nls(LD ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))), 
        data = data.frame(dist=df$dist, LD=df$r2), 
        start = Cstart, 
        control=nls.control(maxiter=100)) 
  }else{
    ((10+rho*df$dist)/((2+rho*df$dist)*(11+rho*df$dist)))*(1+((3+rho*df$dist)*(12+12*rho*df$dist+(rho*df$dist)^2))/(n*(2+rho*df$dist)*(11+rho*df$dist)))
  }
  
}
#### end ####

#### SFStidy ####
SFStidy <- function(path, popfun, type){
  if(type == "obs"){
    colnam <- read_lines(path, skip = 1, n_max = 1) %>% str_split("\t") %>% .[[1]] %>% .[-1]
    body <- read_lines(path, skip = 2)
    rownam <- body %>% str_split("\t", simplify = TRUE) %>% .[,1] 
    mat <- 
      body %>% 
      str_split("\t", simplify = TRUE) %>% 
      .[,2] %>% 
      str_split(., "(?<=[:digit:]) +", simplify = TRUE) %>% 
      as.data.frame %>% 
      mutate_all(as.integer) %>% 
      as.matrix
    
    colnames(mat) <- colnam
    rownames(mat) <- rownam
    
    tmp <- 
      tibble(diff1 = colnames(mat)[col(mat)], diff2 = rownames(mat)[row(mat)], site = c(mat)) %>% 
      mutate(pop1 = str_extract(diff1, "(?<=d).") %>% popfun,
             pop2 = str_extract(diff2, "(?<=d).") %>% popfun) %>% 
      mutate(diff1 = str_remove(diff1, "d._") %>% as.integer,
             diff2 = str_remove(diff2, "d._") %>% as.integer) 
    maxdiff1 <- tmp$diff1 %>% max
    maxdiff2 <- tmp$diff2 %>% max
    
    tmp %<>% mutate(sum=diff1+diff2) %>% filter(sum <= maxdiff1+maxdiff2) %>% select(!sum)
    
  }else{
    tmp <- read_tsv(path) %>% as.data.frame %>% 
      column_to_rownames(var = "...1") %>% 
      as.matrix
    tmp <- tibble(diff1 = colnames(tmp)[col(tmp)], diff2 = rownames(tmp)[row(tmp)], site = c(tmp)) %>% 
      mutate(pop1 = str_extract(diff1, "(?<=d).") %>% popfun,
             pop2 = str_extract(diff2, "(?<=d).") %>% popfun) %>% 
      mutate(diff1 = str_remove(diff1, "d._") %>% as.numeric,
             diff2 = str_remove(diff2, "d._") %>% as.numeric)
  }
 
  # model name
  model <- str_extract(path, "[^/]+(?=/run[:digit:]{1,3}/)")
  
  tmp %>% mutate(model = model, type = type)
  
}





