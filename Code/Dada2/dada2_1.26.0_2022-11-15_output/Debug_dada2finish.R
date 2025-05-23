# Debug dada2.finish


fastq.path = fastq.dir
truncFwd.len = 225
truncRev.len = 200
taxa.db = taxa.db.path
metadata.file = file.path(inDir, metadata.filename)
maxCores = 9
build.tree = TRUE
guide.seqs.file = NULL
alignment.template.file = NULL
fasttree.path = fasttree.path
user.output.path = NULL
paired = TRUE
force = FALSE
ask = TRUE
# ) 

# {
  if (length(names(run.env)) == 0) {
    rlang::abort("Functions 'initiate.pipeline() and 'dada2.upto.qualPlots()' must be run first.")
  }
  if (is.null(user.output.path)) {
    output <- run.env$output.path
  } else {
    output <- user.output.path
  }
  if (!any(file.exists(list.files(path = output, pattern = "qualPlot.pdf", full.names = T)))) {
    rlang::abort("Function 'dada2.upto.qualPlots()' must be run first.")
  }
  if (build.tree) {
    if (length(suppressWarnings(system("which mothur", intern = T))) == 0) {
      rlang::abort(
        "It appears you are trying to build a phylogenetic tree,\nbut mothur is not installed on your system.\nPlease install mothur and try again. If mothur is installed,\ncheck that it is in your PATH [try Sys.getenv('PATH')].\nSee function add2PATH() [?add2PATH] for adding a directory to your PATH."
      )
    }
    if (is.null(fasttree.path) | !file.exists(fasttree.path)) {
      rlang::abort(
        "It appears you are trying to build a phylogenetic tree, but you have not provided a viable path to FastTree."
      )
    }
    if (is.null(guide.seqs.file)) {
      writeLines(guide.seqs.lines, con = "guide_seqs.fasta")
      guide.seqs.file <- "guide_seqs.fasta"
    }
    if (is.null(alignment.template.file)) {
      writeLines(alignment.template.file.lines, con = "template.align")
      alignment.template.file <- "template.align"
    }
    if (!(
      guide.seqs.file %in% list.files() &
      alignment.template.file %in% list.files()
    )) {
      rlang::abort(
        paste(
          "Files", guide.seqs.file, "and", alignment.template.file,
          "must be in your current directory to build a tree."
        )
      )
    }
  }
  smpl.tbl <- read.file(metadata.file)
  if ("data.table" %in% class(smpl.tbl)) {
    smpl.col <- smpl.tbl[
      , sapply(
        .SD,
        function(col) {
          sum(col %in% run.env$sample.names) == length(run.env$sample.names)
        }
      )
    ] %>%
      extract(., .) %>%
      names() %>%
      extract(1)
    smpl.df <- as.data.frame(smpl.tbl)
    row.names(smpl.df) <- smpl.tbl[[smpl.col]]
  } else {
    smpl.df <- smpl.tbl
  }
  matching.names <- identical(
    sort(row.names(smpl.df)),
    sort(run.env$sample.names)
  )
  choice <- 1
  if (!matching.names) {
    if (ask) {
      rlang::inform(
        "The sample names supplied in the `metadata.file` do not match the sample names detected from sequence files during processing."
      )
      rlang::inform(
        "Do you to want to\n\t1) proceed with only those samples found in both the fastq files and metadata.file\n\t2) terminate processing"
      )
      choice <- readline(prompt = "Choice: ")
      while (!(choice %in% as.character(1:2))) {
        proceed <- readline(prompt = "Please choose 1 or 2: ")
      }
      choice <- as.integer(choice)
    } else {
      choice <- 1
    }
    if (choice == 1) {
      run.env$sample.names <- run.env$sample.names[run.env$sample.names %in% row.names(smpl.df)]
      smpl.df <- smpl.df[run.env$sample.names, , drop = F]
    }
  }
  if (choice == 2) {
    rlang::inform("Execution of dada2.finish has been terminated")
  } else {
    filtFs <- file.path(
      fastq.path,
      "Filtered",
      paste0(run.env$sample.names, "_F_filt.fastq.gz")
    )
    names(filtFs) <- run.env$sample.names
    if (paired) {
      filtRs <- file.path(
        fastq.path,
        "Filtered",
        paste0(run.env$sample.names, "_R_filt.fastq.gz")
      )
      names(filtRs) <- run.env$sample.names
    }
    
    out.file <- file.path(output, "filter_and_trim_numbers.rds")
    if (file.exists(out.file) & !force) {
      my.cat("Filtering and trimming completed previously, skipping...")
      out <- readRDS(out.file)
      my.cat("\tDONE:")
      print(head(out))
    } else {
      my.cat("Filtering and trimming...")
      if (paired) {
        out <- filterAndTrim(
          run.env$fnFs, filtFs,
          run.env$fnRs, filtRs,
          truncLen = c(truncFwd.len, truncRev.len),
          maxN = 0,
          maxEE = c(2, 2),
          truncQ = 2,
          rm.phix = TRUE,
          compress = TRUE,
          multithread = maxCores
        )
      } else {
        out <- filterAndTrim(
          run.env$fnFs, filtFs,
          truncLen = truncFwd.len,
          maxN = 0,
          maxEE = 2,
          truncQ = 2,
          rm.phix = TRUE,
          compress = TRUE,
          multithread = maxCores
        )
      }
      my.cat("\tDONE:")
      print(head(out))
      saveRDS(out, file = out.file)
    }
    
    
    filtFs <- filtFs[file.exists(filtFs)]
    if (paired) {
      filtRs <- filtRs[file.exists(filtRs)]
    }
    
    err.files <- file.path(output, "errF.rds")
    err.plot.files <- file.path(output, "errF_plot.pdf")
    if (paired) {
      err.files <- c(err.files, file.path(output, "errR.rds"))
      err.plot.files <- c(err.plot.files, file.path(output, "errR_plot.pdf"))
    }
    
    if (all(file.exists(err.files)) & !force) {
      my.cat("Errors learned previously, skipping...")
      errF <- readRDS(err.files[1])
      if (paired) {errR <- readRDS(err.files[2])}
    } else {
      my.cat("Learning errors and making error plots...")
      errF <- learnErrors(filtFs, multithread = maxCores)
      saveRDS(errF, file = err.files[1])
      errF.plot <- plotErrors(errF, nominalQ = TRUE)
      ggsave(errF.plot, file = err.plot.files[1])
      if (paired) {
        errR <- learnErrors(filtFs, multithread = maxCores)
        saveRDS(errR, file = err.files[2])
        errR.plot <- plotErrors(errR, nominalQ = TRUE)
        ggsave(errR.plot, file = err.plot.files[2])
      }
    }
    my.cat("\tDONE")
    
    dada.files <- file.path(output, "dadaFs.rds")
    if (paired) {
      dada.files <- c(dada.files, file.path(output, "dadaRs.rds"))
    }
    if (all(file.exists(dada.files)) & !force) {
      my.cat("dada-ing done previously, skipping...")
      dadaFs <- readRDS(dada.files[1])
      if (paired) {dadaRs <- readRDS(dada.files[2])}
    } else {
      my.cat("dada-ing...")
      dadaFs <- dada(filtFs, err = errF, multithread = maxCores)
      saveRDS(dadaFs, file = dada.files[1])
      if (paired) {
        dadaRs <- dada(filtRs, err = errR, multithread = maxCores)
        saveRDS(dadaRs, file = dada.files[2])
      }
    }
    my.cat("\tDONE")
    
    seqtab.file <- file.path(output, "seqtab.rds")
    if (paired) {
      mergers.file <- file.path(output, "mergers.rds")
    }
    
    if (file.exists(seqtab.file) & !force) {
      my.cat("Sequence table file exists, skipping...")
      seqtab <- readRDS(seqtab.file)
    } else {
      if (paired) {
        my.cat("Merging pairs...")
        mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
        merge.results <- sapply(mergers, getN)
        my.cat("\tDONE")
      } else {
        merge.results <- 1
      }
      
      if (sum(merge.results == 0) >= 3 | !paired) {
        if (paired) {
          my.cat("Multiple sequence merging fails, proceeding with FORWARD READS ONLY")
        }
        seqtab <- makeSequenceTable(dadaFs)
        saveRDS(seqtab, file = seqtab.file)
      } else {
        saveRDS(mergers, file = mergers.file)
        seqtab <- makeSequenceTable(mergers)
        saveRDS(seqtab, file = seqtab.file)
      }
    }
    
    seqtab.nochim.file <- file.path(output, "seqtab_nochim.rds")
    
    if (file.exists(seqtab.nochim.file) & !force) {
      my.cat("Bimeras previously removed, skipping...")
      seqtab.nochim <- readRDS(seqtab.nochim.file)
      my.cat(paste("\t\tPercent seqs kept:", sum(seqtab.nochim)/sum(seqtab)))
    } else {
      my.cat("Removing bimeras...")
      seqtab.nochim <- removeBimeraDenovo(
        seqtab,
        method = "consensus",
        multithread = maxCores,
        verbose = TRUE
      )
      my.cat("\tDONE")
      my.cat(paste("\t\tPercent seqs kept:", sum(seqtab.nochim)/sum(seqtab)))
      saveRDS(
        seqtab.nochim,
        file = seqtab.nochim.file
      )
    }
    track.file <- file.path(output, "track.rds")
    if (!file.exists(track.file)) {
      if (sum(merge.results == 0) >= 3) {
        track <- cbind(
          out,
          sapply(dadaFs, getN),
          sapply(dadaRs, getN),
          rowSums(seqtab.nochim)
        )
        colnames(track) <- c(
          "input",
          "filtered",
          "denoisedF",
          "denoisedR",
          "nonchimF"
        )
      } else {
        if (paired) {
          track <- cbind(
            out,
            sapply(dadaFs, getN),
            sapply(dadaRs, getN),
            merge.results,
            rowSums(seqtab.nochim)
          )
          colnames(track) <- c(
            "input",
            "filtered",
            "denoisedF",
            "denoisedR",
            "merged",
            "nonchim"
          )
        } else {
          track <- cbind(
            out,
            sapply(dadaFs, getN),
            rowSums(seqtab.nochim)
          )
          colnames(track) <- c(
            "input",
            "filtered",
            "denoisedF",
            "nonchimF"
          )
        }
      }
      rownames(track) <- run.env$sample.names
      saveRDS(track, file = track.file)
    } else {
      track <- readRDS(track.file)
    }
    my.cat(
      "See 'track.rds' in output directory for how many seqs made it through. First 6 samples look as follows:"
    )
    print(head(track))
    
    
    taxa.file <-  file.path(output, "taxa.rds")
    
    if (file.exists(taxa.file) & !force) {
      my.cat("Taxonomy previously assigned, skipping...")
      taxa <- readRDS(taxa.file)
    } else {
      my.cat("Assigning taxonomy...")
      taxa <- assignTaxonomy(
        seqtab.nochim,
        taxa.db,
        multithread = maxCores
      )
      saveRDS(taxa, file = taxa.file)
    }
    my.cat("\tDONE")
    
    ps0 <- phyloseq(
      otu_table(seqtab.nochim, taxa_are_rows = FALSE),
      sample_data(smpl.df),
      tax_table(taxa)
    )
    
    ps1 <- numbered.ASVs(
      ps = ps0,
      # prefix = paste0(proj.name, "_ASV"),
      save.dir = output,
      save.file = "asv_sequences"
    )
    asv.seqs <- readRDS(file.path(output, "asv_sequences.rds"))
    
    seqinr::write.fasta(
      sequences = as.list(asv.seqs),
      names = taxa_names(ps1),
      as.string = TRUE,
      file.out = file.path(output, "asv_sequences.fasta")
    )
    if (build.tree) {
      my.cat("Proceeding with phylogenetic tree:")
      asv.seqs.file <- file.path(output, "asv_sequences.fasta")
      asv.withguides.file  <- file.path(output, "asv_and_guide_seqs.fasta")
      asv.tree.rooted.file <- file.path(output, "asv_NASTaligned_seqs.nwk")
      
      cmd <- paste(
        "cat",
        asv.seqs.file,
        guide.seqs.file,
        ">",
        asv.withguides.file
      )
      system(cmd)
      
      my.cat("Aligning sequences...")
      cmd <- paste0(
        "mothur \"#align.seqs( fasta=",
        asv.withguides.file,
        ", reference=",
        alignment.template.file,
        ", flip=t",
        ", keepdots=t",
        ", processors=", maxCores,
        ", outputdir=",
        output,
        "/ )\""
      )
      system(cmd)
      my.cat("\tDONE")
      
      mothur.output.file <- file.path(output, "asv_and_guide_seqs.align")
      fasttree.log.file <- file.path(output, "fasttree.log")
      fasttree.output.file <- file.path(output, "asv_and_guide_seqs.nwk")
      
      my.cat("Building phylogenetic tree...")
      cmd <- paste0(
        "export OMP_NUM_THREADS=",
        maxCores, "; ",
        paste0("'", fasttree.path, "'"), " -nt -nosupport -quote -gtr -gamma -log ",
        fasttree.log.file,
        " ",
        mothur.output.file,
        " > ",
        fasttree.output.file
        
      )
      system(cmd)
      my.cat("\tDONE")
      asvs.and.guides.tree <- read_tree(fasttree.output.file)
      asvs.and.guides.tree.rooted <- phangorn::midpoint(asvs.and.guides.tree)
      
      guides <- scan(guide.seqs.file, what = "character" ) # Should be: "guide_seqs.fasta"
      
      guide.ids <- guides[stringr::str_detect(guides, ">" )]
      guide.ids <- stringr::str_remove(guide.ids, ">")
      
      asvs.tree.rooted <- ape::drop.tip(asvs.and.guides.tree.rooted, guide.ids)
      write.tree(asvs.tree.rooted, file = asv.tree.rooted.file)
      
      phy_tree(ps1) <- phy_tree(asvs.tree.rooted)
      system(paste("mv mothur*", output))
    }
    saveRDS(ps1, file = file.path(output, "phyloseq.rds"))
    for (file in c("tmp.txt", guide.seqs.file, alignment.template.file)) {
      if (file.exists(file)) { file.remove(file) }
    }
    my.cat("\tDONE and DONE")
  }
  