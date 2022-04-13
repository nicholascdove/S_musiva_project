
#' read qiime2 metadata (.tsv)
#'
#' Loads a qiime2 metadata file wherein the 2nd line contains the #q2:types line dictating the type of variable (categorical/numeric)
#'
#' @param file path to the input file, ex: file="~/data/moving_pictures/table.qza"

#' @return a data.frame wherein the first column is SampleID
#'
#' @examples \dontrun{metadata<-read_q2metadata("q2metadata.tsv")}
#' @export
#'
#'


read_q2metadata <- function(file) {
  if(missing(file)){stop("Path to metadata file not found")}
  
  defline<-suppressWarnings(readLines(file)[2])
  if(!grepl("^#q2:types", defline)){stop("Metadata does not define types (ie second line does not start with #q2:types)")}
  
  defline<-strsplit(defline, split="\t")[[1]]
  
  defline[grep("numeric", defline)]<-"as.numeric"
  defline[grep("categorical|q2:types", defline)]<-"as.factor"
  
  metadata<-subset(read.table(file, header=F, comment.char = "#", sep='\t'), !V1 %in% c("#id","id","sampleid","sample id","sample-id","#SampleID","#Sample ID", "sample_name", "SampleID","Sample ID"))
  colnames(metadata)<-strsplit(suppressWarnings(readLines(file)[1]), "\t")[[1]]
  colnames(metadata)[1]<-"SampleID"
  return(metadata)
}
#' plot the provenance of a QIIME2 artifact (.qza)
#'
#' extracts embedded data and object metadata into an R session
#'
#' @param artifact the object returned by read_qza
#' @return a text summary for the provenance information
#' @examples \dontrun{print_provenance(artifact)}
#' @export
#'
#'
#'


print_provenance<-function(artifact){
  if(missing(artifact)){stop("Artifact not provided...")}
  
  return(list.tree(artifact$provenance, maxcomp=1000, attr.print=FALSE))
  
}

#' read qiime2 artifacts (.qza)
#'
#' extracts embedded data and object metadata into an R session
#'
#' @param file path to the input file, ex: file="~/data/moving_pictures/table.qza"
#' @param tmp a temporary directory that the object will be decompressed to (default="tempdir()")
#' @param rm should the decompressed object be removed at completion of function (T/F default=TRUE)
#' @return a named list of the following objects:
#' \itemize{
#' \item artifact$data - the raw data ex OTU table as matrix or tree in phylo format
#' \item artifact$uuid - the unique identifer of the artifact
#' \item artifact$type - the semantic type of the object (ex FeatureData[Sequence])
#' \item artifact$format - the format of the qiime artifact
#' \item artifact$provenance - information tracking how the object was created
#' \item artifact$contents - a table of all the files contained within the artifact and their file size
#' \item artifact$version - the reported version for the artifact, a warning error may be thrown if a new version is seen
#' }
#'
#' @examples \dontrun{SVs<-read_qza("data/table.qza")}
#' @export
#'
#'


read_qza <- function(file, tmp, rm) {
  
  if(missing(tmp)){tmp <- tempdir()}
  if(missing(file)){stop("Path to artifact (.qza) not provided")}
  if(missing(rm)){rm=TRUE} #remove the decompressed object from tmp
  
  unzip(file, exdir=tmp)
  unpacked<-unzip(file, exdir=tmp, list=TRUE)
  
  artifact<-read_yaml(paste0(tmp,"/", paste0(gsub("/..+","", unpacked$Name[1]),"/metadata.yaml"))) #start by loading in the metadata not assuming it will be first file listed
  artifact$contents<-data.frame(files=unpacked)
  artifact$contents$size=sapply(paste0(tmp, "/", artifact$contents$files), file.size)
  artifact$version=read.table(paste0(tmp,"/",artifact$uuid, "/VERSION"))
  
  #if(sum(artifact$version$V2==c("2","4","2018.4.0"))!=3){warning("Artifact was not generated with Qiime2 2018.4, if data is not successfully imported, please report here github.com/jbisanz/qiime2R/issues")}#check version and throw warning if new format
  
  #get data dependent on format
  if(grepl("BIOMV", artifact$format)){
    suppressWarnings(artifact$data<-as(biom_data(read_biom(paste0(tmp, "/", artifact$uui,"/data/feature-table.biom"))),"matrix")) #suppressing warning about \n characters
  } else if (artifact$format=="NewickDirectoryFormat"){
    artifact$data<-read.tree(paste0(tmp,"/",artifact$uuid,"/data/tree.nwk"))
  } else if (artifact$format=="DistanceMatrixDirectoryFormat") {
    artifact$data<-as.dist(read.table(paste0(tmp,"/", artifact$uuid, "/data/distance-matrix.tsv"), header=TRUE, row.names=1))
  } else if (grepl("StatsDirFmt", artifact$format)) {
    if(paste0(artifact$uuid, "/data/stats.csv") %in% artifact$contents$files.Name){artifact$data<-read.csv(paste0(tmp,"/", artifact$uuid, "/data/stats.csv"), header=TRUE, row.names=1)}
    if(paste0(artifact$uuid, "/data/stats.tsv") %in% artifact$contents$files.Name){artifact$data<-read.table(paste0(tmp,"/", artifact$uuid, "/data/stats.tsv"), header=TRUE, row.names=1, sep='\t')} #can be tsv or csv
  } else if (artifact$format=="TSVTaxonomyDirectoryFormat"){
    artifact$data<-read.table(paste0(tmp,"/", artifact$uuid, "/data/taxonomy.tsv"), sep='\t', header=TRUE, quote="")
  } else if (artifact$format=="OrdinationDirectoryFormat"){
    
    linesplit<-suppressWarnings(readLines(paste0(tmp,"/", artifact$uuid, "/data/ordination.txt")))
    linesplit<-linesplit[sapply(linesplit, function(x) x!="")]
    
    for (i in 1:length(linesplit)){
      if(grepl("^Eigvals\\t|^Proportion explained\\t|^Species\\t|^Site\\t|^Biplot\\t|^Site constraints\\t", linesplit[i])){
        curfile=strsplit(linesplit[i],"\t")[[1]][1]
      } else {
        write(linesplit[i], paste0(tmp,"/", artifact$uuid, "/data/",curfile,".tmp"), append=TRUE)
      }
    }
    
    for (outs in list.files(paste0(tmp,"/", artifact$uuid,"/data"), full.names = TRUE, pattern = "\\.tmp")){
      NewLab<-gsub(" ", "", toTitleCase(gsub("\\.tmp", "", basename(outs))))
      artifact$data[[NewLab]]<-read.table(outs,sep='\t', header=FALSE)
      if(NewLab %in% c("Eigvals","ProportionExplained")){colnames(artifact$data[[NewLab]])<-paste0("PC",1:ncol(artifact$data[[NewLab]]))}
      if(NewLab %in% c("Site","SiteConstraints")){colnames(artifact$data[[NewLab]])<-c("SampleID", paste0("PC",1:(ncol(artifact$data[[NewLab]])-1)))}
      if(NewLab %in% c("Species")){colnames(artifact$data[[NewLab]])<-c("FeatureID", paste0("PC",1:(ncol(artifact$data[[NewLab]])-1)))}
    }
    
    artifact$data$Vectors<-artifact$data$Site #Rename Site to Vectors so this matches up with the syntax used in the tutorials
    artifact$data$Site<-NULL
    
  } else if (artifact$format=="DNASequencesDirectoryFormat") {
    artifact$data<-readDNAStringSet(paste0(tmp,"/",artifact$uuid,"/data/dna-sequences.fasta"))
  } else if (artifact$format=="AlignedDNASequencesDirectoryFormat") {
    artifact$data<-readDNAMultipleAlignment(paste0(tmp,"/",artifact$uuid,"/data/aligned-dna-sequences.fasta"))
  } else if (grepl("EMPPairedEndDirFmt|EMPSingleEndDirFmt|FastqGzFormat|MultiplexedPairedEndBarcodeInSequenceDirFmt|MultiplexedSingleEndBarcodeInSequenceDirFmt|PairedDNASequencesDirectoryFormat|SingleLanePerSamplePairedEndFastqDirFmt|SingleLanePerSampleSingleEndFastqDirFmt", artifact$format)) {
    artifact$data<-data.frame(files=list.files(paste0(tmp,"/", artifact$uuid,"/data")))
    artifact$data$size<-format(sapply(artifact$data$files, function(x){file.size(paste0(tmp,"/",artifact$uuid,"/data/",x))}, simplify = TRUE))
  } else if (artifact$format=="AlphaDiversityDirectoryFormat") {
    artifact$data<-read.table(paste0(tmp, "/", artifact$uuid, "/data/alpha-diversity.tsv"))
  } else {
    message("Format not supported, only a list of internal files and provenance is being imported.")
    artifact$data<-list.files(paste0(tmp,"/",artifact$uuid, "/data"))
  }
  
  pfiles<-paste0(tmp,"/", grep("..+provenance/..+action.yaml", unpacked$Name, value=TRUE))
  artifact$provenance<-lapply(pfiles, read_yaml)
  names(artifact$provenance)<-grep("..+provenance/..+action.yaml", unpacked$Name, value=TRUE)
  if(rm==TRUE){unlink(paste0(tmp,"/", artifact$uuid), recursive=TRUE)}
  return(artifact)
}
#' generates a phyloseq object from .qza artifacts
#'
#' Construct a phyloseq object from multiple qiime2 artifacts (.qza). Embedded metadata for provenance is not maintained in this function and instead read_qza() should be used.
#'
#' @param features file path for artifact containing a feature (OTU/SV) table
#' @param tree file path for  artifact containing a tree
#' @param taxonomy file path for artifact containg taxonomy
#' @param metadata file path for a qiime2-compliant TSV metadata file
#' @param tmp a temporary directory that the object will be decompressed to (default="/tmp")
#' @return a phyloseq object
#'
#' @examples \dontrun{physeq<-qza_to_phyloseq(features="data/table.qza", tree="data/rooted-tree.qza", taxonomy="data/taxonomy.qza", metdata="data/sample-metadata.qza")}
#' @export
#'
#'
#'

qza_to_phyloseq<-function(features,tree,taxonomy,metadata, tmp){
  
  if(missing(features) & missing(tree) & missing(taxonomy) & missing(metadata)){
    stop("At least one required artifact is needed (features/tree/taxonomy/) or the metadata.")
  }
  
  if(missing(tmp)){tmp <- tempdir()}
  
  
  
  
  argstring<-""
  
  if(!missing(features)){
    features<-read_qza(features, tmp=tmp)$data
    argstring<-paste(argstring, "otu_table(features, taxa_are_rows=T),")
  }
  
  if(!missing(taxonomy)){
    taxonomy<-read_qza(taxonomy, tmp=tmp)$data
    taxt<-strsplit(as.character(taxonomy$Taxon),"\\; ")
    taxt<-lapply(taxt, function(x){length(x)=7;return(x)})
    taxt<-do.call(rbind, taxt)
    taxt<-apply(taxt,2, function(x) replace(x, grepl("^[kpcofgs]__$", x), NA))
    rownames(taxt)<-taxonomy$Feature.ID
    colnames(taxt)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    argstring<-paste(argstring, "tax_table(taxt),")
  }
  
  if(!missing(tree)){
    tree<-read_qza(tree, tmp=tmp)$data
    argstring<-paste(argstring, "phy_tree(tree),")
  }
  
  if(!missing(metadata)){
    
    defline<-suppressWarnings(readLines(metadata)[2])
    if(grepl("^#q2:types", defline)){
      metadata<-read_q2metadata(metadata)
      rownames(metadata)<-metadata$SampleID
      metadata$SampleID<-NULL
    } else{
      metadata<-read.table(metadata, row.names=1, sep='\t', quote="", header=TRUE)
    }
    argstring<-paste(argstring, "sample_data(metadata),")
    sample_data(metadata)
  }
  
  argstring<-gsub(",$","", argstring) #remove trailing ","
  
  physeq<-eval(parse(text=paste0("phyloseq(",argstring,")")))
  
  return(physeq)
}



qza_to_phyloseq_ITS<-function(features,tree,taxonomy,metadata, tmp){
  
  if(missing(features) & missing(tree) & missing(taxonomy) & missing(metadata)){
    stop("At least one required artifact is needed (features/tree/taxonomy/) or the metadata.")
  }
  
  if(missing(tmp)){tmp <- tempdir()}
  
  
  
  
  argstring<-""
  
  if(!missing(features)){
    features<-read_qza(features, tmp=tmp)$data
    argstring<-paste(argstring, "otu_table(features, taxa_are_rows=T),")
  }
  
  if(!missing(taxonomy)){
    taxonomy<-read_qza(taxonomy, tmp=tmp)$data
    taxt<-strsplit(as.character(taxonomy$Taxon),"\\;")
    taxt<-lapply(taxt, function(x){length(x)=7;return(x)})
    taxt<-do.call(rbind, taxt)
    taxt<-apply(taxt,2, function(x) replace(x, grepl("^[kpcofgs]__$", x), NA))
    rownames(taxt)<-taxonomy$Feature.ID
    colnames(taxt)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
    argstring<-paste(argstring, "tax_table(taxt),")
  }
  
  if(!missing(tree)){
    tree<-read_qza(tree, tmp=tmp)$data
    argstring<-paste(argstring, "phy_tree(tree),")
  }
  
  if(!missing(metadata)){
    
    defline<-suppressWarnings(readLines(metadata)[2])
    if(grepl("^#q2:types", defline)){
      metadata<-read_q2metadata(metadata)
      rownames(metadata)<-metadata$SampleID
      metadata$SampleID<-NULL
    } else{
      metadata<-read.table(metadata, row.names=1, sep='\t', quote="", header=TRUE)
    }
    argstring<-paste(argstring, "sample_data(metadata),")
    sample_data(metadata)
  }
  
  argstring<-gsub(",$","", argstring) #remove trailing ","
  
  physeq<-eval(parse(text=paste0("phyloseq(",argstring,")")))
  
  return(physeq)
}