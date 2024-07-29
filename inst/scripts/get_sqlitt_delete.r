
setwd("E:\\enrichplot_export\\HPO.db")
packagedir <- getwd()
sqlite_path <- paste(packagedir, sep=.Platform$file.sep, "inst", "extdata")
if(!dir.exists(sqlite_path)){dir.create(sqlite_path,recursive = TRUE)}
dbfile <- file.path(sqlite_path, "HPO.sqlite")
unlink(dbfile)

###################################################
### create database
###################################################
## Create the database file
library(RSQLite)
drv <- dbDriver("SQLite")
db <- dbConnect(drv, dbname=dbfile)
## dbDisconnect(db)

source("E:\\enrichplot_export\\HPO.db\\inst\\extdata\\parse-obo.R")
library(dplyr)
## DOTERM
source("E:\\enrichplot_export\\HPO.db\\inst\\extdata\\parse-obo.R")
doobo <- parse_do("E:\\南方医科大学\\DOSE新paper\\hp.obo")
doterm <- doobo$doinfo[, c(1,2)]
colnames(doterm) <- c("hpoid", "term")
HPOTERM <- doterm
dbWriteTable(conn = db, "do_term", HPOTERM, row.names=FALSE, overwrite = TRUE)


## ALIAS
alias <- doobo$alias
colnames(alias) <- c("hpoid", "alias")
ALIAS <- alias
dbWriteTable(conn = db, "do_alias", ALIAS, row.names=FALSE, overwrite = TRUE)

## SYNONYM
synonym <- doobo$synonym
colnames(synonym) <- c("hpoid", "synonym")
SYNONYM <- synonym
dbWriteTable(conn = db, "do_synonym", SYNONYM, row.names=FALSE, overwrite = TRUE)

## DOPARENTS
HPOPARENTS <- doobo$rel
colnames(HPOPARENTS) <- c("hpoid", "parent")
dbWriteTable(conn = db, "do_parent", HPOPARENTS, row.names=FALSE)


## DOCHILDREN
HPOCHILDREN <- doobo$rel[, c(2,1)]
HPOCHILDREN <- HPOCHILDREN[order(HPOCHILDREN[, 1]), ]
colnames(HPOCHILDREN) <- c("hpoid", "children")
dbWriteTable(conn = db, "do_children", HPOCHILDREN, row.names=FALSE)

## DOANCESTOR
ancestor_list <- split(HPOPARENTS[, 2], HPOPARENTS[, 1])
getAncestor <- function(id) {
    ans_temp <- which(HPOPARENTS[, 1] %in% ancestor_list[[id]])
    ids <- HPOPARENTS[ans_temp, 2]
    content <- c(ancestor_list[[id]], ids)
    while(!all(is.na(ids))) {
        ans_temp <- which(HPOPARENTS[, 1] %in% ids)
        ids <- HPOPARENTS[ans_temp, 2]
        content <- c(content, ids)
    }
    content[!is.na(content)]
}

for (id in names(ancestor_list)) {
    ancestor_list[[id]] <- getAncestor(id)
}
ancestordf <- stack(ancestor_list)[, c(2, 1)]
ancestordf[, 1] <- as.character(ancestordf[, 1])
ancestordf <- unique(ancestordf)
HPOANCESTOR <- ancestordf
colnames(HPOANCESTOR) <- c("hpoid", "ancestor")
dbWriteTable(conn = db, "do_ancestor", HPOANCESTOR, row.names=FALSE)


# DOOFFSPRING
HPOOFFSPRING <- ancestordf[, c(2, 1)]
HPOOFFSPRING <- HPOOFFSPRING[order(HPOOFFSPRING[, 1]), ]
colnames(HPOOFFSPRING) <- c("hpoid", "offspring")
dbWriteTable(conn = db, "do_offspring", HPOOFFSPRING, row.names=FALSE, overwrite = TRUE)


## HPO2gene
file1 <- read.table("E:\\南方医科大学\\DOSE新paper\\HPO\\genes_to_phenotype.txt", 
    sep = "\t", header = TRUE)
file2 <- read.table("E:\\南方医科大学\\DOSE新paper\\HPO\\phenotype_to_genes.txt", 
    sep = "\t", header = TRUE)
file1 <- file1[, c("hpo_id", "gene_symbol")]
file2 <- file2[, c("hpo_id", "gene_symbol")]
hpo2gene <- unique(rbind(file1, file2))
library(clusterProfiler)
library(org.Hs.eg.db)
gene_bitr <- bitr(geneID = hpo2gene[,2], "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)
hpo2gene$geneId <- gene_bitr[match(hpo2gene$gene_symbol, gene_bitr[, 1]), 2]
hpo2gene <- unique(hpo2gene[, c("hpo_id", "geneId")])

HPOGENE <- na.omit(hpo2gene)
colnames(HPOGENE) <- c("hpoid", "gene")
dbWriteTable(conn = db, "hpoid_gene", HPOGENE, row.names=FALSE, overwrite = TRUE)

## HPO2DO
setwd("E:\\南方医科大学\\DOSE新paper\\mh_mapping_initiative-master\\mh_mapping_initiative-master")
library(data.table)
file1 <- fread("mappings\\hp_hpoid_pistoia.sssom.tsv", sep = "\t")
class(file1) <- "data.frame"
file1 <- file1[, c(1, 3)]
colnames(file1) <- c("HP", "DOID")
# HP2OMIM 迂回到DO
phenotype <- fread("E:\\南方医科大学\\DOSE新paper\\MGI\\phenotype.hpoa", sep = "\t", header = TRUE)
class(phenotype) <- "data.frame"
phenotype <- phenotype[, c(4, 5)]
colnames(phenotype) <- c("HP", "OMIM")

# OMIM2DO
xrefs <- fread("E:\\南方医科大学\\DOSE新paper\\HumanDiseaseOntology-main\\src\\DOreports\\xrefs_in_DO.tsv", sep = "\t", header = TRUE)
class(xrefs) <- "data.frame"
xrefs <- unique(xrefs[grep("OMIM", xrefs[, 3]), c(1, 3)])
colnames(xrefs) <- c("DOID", "OMIM")
library(dplyr)
HP2DO <- inner_join(phenotype, xrefs, "OMIM")
HP2DO <- unique(HP2DO[, c(1, 3)])
HP2DO <- unique(rbind(file1, HP2DO))
HPODO <- na.omit(HP2DO)
colnames(HPODO) <- c("hpoid", "doid")
dbWriteTable(conn = db, "hpoid_doid", HPODO, row.names=FALSE, overwrite = TRUE)

## HPO2HPO
mp_hp_eye_impc <- fread("mappings\\mp_hp_eye_impc.sssom.tsv", sep = "\t", header = TRUE)
class(mp_hp_eye_impc) <- "data.frame"
mp_hp_eye_impc <- mp_hp_eye_impc[, c(1,4)]

mp_hp_hwt_impc <- fread("mappings\\mp_hp_hwt_impc.sssom.tsv", sep = "\t", header = TRUE)
class(mp_hp_hwt_impc) <- "data.frame"
mp_hp_hwt_impc <- mp_hp_hwt_impc[, c(1,4)]

mp_hp_mgi_all <- fread("mappings\\mp_hp_mgi_all.sssom.tsv", sep = "\t", header = TRUE)
class(mp_hp_mgi_all) <- "data.frame"
mp_hp_mgi_all <- mp_hp_mgi_all[, c(5,1)]

mp_hp_owt_impc <- fread("mappings\\mp_hp_owt_impc.sssom.tsv", sep = "\t", header = TRUE)
class(mp_hp_owt_impc) <- "data.frame"
mp_hp_owt_impc <- mp_hp_owt_impc[, c(1,4)]
mp_hp_owt_impc <- mp_hp_owt_impc[grep("HP:", mp_hp_owt_impc[, 2]), ]

mp_hp_pat_impc <- fread("mappings\\mp_hp_pat_impc.sssom.tsv", sep = "\t", header = TRUE)
class(mp_hp_pat_impc) <- "data.frame"
mp_hp_pat_impc <- mp_hp_pat_impc[, c(1,4)]
mp_hp_pat_impc <- mp_hp_pat_impc[grep("HP:", mp_hp_pat_impc[, 2]), ]

mp_hp_pistoia <- fread("mappings\\mp_hp_pistoia.sssom.tsv", sep = "\t", header = TRUE)
class(mp_hp_pistoia) <- "data.frame"
mp_hp_pistoia <- mp_hp_pistoia[, c(1,3)]
mp_hp_pistoia <- mp_hp_pistoia[grep("HP:", mp_hp_pistoia[, 2]), ]

colnames(mp_hp_eye_impc) <- colnames(mp_hp_hwt_impc) <- colnames(mp_hp_mgi_all) <- 
    colnames(mp_hp_owt_impc) <- colnames(mp_hp_pat_impc) <- colnames(mp_hp_pistoia) <- c("MP", "HP")
mp2hp <- do.call(rbind, list(mp_hp_eye_impc, mp_hp_hwt_impc, mp_hp_mgi_all, mp_hp_owt_impc, mp_hp_pat_impc, mp_hp_pistoia))
mp2hp <- mp2hp[grep("HP:", mp2hp[, 2]), ]
mp2hp <- unique(mp2hp)
HPOMPO <- na.omit(mp2hp[, c(2,1)])
colnames(HPOMPO) <- c("hpoid", "mpoid")
dbWriteTable(conn = db, "hpoid_mpoid", HPOMPO, row.names=FALSE, overwrite = TRUE)

metadata <-rbind(c("DBSCHEMA","HPO_DB"),
        c("DBSCHEMAVERSION","1.0"),
        c("HPOSOURCENAME","Human Phenotype Ontology"),
        c("HPOSOURCURL","https://github.com/DiseaseOntology/HumanDiseaseOntology/blob/main/src/ontology/HumanDO.obo"),
        c("HPOSOURCEDATE","20230405"),
        c("Db type", "HPODb"))

metadata <- as.data.frame(metadata)
colnames(metadata) <- c("name", "value") 
dbWriteTable(conn = db, "metadata", metadata, row.names=FALSE, overwrite = TRUE)



map.counts<-rbind(c("TERM", nrow(HPOTERM)),
        # c("OBSOLETE","$obsolete_counts"),
        c("CHILDREN", nrow(HPOCHILDREN)),
        c("PARENTS", nrow(HPOPARENTS)),
        c("ANCESTOR", nrow(HPOANCESTOR)),
        c("OFFSPRING", nrow(HPOOFFSPRING)),
        c("GENE", nrow(HPOGENE)),
        c("DO", nrow(HPODO)),
        c("MPO", nrow(HPOMPO))
        )


map.counts <- as.data.frame(map.counts)
colnames(map.counts) <- c("map_name","count")
# dbWriteTable(conn = db, "map.counts", map.counts, row.names=FALSE)
dbWriteTable(conn = db, "map_counts", map.counts, row.names=FALSE, overwrite = TRUE)

dbListTables(db)
dbListFields(conn = db, "metadata")
dbReadTable(conn = db,"metadata")


map.metadata <- rbind(c("TERM", "Human Phenotype Ontology", "https://hpo.jax.org/app/about","20230405"),
            c("CHILDREN", "Human Phenotype Ontology", "https://hpo.jax.org/app/about","20230405"),
            c("PARENTS", "Human Phenotype Ontology", "https://hpo.jax.org/app/about","20230405"),
            c("ANCESTOR", "Human Phenotype Ontology", "https://hpo.jax.org/app/about","20230405"),
            c("OFFSPRING", "Human Phenotype Ontology", "https://hpo.jax.org/app/about","20230405"),
            c("GENE", "Human Phenotype Ontology", "https://hpo.jax.org/app/about","20230405"),
            c("DO", "Mouse-Human Ontology Mapping Initiative (MHMI)", 
                "https://github.com/mapping-commons/mh_mapping_initiative",
                "20230201"),
            c("MPO", "Mouse-Human Ontology Mapping Initiative (MHMI)", 
                "https://github.com/mapping-commons/mh_mapping_initiative",
                "20230201")
            )	
map.metadata <- as.data.frame(map.metadata)
colnames(map.metadata) <- c("map_name","source_name","source_url","source_date")
dbWriteTable(conn = db, "map_metadata", map.metadata, row.names=FALSE, overwrite = TRUE)


dbListTables(db)
dbListFields(conn = db, "map_metadata")
dbReadTable(conn = db,"map_metadata")
dbDisconnect(db)
