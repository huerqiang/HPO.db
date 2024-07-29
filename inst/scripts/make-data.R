setwd("E:\\enrichplot_export\\DOSE数据更新\\HPO.db数据更新_20240725")
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
obo <- ontologyIndex::get_ontology("E:\\enrichplot_export\\DOSE数据更新\\HPO.db数据更新_20240725\\hp.obo", extract_tags = "everything")
# HPOTERM
HPOTERM <- data.frame(hpid = names(obo$name), term = obo$name)
## 筛选掉is_obsolete
not_obsolete <- names(obo$obsolete)[obo$obsolete == FALSE] |> intersect(HPOTERM$hpid)
# just keep HPO:
not_obsolete <- grep("^HP:", not_obsolete, value = TRUE)
HPOTERM <- HPOTERM[HPOTERM[, 1] %in% not_obsolete, ]
colnames(HPOTERM) <- c("hpid", "term")
dbWriteTable(conn = db, "hpo_term", HPOTERM, row.names=FALSE, overwrite = TRUE)


# ALIAS 
ALIAS <- stack(obo$alt_id)[, c(2, 1)]
colnames(ALIAS) <- c("hpoid", "alias")
ALIAS <- ALIAS[ALIAS[, 1] %in% not_obsolete, ]
dbWriteTable(conn = db, "hpo_alias", ALIAS, row.names=FALSE, overwrite = TRUE)
# SYNONYM
SYNONYM <- stack(obo$synonym)[, c(2, 1)]
colnames(SYNONYM) <- c("hpoid", "synonym")
SYNONYM <- SYNONYM[SYNONYM[, 1] %in% not_obsolete, ]
dbWriteTable(conn = db, "hpo_synonym", SYNONYM, row.names=FALSE, overwrite = TRUE)


# HPOPARENTS
HPOPARENTS <- stack(obo$parents)[, c(2, 1)]
colnames(HPOPARENTS) <- c("hpoid", "parent")
HPOPARENTS <- HPOPARENTS[HPOPARENTS[, 1] %in% not_obsolete, ]
dbWriteTable(conn = db, "hpo_parent", HPOPARENTS, row.names=FALSE, overwrite = TRUE)


# HPOCHILDREN
HPOCHILDREN <- stack(obo$children)[, c(2, 1)]
colnames(HPOCHILDREN) <- c("hpoid", "children")
HPOCHILDREN <- HPOCHILDREN[HPOCHILDREN[, 1] %in% not_obsolete, ]
dbWriteTable(conn = db, "hpo_children", HPOCHILDREN, row.names=FALSE, overwrite = TRUE)
# HPOANCESTOR
HPOANCESTOR <- stack(obo$ancestors)[, c(2, 1)]
HPOANCESTOR <- HPOANCESTOR[HPOANCESTOR[, 1] != HPOANCESTOR[, 2], ]
colnames(HPOANCESTOR) <- c("hpoid", "ancestor")
HPOANCESTOR <- HPOANCESTOR[HPOANCESTOR[, 1] %in% not_obsolete, ]
dbWriteTable(conn = db, "hpo_ancestor", HPOANCESTOR, row.names=FALSE, overwrite = TRUE)
# HPOOFFSPRING
HPOOFFSPRING <- HPOANCESTOR[, c(2, 1)]
colnames(HPOOFFSPRING) <- c("hpoid", "offspring")
HPOOFFSPRING <- HPOOFFSPRING[HPOOFFSPRING[, 1] %in% not_obsolete, ]
dbWriteTable(conn = db, "hpo_offspring", HPOOFFSPRING, row.names=FALSE, overwrite = TRUE)


## HPO2gene
# download from https://hpo.jax.org/data/annotations
file1 <- read.table("E:\\enrichplot_export\\DOSE数据更新\\HPO.db数据更新_20240725\\genes_to_phenotype.txt", 
    sep = "\t", header = TRUE)
file2 <- read.table("E:\\enrichplot_export\\DOSE数据更新\\HPO.db数据更新_20240725\\phenotype_to_genes.txt", 
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
dbWriteTable(conn = db, "hpo_gene", HPOGENE, row.names=FALSE, overwrite = TRUE)

## HPO2DO
# download from https://github.com/mapping-commons/mh_mapping_initiative
setwd("E:\\enrichplot_export\\DOSE数据更新\\MPO.db数据更新_20240725\\mh_mapping_initiative-master")
library(data.table)
file1 <- fread("mappings\\hp_doid_pistoia.sssom.tsv", sep = "\t")
class(file1) <- "data.frame"
file1 <- file1[, c(1, 3)]
colnames(file1) <- c("HP", "DOID")
# HP2OMIM 迂回到DO
phenotype <- fread("E:\\enrichplot_export\\DOSE数据更新\\MPO.db数据更新_20240725\\phenotype.hpoa", sep = "\t", header = TRUE)
class(phenotype) <- "data.frame"
phenotype <- phenotype[, c("hpo_id", "database_id")]
colnames(phenotype) <- c("HP", "OMIM")

# OMIM2DO
# xrefs <- fread("E:\\南方医科大学\\DOSE新paper\\HumanDiseaseOntology-main\\src\\DOreports\\xrefs_in_DO.tsv", sep = "\t", header = TRUE)
xrefs <- fread("E:\\enrichplot_export\\DOSE数据更新\\HPO.db数据更新_20240725\\HumanDiseaseOntology-main\\DOreports\\allXREFinDO.tsv", sep = "\t", header = TRUE)
class(xrefs) <- "data.frame"
# xrefs <- unique(xrefs[grep("OMIM", xrefs[, 3]), c(1, 3)])
xrefs <- unique(xrefs[grep("MIM", xrefs[, 3]), c(1, 3)])
colnames(xrefs) <- c("DOID", "OMIM")
##
xrefs[, 2] <- gsub("^MIM", "OMIM", xrefs[, 2]) 
##

library(dplyr)
HP2DO <- inner_join(phenotype, xrefs, "OMIM")
HP2DO <- unique(HP2DO[, c(1, 3)])
HP2DO <- unique(rbind(file1, HP2DO))
HPODO <- na.omit(HP2DO)
colnames(HPODO) <- c("hpoid", "doid")
dbWriteTable(conn = db, "hpo_do", HPODO, row.names=FALSE, overwrite = TRUE)

## HPO2HPO
mp_hp_eye_impc <- fread("mappings\\mp_hp_eye_impc.sssom.tsv", sep = "\t", header = TRUE)
class(mp_hp_eye_impc) <- "data.frame"
mp_hp_eye_impc <- mp_hp_eye_impc[, c(1,4)]

mp_hp_hwt_impc <- fread("mappings\\mp_hp_hwt_impc.sssom.tsv", sep = "\t", header = TRUE)
class(mp_hp_hwt_impc) <- "data.frame"
mp_hp_hwt_impc <- mp_hp_hwt_impc[, c(1,4)]

mp_hp_mgi_all <- read.table("mappings\\mp_hp_mgi_all.sssom.tsv", sep = "\t", header = TRUE, fill = TRUE, quote = "")
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
dbWriteTable(conn = db, "hpo_mpo", HPOMPO, row.names=FALSE, overwrite = TRUE)

metadata <-rbind(c("DBSCHEMA","HPO_DB"),
        c("DBSCHEMAVERSION","2.0"),
        c("HPOSOURCENAME","Human Phenotype Ontology"),
        c("HPOSOURCURL","https://github.com/DiseaseOntology/HumanDiseaseOntology/blob/main/src/ontology/HumanDO.obo"),
        c("HPOSOURCEDATE","20240725"),
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


map.metadata <- rbind(c("TERM", "Human Phenotype Ontology", "https://hpo.jax.org/app/about","20240725"),
            c("CHILDREN", "Human Phenotype Ontology", "https://hpo.jax.org/app/about","20240725"),
            c("PARENTS", "Human Phenotype Ontology", "https://hpo.jax.org/app/about","20240725"),
            c("ANCESTOR", "Human Phenotype Ontology", "https://hpo.jax.org/app/about","20240725"),
            c("OFFSPRING", "Human Phenotype Ontology", "https://hpo.jax.org/app/about","20240725"),
            c("GENE", "Human Phenotype Ontology", "https://hpo.jax.org/app/about","20240725"),
            c("DO", "Mouse-Human Ontology Mapping Initiative (MHMI)", 
                "https://github.com/mapping-commons/mh_mapping_initiative",
                "20240725"),
            c("MPO", "Mouse-Human Ontology Mapping Initiative (MHMI)", 
                "https://github.com/mapping-commons/mh_mapping_initiative",
                "20240725")
            )	
map.metadata <- as.data.frame(map.metadata)
colnames(map.metadata) <- c("map_name","source_name","source_url","source_date")
dbWriteTable(conn = db, "map_metadata", map.metadata, row.names=FALSE, overwrite = TRUE)


dbListTables(db)
dbListFields(conn = db, "map_metadata")
dbReadTable(conn = db,"map_metadata")
dbDisconnect(db)

