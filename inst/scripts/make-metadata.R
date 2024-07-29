meta <- data.frame(
     Title = "Human Phenotype Ontology database",
     Description = paste0("A set of annotation maps describing the Human Phenotype Ontology, ",
         "it mainly contains four kinds of annotation information: ",
         "Human Phenotype Ontology data,",
         "Human gene - phenotype association data, ",
         "Human phenotype - mouse phenotype association data",
         "Human phenotype - human disease association data"
     ),
     BiocVersion = "3.19",
     Genome = NA,
     SourceType = "Multiple",
     SourceUrl = paste("https://hpo.jax.org/app/data/ontology", 
         "https://github.com/mapping-commons/mh_mapping_initiative",
         sep = ","),
     SourceVersion = "2",
     Species = NA,
     TaxonomyId = NA,
     Coordinate_1_based = TRUE,
     DataProvider = "The Human Phenotype Ontology",
     Maintainer = "Erqiang Hu <13766876214@163.com>",
     RDataClass = "SQLite",
     # DispatchClass = "SQLiteFile",
     DispatchClass = "FilePath",
     RDataPath = "HPO.db/V2/HPO.sqlite",
     ResourceName = "HPO.sqlite",
     Tags = "Annotation"
)
write.csv(meta, file="inst/extdata/metadata_v2.csv", row.names=FALSE)


