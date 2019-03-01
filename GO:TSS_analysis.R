#PART 1 
#Part A

library(RMySQL)

#connect to MySQL
con <- dbConnect(MySQL(),
                 user = 'bhm263',
                 password = 'FJKVC9su5LwFLMaK',
                 host = 'orion.bio.nyu.edu',
                 dbname='bhm263')


#Create association table
dbSendQuery(con, "DROP TABLE IF EXISTS `association`;")

dbSendQuery(con, "CREATE TABLE `association` (
            `id` int(11) NOT NULL AUTO_INCREMENT,
            `term_id` int(11) NOT NULL,
            `gene_product_id` int(11) NOT NULL,
            `is_not` int(11) DEFAULT NULL,
            `role_group` int(11) DEFAULT NULL,
            `assocdate` int(11) DEFAULT NULL,
            `source_db_id` int(11) DEFAULT NULL,
            PRIMARY KEY (`id`),
            UNIQUE KEY `a0` (`id`),
            KEY `source_db_id` (`source_db_id`),
            KEY `a1` (`term_id`),
            KEY `a2` (`gene_product_id`),
            KEY `a3` (`term_id`,`gene_product_id`),
            KEY `a4` (`id`,`term_id`,`gene_product_id`),
            KEY `a5` (`id`,`gene_product_id`),
            KEY `a6` (`is_not`,`term_id`,`gene_product_id`),
            KEY `a7` (`assocdate`)
) ENGINE=MyISAM AUTO_INCREMENT=258890627 DEFAULT CHARSET=latin1;")

#upload association.txt file data into table
dbSendQuery(con, " LOAD DATA LOCAL INFILE 
            '/home/bhm263/Midterm/go_monthly-assocdb-tables/association.txt' 
            INTO TABLE association
            LINES TERMINATED BY '\n';")

#Create gene_procuct table
dbSendQuery(con, "DROP TABLE IF EXISTS `gene_product`;")

dbSendQuery(con, "CREATE TABLE `gene_product` (
            `id` int(11) NOT NULL AUTO_INCREMENT,
            `symbol` varchar(128) NOT NULL,
            `dbxref_id` int(11) NOT NULL,
            `species_id` int(11) DEFAULT NULL,
            `type_id` int(11) DEFAULT NULL,
            `full_name` text,
            PRIMARY KEY (`id`),
            UNIQUE KEY `dbxref_id` (`dbxref_id`),
            UNIQUE KEY `g0` (`id`),
            KEY `type_id` (`type_id`),
            KEY `g1` (`symbol`),
            KEY `g2` (`dbxref_id`),
            KEY `g3` (`species_id`),
            KEY `g4` (`id`,`species_id`),
            KEY `g5` (`dbxref_id`,`species_id`),
            KEY `g6` (`id`,`dbxref_id`),
            KEY `g7` (`id`,`species_id`),
            KEY `g8` (`id`,`dbxref_id`,`species_id`)
) ENGINE=MyISAM AUTO_INCREMENT=43505133 DEFAULT CHARSET=latin1;")

#upload txt file data into table
dbSendQuery(con, " LOAD DATA LOCAL INFILE 
            '/home/bhm263/Midterm/go_monthly-assocdb-tables/gene_product.txt' 
            INTO TABLE gene_product
            LINES TERMINATED BY '\n';")

#Create species table
dbSendQuery(con, "DROP TABLE IF EXISTS `species`;")

dbSendQuery(con, "CREATE TABLE `species` (
            `id` int(11) NOT NULL AUTO_INCREMENT,
            `ncbi_taxa_id` int(11) DEFAULT NULL,
            `common_name` varchar(255) DEFAULT NULL,
            `lineage_string` text,
            `genus` varchar(55) DEFAULT NULL,
            `species` varchar(255) DEFAULT NULL,
            `parent_id` int(11) DEFAULT NULL,
            `left_value` int(11) DEFAULT NULL,
            `right_value` int(11) DEFAULT NULL,
            `taxonomic_rank` varchar(255) DEFAULT NULL,
            PRIMARY KEY (`id`),
            UNIQUE KEY `sp0` (`id`),
            UNIQUE KEY `ncbi_taxa_id` (`ncbi_taxa_id`),
            KEY `sp1` (`ncbi_taxa_id`),
            KEY `sp2` (`common_name`),
            KEY `sp3` (`genus`),
            KEY `sp4` (`species`),
            KEY `sp5` (`genus`,`species`),
            KEY `sp6` (`id`,`ncbi_taxa_id`),
            KEY `sp7` (`id`,`ncbi_taxa_id`,`genus`,`species`),
            KEY `sp8` (`parent_id`),
            KEY `sp9` (`left_value`),
            KEY `sp10` (`right_value`),
            KEY `sp11` (`left_value`,`right_value`),
            KEY `sp12` (`id`,`left_value`),
            KEY `sp13` (`genus`,`left_value`,`right_value`)
) ENGINE=MyISAM AUTO_INCREMENT=1547882 DEFAULT CHARSET=latin1;")

#upload txt file data into table
dbSendQuery(con, " LOAD DATA LOCAL INFILE 
            '/home/bhm263/Midterm/go_monthly-assocdb-tables/species.txt' 
            INTO TABLE species
            LINES TERMINATED BY '\n';")

#Create term table
dbSendQuery(con, "DROP TABLE IF EXISTS `term`;")

dbSendQuery(con, "CREATE TABLE `term` (
            `id` int(11) NOT NULL AUTO_INCREMENT,
            `name` varchar(255) NOT NULL DEFAULT '',
            `term_type` varchar(55) NOT NULL,
            `acc` varchar(255) NOT NULL,
            `is_obsolete` int(11) NOT NULL DEFAULT '0',
            `is_root` int(11) NOT NULL DEFAULT '0',
            `is_relation` int(11) NOT NULL DEFAULT '0',
            PRIMARY KEY (`id`),
            UNIQUE KEY `acc` (`acc`),
            UNIQUE KEY `t0` (`id`),
            KEY `t1` (`name`),
            KEY `t2` (`term_type`),
            KEY `t3` (`acc`),
            KEY `t4` (`id`,`acc`),
            KEY `t5` (`id`,`name`),
            KEY `t6` (`id`,`term_type`),
            KEY `t7` (`id`,`acc`,`name`,`term_type`)
) ENGINE=MyISAM AUTO_INCREMENT=46049 DEFAULT CHARSET=latin1;")

#upload txt file data into table
dbSendQuery(con, " LOAD DATA LOCAL INFILE 
            '/home/bhm263/Midterm/go_monthly-assocdb-tables/term.txt' 
            INTO TABLE term
            FIELDS TERMINATED BY '\t'
            LINES TERMINATED BY '\n';")

#Create term2term table
dbSendQuery(con, "DROP TABLE IF EXISTS `term2term`;")

dbSendQuery(con, "CREATE TABLE `term2term` (
            `id` int(11) NOT NULL AUTO_INCREMENT,
            `relationship_type_id` int(11) NOT NULL,
            `term1_id` int(11) NOT NULL,
            `term2_id` int(11) NOT NULL,
            `complete` int(11) NOT NULL DEFAULT '0',
            PRIMARY KEY (`id`),
            UNIQUE KEY `term1_id` (`term1_id`,`term2_id`,`relationship_type_id`),
            KEY `tt1` (`term1_id`),
            KEY `tt2` (`term2_id`),
            KEY `tt3` (`term1_id`,`term2_id`),
            KEY `tt4` (`relationship_type_id`)
) ENGINE=MyISAM AUTO_INCREMENT=93714 DEFAULT CHARSET=latin1;")

#upload txt file data into table
dbSendQuery(con, " LOAD DATA LOCAL INFILE 
            '/home/bhm263/Midterm/go_monthly-assocdb-tables/term2term.txt' 
            INTO TABLE term2term
            LINES TERMINATED BY '\n'
            IGNORE 40 ROWS;")

#Part B

##1

#function that puts variable into proper format for query
addVals <- function(x) paste0("'", x, "'",  collapse = ",")

#INPUT 'genus', 'species'
#finds the gene_product.id from the species_id (use this id to obtain gene_product.symbol)
#matches the gene_product.id with the association.gene_product_id to provide the association.term_id
#matches association.term_id with term.id to provide term.acc or the GO term
#OUTPUT a file with two columns 1 go term 2 gene symbol

gene_assoc_by_species <- function(string1, string2) {
  queryform1 <- "SELECT id FROM species WHERE genus= %s AND species= %s;"
  query1 <- sprintf(queryform1, addVals(string1), addVals(string2))
  z<- dbGetQuery(con, query1)
  queryform2<- "SELECT term.acc from term 
              INNER JOIN association ON term.id= association.term_id
              WHERE term.id IN (SELECT term_id FROM association
              WHERE gene_product_id IN 
              (SELECT id FROM gene_product WHERE species_id= (%s)))
              ORDER BY association.id;"
  query2<- sprintf(queryform2, addVals(z))
  y<- dbGetQuery(con, query2)
  queryform3<- "SELECT gene_product.symbol from gene_product 
              INNER JOIN association ON gene_product.id= association.gene_product_id
              WHERE gene_product_id IN 
              (SELECT id FROM gene_product WHERE species_id= (%s))
              ORDER BY association.id;"
  query3<- sprintf(queryform3, addVals(z))
  x<- dbGetQuery(con, query3)
  file1<- data.frame(y, x)
  write.table(file1, file='B1_bhm263.txt)', quote=FALSE, sep='\t', col.names = NA)
  }

gene_assoc_by_species('Homo', 'sapiens')

##2

Z<- dbGetQuery(con, "SELECT term.acc FROM term
               INNER JOIN term2term ON term.id=term2term.term1_id
               WHERE term.id IN (SELECT term1_id FROM term2term)
               GROUP BY term2term.id; ")


X<- dbGetQuery(con, "SELECT term.acc FROM term
               INNER JOIN term2term ON term.id=term2term.term2_id
               WHERE term.id IN (SELECT term2_id FROM term2term)
               GROUP BY term2term.id; ")

parent2child<- data.frame(Z, X)

write.table(parent2child, file='B2_bhm263.txt)', quote=FALSE, sep='\t', col.names = NA)

##3
#function that puts a value into format for SQL query
addVals <- function(x) paste0("'", x, "'",  collapse = ",")

#INPUT go term 'GO:0000000'
#finds the term.id from term.acc (go term inputed)
#compares term.id to term2term.term2_id to find any matches to child 'node'
#returns corresponding term2term.term1_id to provide parent node
#matches term2term.term1_id to term.id to obtain term.acc or Go term
#Output: parent go-terms
go_term2parent <- function(string) {
  queryform8 <- "SELECT acc from term WHERE id IN 
                (SELECT term1_id from term2term WHERE term2_id IN
               (SELECT id FROM term where acc = %s)) ;"
  query8 <- sprintf(queryform8, addVals(string))
  z<- dbGetQuery(con, query8)}

r<- go_term2parent('GO:0000256')
r

#PART 2

set.seed(4)
library(Repitools)
library(GenomicRanges)

load("~/Midterm/hg38_gene.anno.RData")

genelist<- BAM2GRangesList(c('GSM409312_UCSD_H3K36me3.bam', 'GSM409308_UCSD_H3K4me3.bam', 'GSM409307_UCSD_H3K4me1.bam'))

##1

Scores <- featureScores(genelist, gene.anno,
                      up = 5000, down = 2000, freq = 1000, s.width = 500)

scores <- featureScores(genelist, gene.anno, up = 5000, down = 2000, dist = "base",
                      freq = 200, s.width = 500)
cluster_plot <- clusterPlots(scores, scale = function(x) sqrt(x), plot.type = "line",
                   t.name = "Data", n.clusters = 1 )

##2
kmean_cluster<- clusterPlots(scores, scale = function(x) sqrt(x), plot.type = "line",
                          t.name = "Data", n.clusters = 4)

##3
#The H3K4me3 modification is localized around the TSS.
#This histone modification appears concentrated to a particular gene subset, not the entire genome

##4
cp <- clusterPlots(scores, scale = function(x) sqrt(x), plot.type = "heatmap",
                   t.name = "Data", n.clusters = 4)
#The heatmap and line plot are consistent

##5

#filter out for human genes only using gene_product.species_id
id<- dbGetQuery(con, "SELECT id from species WHERE genus= 'Homo' and species= 'sapiens';")


#obtain gene_product.id
#match gene_product.id to association.gene_product_id to obtain association.term_id
#match association.term_id to term.id and filter for term_type contains transcription factor activity or binding
r<-dbGetQuery(con, "SELECT distinct(term.id) from term 
              INNER JOIN association ON term.id= association.term_id
              WHERE term.id IN (SELECT term_id FROM association
              WHERE gene_product_id IN 
              (SELECT id FROM gene_product WHERE species_id= '462613'))
              AND term.name LIKE '%transcription factor binding%'
              OR term.name LIKE '%transcription factor activity%' 
              ORDER BY association.id;")

#function that puts variable into proper format for query
addVals <- function(x) paste0("'", x, "'",  collapse = ",")

#matches selected term_ids(variable left open for them) to gene_product ids and provides symbols
queryform3<- "SELECT symbol from gene_product
              WHERE gene_product.id IN
              (SELECT gene_product_id from association 
              WHERE term_id IN (%s)
              AND is_not= '0')
              GROUP BY symbol;"

#converts sql data into a vector for addValues function
x<- c()
for (i in 1:length(r$id)) {
  x[i]<- r$id[i]}

#function to put variable into format
addVals <- function(x) paste0("'", x, "'",  collapse = ",")

#puts function and formatted variables together
query8 <- sprintf(queryform3, addVals(x))
           
#resulting symbol list
Z <- dbGetQuery(con, query8)

#idea is to match a symbol from Z(gene list) with a symbol in gene.anno table
which.loss <- which(gene.anno$symbol %in% as.character(Z$symbol))

#create enrichment plots
#There appears to be enrichment on the second and third histone plots
profilePlots(Scores, gene.lists = list("Genes of Interest"= which.loss) ,cols = "red")


##6

#download knownGene table select for 5'UTR
#inlcude snapshots etc
#include only 5’UTRS that are longer than 200 bases 


##7

#read bedfile as dataframe and only take the first 4000 5'UTRs
dat = read.table(file="known_gene.txt", header=F, stringsAsFactors=F)

sub_dat<- subset(dat, dat$V3 - dat$V2 > 200)

utr5 <- sub_dat[1:4000,]

# write data.frame
write.table(utr5,file="hg38_utr5_out.bed",row.names=F,col.names=F,quote = F,sep="\t")

library(metagene)
region <- c('hg38_utr5_out.bed')

bam_file <- c('GSM409307_UCSD_H3K4me1.bam', 'GSM409308_UCSD_H3K4me3.bam', 'GSM409312_UCSD_H3K36me3.bam')
mg <- metagene$new(regions = region, bam_files = bam_file)
dev.off()
df <- mg$plot()

