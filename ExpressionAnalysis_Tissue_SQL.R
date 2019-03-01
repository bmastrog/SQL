#Part 1  
library(RMySQL) #call mysql package

#connect to databases on orion
con1 <- dbConnect(MySQL(),
                 user = 'ga1009',
                 password = 'mkatari@nyu',
                 host = 'orion.bio.nyu.edu',
                 dbname='hg19')

#isolate cancer gene names
names<- read.delim('cancerGenes.names.txt', header=FALSE) 
names1 <- names[,1]

#creates a query where a variable can be inserted
queryform <- "SELECT name, name2, chrom, strand, txStart,txEnd FROM refGene WHERE name2 IN (%s)"

#function that puts r objects into format that can be inserted into RMySQL query
addVals <- function(x) paste0("'", x, "'", collapse = ",")

#join query and formatted variables
query1 <- sprintf(queryform, addVals(names1))

#execute query
result<- dbGetQuery(con1, query1)

#connect to my database
con2 <- dbConnect(MySQL(),
                  user = 'bhm263',
                  password = 'FJKVC9su5LwFLMaK',
                  host = 'orion.bio.nyu.edu',
                  dbname='bhm263')

#write the results of hg19 query into a table in my database
dbWriteTable(conn = con2, name = 'refGenesCoor', value = as.data.frame(result))

#Part 2

brain_output <- numeric(nrow(result)) #create empty vector to put query results into

#loops through the r object containing refGenesCoor info 
#provides gene hits within a certain range for each gene for this tissue
for(row in seq_len(nrow(result))){
  query2 <- paste0("SELECT count(*) FROM burgeRnaSeqGemMapperAlignBrain	
  WHERE chrom= '", result$chrom[row],"'
  AND chromStart > '", result$txStart[row],"' 
  AND chromEnd < '", result$txEnd[row],"'")
  brain_output[row]<- dbGetQuery(con1, query2)}

breast_output <- numeric(nrow(result))#create empty vector to put query results into

#loops through the r object containing refGenesCoor info 
#provides gene hits within a certain range for each gene for this tissue
for(row in seq_len(nrow(result))){
  query3 <- paste0("SELECT count(*) FROM burgeRnaSeqGemMapperAlignBreast	WHERE chrom= '", result$chrom[row],"'
  AND chromStart > '", result$txStart[row],"' 
  AND chromEnd < '", result$txEnd[row],"'")
  breast_output[row]<- dbGetQuery(con1, query3)}

colon_output <- numeric(nrow(result)) #create empty vector to put query results into

#loops through the r object containing refGenesCoor info 
#provides gene hits within a certain range for each gene for this tissue
for(row in seq_len(nrow(result))){
  query4 <- paste0("SELECT count(*) FROM burgeRnaSeqGemMapperAlignColon	WHERE chrom= '", result$chrom[row],"'
  AND chromStart > '", result$txStart[row],"' 
  AND chromEnd < '", result$txEnd[row],"'")
  colon_output[row]<- dbGetQuery(con1, query4)}

heart_output <- numeric(nrow(result)) #create empty vector to put query results into

#loops through the r object containing refGenesCoor info 
#provides gene hits within a certain range for each gene for this tissue
for(row in seq_len(nrow(result))){
  query5 <- paste0("SELECT count(*) FROM burgeRnaSeqGemMapperAlignHeart	WHERE chrom= '", result$chrom[row],"'
  AND chromStart > '", result$txStart[row],"' 
  AND chromEnd < '", result$txEnd[row],"'")
  heart_output[row]<- dbGetQuery(con1, query5)}

liver_output <- numeric(nrow(result)) #create empty vector to put query results into

#loops through the r object containing refGenesCoor info 
#provides gene hits within a certain range for each gene for this tissue
for(row in seq_len(nrow(result))){
  query6 <- paste0("SELECT count(*) FROM burgeRnaSeqGemMapperAlignLiver	WHERE chrom= '", result$chrom[row],"'
  AND chromStart > '", result$txStart[row],"' 
  AND chromEnd < '", result$txEnd[row],"'")
  liver_output[row]<- dbGetQuery(con1, query6)}

lymphnode_output <- numeric(nrow(result)) #create empty vector to put query results into

#loops through the r object containing refGenesCoor info 
#provides gene hits within a certain range for each gene for this tissue
for(row in seq_len(nrow(result))){
  query7 <- paste0("SELECT count(*) FROM burgeRnaSeqGemMapperAlignLymphNode	WHERE chrom= '", result$chrom[row],"'
  AND chromStart > '", result$txStart[row],"' 
  AND chromEnd < '", result$txEnd[row],"'")
  lymphnode_output[row]<- dbGetQuery(con1, query7)}

skelmuscle_output <- numeric(nrow(result)) #create empty vector to put query results into

#loops through the r object containing refGenesCoor info 
#provides gene hits within a certain range for each gene for this tissue
for(row in seq_len(nrow(result))){
  query8 <- paste0("SELECT count(*) FROM burgeRnaSeqGemMapperAlignSkelMuscle	WHERE chrom= '", result$chrom[row],"'
  AND chromStart > '", result$txStart[row],"' 
  AND chromEnd < '", result$txEnd[row],"'")
  skelmuscle_output[row]<- dbGetQuery(con1, query8)}

testes_output <- numeric(nrow(result)) #create empty vector to put query results into

#loops through the r object containing refGenesCoor info 
#provides gene hits within a certain range for each gene for this tissue
for(row in seq_len(nrow(result))){
  query9 <- paste0("SELECT count(*) FROM burgeRnaSeqGemMapperAlignTestes	WHERE chrom= '", result$chrom[row],"'
  AND chromStart > '", result$txStart[row],"' 
  AND chromEnd < '", result$txEnd[row],"'")
  testes_output[row]<- dbGetQuery(con1, query9)}

#connect to my database
con2 <- dbConnect(MySQL(),
                  user = 'bhm263',
                  password = 'FJKVC9su5LwFLMaK',
                  host = 'orion.bio.nyu.edu',
                  dbname='bhm263')

dbSendQuery(con2, "DROP TABLE IF EXISTS `GeneRNAseq`;")

#create table to input RNA seq data 
dbSendQuery(con2, "
CREATE TABLE GeneRNAseq (
name_id INT PRIMARY KEY,
name VARCHAR(50),
brain FLOAT,
breast FLOAT,
colon FLOAT,
heart FLOAT,
liver FLOAT,
lymphnode FLOAT,
skelmuscle FLOAT,
testes FLOAT);")

name_id<- seq(1, length(result$name2)) #create primary key vector
gene_name <- result$name2 #create gene name vector

#loops through data vectors for gene name and values from each tissue type 
#values go into table created to link gene name and RNA seq data
for(row in seq_len(nrow(result))){
  table <- paste("INSERT INTO GeneRNAseq VALUES(",name_id[row]," , '", gene_name[row],"' , 
  ",brain_output[row]," , ",breast_output[row], " ,
  ",colon_output[row]," , ",heart_output[row]," , 
  ",liver_output[row]," , ",lymphnode_output[row]," , 
  ",skelmuscle_output[row],", ",testes_output[row],")")
  dbGetQuery(con2, table)}


#Part 3

#create empty vectors to enter normalized data into
brain_data <- c()
breast_data <- c()
colon_data <- c() 
heart_data <- c()
liver_data <- c()
lymphnode_data <- c()
skelmuscle_data <- c()
testes_data <- c()

#creates r object from GeneRNAseq table to call on data
genernaseq<- dbGetQuery(con2, "SELECT * FROM GeneRNAseq" )

#loops go through each row for a particular tissue type and 
#do first step of normalization calculation
for(row in seq_len(nrow(result))){
   brain_data[row]<- genernaseq$brain[row]/(result$txEnd[row]- result$txStart[row]+1)} 
   
for(row in seq_len(nrow(result))){
  breast_data[row]<- genernaseq$breast[row]/(result$txEnd[row]- result$txStart[row]+1)} 

for(row in seq_len(nrow(result))){
  colon_data[row]<- genernaseq$colon[row]/(result$txEnd[row]- result$txStart[row]+1)}

for(row in seq_len(nrow(result))){
  heart_data[row]<- genernaseq$heart[row]/(result$txEnd[row]- result$txStart[row]+1)}

for(row in seq_len(nrow(result))){
  liver_data[row]<- genernaseq$liver[row]/(result$txEnd[row]- result$txStart[row]+1)}

for(row in seq_len(nrow(result))){
  lymphnode_data[row]<- genernaseq$lymphnode[row]/(result$txEnd[row]- result$txStart[row]+1)}

for(row in seq_len(nrow(result))){
  skelmuscle_data[row]<- genernaseq$skelmuscle[row]/(result$txEnd[row]- result$txStart[row]+1)}

for(row in seq_len(nrow(result))){
  testes_data[row]<- genernaseq$testes[row]/(result$txEnd[row]- result$txStart[row]+1)}

#value calculation for second step of normalization calculation
b_val<- sum(brain_data)/1000000
brst_val<- sum(breast_data)/1000000 
c_val<- sum(colon_data)/1000000  
h_val<- sum(heart_data)/1000000 
liv_val<- sum(liver_data)/1000000  
lymph_val<- sum(lymphnode_data)/1000000  
skel_val<- sum(skelmuscle_data)/1000000  
t_val<- sum(testes_data)/1000000

#create empty vectors to put normalized data into
brain_norm <- c()
breast_norm <- c()
colon_norm <- c() 
heart_norm <- c()
liver_norm <- c()
lymphnode_norm <- c()
skelmuscle_norm <- c()
testes_norm <- c()

#loops go through each row for a particular tissue type and 
#do final step of normalization calculation
for(row in seq_len(nrow(result))){
  brain_norm[row]<- genernaseq$brain[row]/b_val} 

for(row in seq_len(nrow(result))){
  breast_norm[row]<- genernaseq$breast[row]/brst_val} 

for(row in seq_len(nrow(result))){
  colon_norm[row]<- genernaseq$colon[row]/c_val}

for(row in seq_len(nrow(result))){
  heart_norm[row]<- genernaseq$heart[row]/h_val}

for(row in seq_len(nrow(result))){
  liver_norm[row]<- genernaseq$liver[row]/liv_val}

for(row in seq_len(nrow(result))){
  lymphnode_norm[row]<- genernaseq$lymphnode[row]/lymph_val}

for(row in seq_len(nrow(result))){
  skelmuscle_norm[row]<- genernaseq$skelmuscle[row]/skel_val}

for(row in seq_len(nrow(result))){
  testes_norm[row]<- genernaseq$testes[row]/t_val}

names2<- c(result$name2) #creates names vector

#create matrix of normalized data
mat<- cbind(name_id, names2, brain_norm , breast_norm, colon_norm, heart_norm, 
                      liver_norm, lymphnode_norm, skelmuscle_norm, testes_norm, row.names= name_id)

genenorm<- data.frame(mat) #convert to data frame so can be entered into MYSQL table

#reconnect to my database
con2 <- dbConnect(MySQL(),
                  user = 'bhm263',
                  password = 'FJKVC9su5LwFLMaK',
                  host = 'orion.bio.nyu.edu',
                  dbname='bhm263')

dbSendQuery(con2, "DROP TABLE IF EXISTS `GeneRNAseqNorm`;")
dbWriteTable(con2, 'GeneRNAseqNorm', genenorm) #create table of normalized data

#Part 4

#calculate the averages of each gene 
x<- dbGetQuery(con2, "SELECT names2, (brain_norm + breast_norm +colon_norm + heart_norm +
                      liver_norm + lymphnode_norm + skelmuscle_norm + testes_norm) / 8 as row_avg FROM GeneRNAseqNorm 
                      ORDER BY row_avg DESC;")

x

#Part 5
#HUWE1 was the highest value from the table in part4, so use that for the query

#queries the normalized data for the gene that can be entered into barplot
bar_data<- dbGetQuery(con2, "SELECT brain_norm, breast_norm, colon_norm, heart_norm,
                      liver_norm, lymphnode_norm, skelmuscle_norm, testes_norm FROM GeneRNAseqNorm 
                      WHERE names2 = 'HUWE1';")

barplot(as.matrix(bar_data), cex.names = 0.5) #makes barplot

#Part 6

prot_data<- read.delim('BIOGRID.txt', header=FALSE)

#isolate columns of interest and put into dataframe
interactorA<- as.character(prot_data$V8) 
interactorB<- as.character(prot_data$V9)
z<- cbind(interactorA, interactorB)
z1<- as.data.frame(z)

#create an abridged table for query
dbSendQuery(con2, "DROP TABLE IF EXISTS `Interaction_data`;")
prot_table<- dbWriteTable(con2, 'Interaction_data', value= z1)

#creates a query that variable(gene names) can be added into
interaction_query <- "SELECT interactorA, interactorB, COUNT(*) from Interaction_data 
                    WHERE interactorA IN (%s) OR interactorB IN (%s)
                    GROUP BY interactorA
                    ORDER BY COUNT(*) DESC"

#puts formatted variables into query utilizing earlier function
int_query <- sprintf(interaction_query, addVals(names1), addVals(names1))


int_val<- dbGetQuery(con2, int_query) #execute query
int_val #Table shows that TP53 had the most connections

count_table<- dbWriteTable(con2, 'count_table', value= int_val)
count_valA<- dbGetQuery(con2, "SELECT interactorA, COUNT(interactorA) FROM count_table 
                       GROUP BY interactorA")

count_valA

count_valB<- dbGetQuery(con2, "SELECT interactorB, COUNT(interactorB) FROM count_table 
                       GROUP BY interactorB
                        ORDER BY COUNT(interactorB) DESC")

count_valB #order shows TP53 has highest count

