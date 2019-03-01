#Homework2
#Part B
expvalues <- read.table("expvalues.txt", header = TRUE)
a <- colnames(expvalues)
Experiment <- data.frame(expid=1:length(a),expname=a) #create table of values for Experiment SQLite table
b<- rownames(expvalues)
probeid <- data.frame(probeid=1:length(b), probename=b) #create table of values for Probe SQLite table

library("RSQLite")
drv<-dbDriver("SQLite") #create driver
con<-dbConnect(drv, "B_bhm263.sqlite") #set up connection

dbWriteTable(con, "Experiment", Experiment) #create experiment table
dbWriteTable(con, "Probes", probeid) #create probe table

if (dbExistsTable(con, "data")) dbSendQuery(con, "drop table data") #delete table data if already exists

#Creates the empty table for data to be entered at a later point
dbSendQuery(con, "CREATE TABLE data (dataid INTEGER PRIMARY KEY NOT NULL ,
expid INTEGER NOT NULL,
probeid INTEGER NOT NULL,
expvalue FLOAT);")

dataid_val<- 1:length(expid_val) #Creates the primary key
z<-1:6
expid_val<- rep(z, each=263) #Creates a vector for FK of experimental id, I am entering by each experiment
probeid_val<-1:length(exp_val$Control1) #creates a vector of probe id values
#creates vector of experimental values, grouped by experiment
expvalue_val<- c(expvalues$Control1, expvalues$Control2, expvalues$Control3,  expvalues$Treatment1,  expvalues$Treatment2, expvalues$Treatment3)

data1<- data.frame(dataid_val, expid_val, probeid_val, expvalue_val) #dataframe that matches how values should appear in SQLite table

#loop that iterates through complete dataframe and enters into sqlite
for(i in 1:nrow(data1)) {
  insert_statement = paste("INSERT INTO data (dataid,
  expid,
  probeid,
  expvalue)
  VALUES (\"",data1[i,1],"\",\"",
                           data1[i,2],"\",\"",
                           data1[i,3],"\",\"",
                           data1[i, 4],"\")", sep="")
  dbSendQuery(con, insert_statement)
}


#Part C
dbGetQuery (con, "SELECT probename, 
            avg(expvalue) from Probe, 
            data WHERE data.probeid=Probe.probeid group by data.probeid")

dbDisconnect(con)
