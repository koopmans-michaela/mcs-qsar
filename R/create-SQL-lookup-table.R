lookupdb <- dbConnect(RSQLite::SQLite(),"")
dbWriteTable(lookupdb, "testlookup", testlookup)
dbListTables(lookupdb)
testlookup2 <- testlookup
testlookup2$ID <- seq.int(nrow(testlookup2))
lookupdb2 <- dbConnect(RSQLite::SQLite(),"")
dbWriteTable(lookupdb2, "testlookup2", testlookup2)
