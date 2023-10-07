

list1 = vector(mode = "list", 10L) #generate empty list with 10 entries
lists = vector(mode = "list", 6L) #generate empty list with 6 entries

a = 1
repeat{
  lists = lapply(lists, function(x) matrix(runif(18),nrow = 6, ncol = 3)) #add matrix with random numbers to each list entry
  list1[[a]] = lists #add list with matrices to another list
  if(a == length(list1)) break #stop condition if every entry in list1 is filled
  a = a + 1
}

#list1 now contains 10 elements, which each contain again another list with 6 elements

# Add two columns to each matrix: the first column should contain the index number of list1
# the second column should contain the index number of each matrix within that list 

for (i in 1:length(list1)){
  for(j in 1:length(list1[[i]])){
    list1[[i]][[j]]=cbind(list1[[i]][[j]],i,j)
  }
}


for(i in 1:length(list1))list1[[i]]= do.call(rbind,list1[[i]])
list1
result=do.call(rbind,list1)
write.table(result,"test1.csv",quote=FALSE,sep=";")
getwd()
?do.call
