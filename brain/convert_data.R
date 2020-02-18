

for (i in 1:42) 
{
  name = paste("brain_data/data_column",toString(i),sep="")
  write(t(tensorA[,,i]), file=name, ncolumns=68)
}