# This prepare positions which are not in the end file.
# The last portion in the split function is usually not this length and need to be prepared seperately

setwd("~/lustre/TCGA/MACH_output_EA")

pos1 = 1:110001
write.table(pos1, file = "pos1", quote =F, col.names = F, row.names =F)

for (p in 2:34){
    print(c(((p-1)*100000-10000 +1), (p*100000+10000 + 1)))
    pos = ((p-1)*100000-10000 +1):(p*100000+10000 + 1)
    write.table(pos, file = sprintf("pos%d", p),quote =F, col.names = F, row.names =F)
}
