codon<-tolower(c("TTT","TTC","TTA","TTG","TCT","TCC","TCA","TCG","TAT","TAC","TAA",
         "TAG","TGT","TGC","TGA","TGG","CTT","CTC","CTA","CTG","CCT","CCC","CCA","CCG",
         "CAT","CAC","CAA","CAG","CGT","CGC","CGA","CGG","ATT","ATC","ATA","ATG","ACT",
         "ACC","ACA","ACG","AAT","AAC","AAA","AAG","AGT","AGC","AGA","AGG","GTT","GTC",
         "GTA","GTG","GCT","GCC","GCA","GCG","GAT","GAC","GAA","GAG","GGT","GGC","GGA","GGG"))

amino<-c("F","F","L","L","S","S","S","S","Y","Y","X","X","C","C","X","W",
         "L","L","L","L","P","P","P","P","H","H","Q","Q","R","R","R","R","I","I","I","M",
         "T","T","T","T","N","N","K","K","S","S","R","R","V","V","V","V","A","A","A","A",
         "D","D","E","E","G","G","G","G")

trans_from <- c("a", "t", "g", "c")

trans_to   <- c("t", "a", "c", "g")

