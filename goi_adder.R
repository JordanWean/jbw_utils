df = read.csv('genesofinterest.csv')

newline = c(gene = 'Pmch', 
            name = 'Pro-melanin concentrating hormone',
            notes = 'MCH neurons in the hypothalamus project numerous places and have many effects on feeding, sleep, arousal')

if (!newline['gene'] %in% df$gene){
  df = rbind(df, newline)
}

df = df[order(df$gene), ]

write.csv(x = df, file = 'genesofinterest.csv')