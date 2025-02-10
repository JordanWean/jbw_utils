library(org.Mm.eg.db)

df = read.csv('genesofinterest.csv')

newline = c(gene = '',
            genename = '',
            notes = '')

if (!newline['gene'] %in% df$gene) {
  df = rbind(df, newline)
}

df = df[order(df$gene), ]

# fix any duplications
df = df[!duplicated(df$gene), ]

annotations = select(
  org.Mm.eg.db,
  keys = df$gene,
  columns = "GENENAME",
  keytype = "SYMBOL"
)

df$genename = annotations$GENENAME

write.csv(x = df,
          file = 'genesofinterest.csv',
          row.names = F)
