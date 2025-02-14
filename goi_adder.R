library(org.Mm.eg.db)
library(AnnotationDbi)

df = read.csv('genesofinterest.csv')

newline = c(gene = '',
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