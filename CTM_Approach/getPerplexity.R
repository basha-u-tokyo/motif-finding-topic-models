library("NLP")
library("tm")
library("lda")
library("topicmodels")

vocabFileName <- "ctm_vocab.txt"
docsFileName <- "ctm_seqs.txt"
k <- 10

vocab <- read.vocab(vocabFileName)
docs <- read.documents(docsFileName)
    
dtm <- ldaformat2dtm(docs, vocab, omit_empty = TRUE)
    
SEED <- 2010
    
CTM = CTM(dtm, k = k, control = list(seed = SEED, var = list(tol = 10^-4),
                                         em = list(tol = 10^-3)))
    
perp <- perplexity(CTM)

perp