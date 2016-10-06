library("NLP")
library("tm")
library("lda")
library("topicmodels")

vocabFileName <- "ctm_vocab.txt"
docsFileName <- "ctm_seqs.txt"
resCTMFileName <- "result_CTM.txt"
k <- 10

args <- commandArgs(trailingOnly = TRUE)
numResults <- as.numeric(args[1])

vocab <- read.vocab(vocabFileName)
docs <- read.documents(docsFileName)

dtm <- ldaformat2dtm(docs, vocab, omit_empty = TRUE)

SEED <- 2010

CTM = CTM(dtm, k = k, control = list(seed = SEED, var = list(tol = 10^-4),
                                     em = list(tol = 10^-3)))

top.words <- terms(CTM, numResults)
write(top.words[,1], file = resCTMFileName, sep = "\t")

perp <- perplexity(CTM)

perp