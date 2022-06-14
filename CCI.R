###### Collagen interation network
### CCI analysis
library(igraph)
library(ggraph)
library(stringr)

CCI_list <- read.csv("CCI_list-col.csv")
Chemo = CCI_list
uni.chemo = list(Chemo[,1], Chemo[,2])
uni.chemo = unique(unlist(uni.chemo))

BREAST2 <- BREAST
Idents(BREAST2) <- BREAST@meta.data$CLASS

CCImerge = subset(BREAST2, subset = Tissue == "MammaryGland.Pregnancy")
CCImerge = CCImerge[uni.chemo,]
CCImerge.store = CCImerge

unigenes = rownames(CCImerge)
mel = list()
for(i in 1: length(unigenes)){

  expr <- FetchData(object = CCImerge, vars = unigenes[i])

  if(any(expr > 2) >= 1){

    ge = CCImerge[, which(x = expr > 2)]
    ge = table(ge@active.ident)
    ho = melt(ge)
    ho = mutate(ho, Gene = unigenes[i])
    mel = rbind(mel, ho)

  }else{
    print(unigenes[i])
  }
  
}

tocell = table(CCImerge@active.ident)
Ratiomel = list()

for(i in 1:length(tocell)){

  cell = names(tocell[i])
  LS = mel[mel[,1] == cell,]
  LS[,2] = LS[,2]/tocell[i]
  Ratiomel = rbind(Ratiomel, LS)

}
Ratiomel = dcast(Ratiomel, Var1 ~ Gene, value.var = "value")
write.table(Ratiomel, file = "CCIcol-RatioPregMAT.csv", sep = ",")

A = list()
B = list()
CCItotal = list()
CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(mel[,3])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(mel[,3])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "|"))

ratiomelt = na.omit(melt(Ratiomel))
ratiomelt = ratiomelt[ratiomelt$value > 0.1,]

CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "&"))
a = as.character(unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[CCI_perlist$merge %in% unique(CCI_perlist$merge),]
CCI_perlist = subset(CCI_perlist, CCI_perlist$merge %in% a)
a = unique(CCI_perlist$merge)
aa = str_split(a, "&", n = 2, simplify = TRUE)
CCI_perlist = data.frame(L = aa[,1], R = aa[,2], merge = unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[!(as.character(CCI_perlist$L) == as.character(CCI_perlist$R)),]

CCItotal = list()
for(i in 1:nrow(CCI_perlist)){

  geneset = CCI_perlist[i,]
  Ligand = as.character(unlist(geneset[1]))
  Receptor = as.character(unlist(geneset[2]))

  LigCCI = ratiomelt[ratiomelt$variable == Ligand, ]
  ResCCI = ratiomelt[ratiomelt$variable == Receptor, ]

  for(k in 1:nrow(LigCCI)){

    Lig = LigCCI[k,3]
    score = lapply(ResCCI$value, function(x){x * Lig})
    Res = data.frame( CCIscore = unlist(score), Recepter.Cell = ResCCI$Var1)
    Res = mutate(Res, Ligand.Cell = LigCCI[k,1], CCI = geneset)

    CCItotal = bind_rows(CCItotal, Res)
  }

}
CCItotal = CCItotal[!(as.character(CCItotal$Recepter.Cell) == as.character(CCItotal$Ligand.Cell)),]
write.table(CCItotal, file = "CCIcol-total-Preg.csv", sep = ",")

CCItop = CCItotal[CCItotal$CCIscore > 0.3,]
write.table(CCItop, file = "CCIcol-Pregover0.3.csv", sep = ",")

Lig = unlist(as.character(unique(CCItotal$Ligand.Cell)))
CellCell.CCI = list()
for(i in 1:length(Lig)){

  G <- Lig[i]
  set = CCItotal[CCItotal$Ligand.Cell == G, ]

  CELLset = unique(set$Recepter.Cell)

  for(k in 1:length(CELLset)){

    R <- CELLset[k]
    REset <- set[set$Recepter.Cell == R, ]
    CCIsum <- sum(REset$CCIscore)

    DF = data.frame(LigandCell = G, ReceptorCell = R, CCIscore.sum = CCIsum)

    CellCell.CCI = rbind(CellCell.CCI, DF)

  }

}
PregCellCCI = CellCell.CCI
write.table(CellCell.CCI, file = "CCIcolCELL-Preg.csv", sep = ",")

#######################
CCImerge = subset(BREAST2, subset = Tissue == "MammaryGland.Involution")
CCImerge = CCImerge[uni.chemo,]
CCImerge.store = CCImerge

unigenes = rownames(CCImerge)
mel = list()
for(i in 1: length(unigenes)){

  expr <- FetchData(object = CCImerge, vars = unigenes[i])

  if(any(expr > 2) >= 1){

    ge = CCImerge[, which(x = expr > 2)]
    ge = table(ge@active.ident)
    ho = melt(ge)
    ho = mutate(ho, Gene = unigenes[i])
    mel = rbind(mel, ho)

  }else{
    print(unigenes[i])
  }
  
}

tocell = table(CCImerge@active.ident)
Ratiomel = list()
for(i in 1:length(tocell)){

  cell = names(tocell[i])
  LS = mel[mel[,1] == cell,]
  LS[,2] = LS[,2]/tocell[i]
  Ratiomel = rbind(Ratiomel, LS)

}
Ratiomel = dcast(Ratiomel, Var1 ~ Gene, value.var = "value")

write.table(Ratiomel, file = "CCcol-RatioInvMAT.csv", sep = ",")

A = list()
B = list()
CCItotal = list()
CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(mel[,3])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(mel[,3])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "|"))

ratiomelt = na.omit(melt(Ratiomel))
ratiomelt = ratiomelt[ratiomelt$value > 0.1,]

CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "&"))
a = as.character(unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[CCI_perlist$merge %in% unique(CCI_perlist$merge),]
CCI_perlist = subset(CCI_perlist, CCI_perlist$merge %in% a)
a = unique(CCI_perlist$merge)
aa = str_split(a, "&", n = 2, simplify = TRUE)
CCI_perlist = data.frame(L = aa[,1], R = aa[,2], merge = unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[!(as.character(CCI_perlist$L) == as.character(CCI_perlist$R)),]

CCItotal = list()
for(i in 1:nrow(CCI_perlist)){

  geneset = CCI_perlist[i,]
  Ligand = as.character(unlist(geneset[1]))
  Receptor = as.character(unlist(geneset[2]))

  LigCCI = ratiomelt[ratiomelt$variable == Ligand, ]
  ResCCI = ratiomelt[ratiomelt$variable == Receptor, ]

  for(k in 1:nrow(LigCCI)){

    Lig = LigCCI[k,3]
    score = lapply(ResCCI$value, function(x){x * Lig})
    Res = data.frame( CCIscore = unlist(score), Recepter.Cell = ResCCI$Var1)
    Res = mutate(Res, Ligand.Cell = LigCCI[k,1], CCI = geneset)

    CCItotal = bind_rows(CCItotal, Res)
  }

}
CCItotal = CCItotal[!(as.character(CCItotal$Recepter.Cell) == as.character(CCItotal$Ligand.Cell)),]
write.table(CCItotal, file = "CCIcol-total-Inv.csv", sep = ",")

CCItop = CCItotal[CCItotal$CCIscore > 0.3,]
write.table(CCItop, file = "CCIcol-Invover0.3.csv", sep = ",")

Lig = unlist(as.character(unique(CCItotal$Ligand.Cell)))
CellCell.CCI = list()
for(i in 1:length(Lig)){

  G <- Lig[i]
  set = CCItotal[CCItotal$Ligand.Cell == G, ]

  CELLset = unique(set$Recepter.Cell)

  for(k in 1:length(CELLset)){

    R <- CELLset[k]
    REset <- set[set$Recepter.Cell == R, ]
    CCIsum <- sum(REset$CCIscore)

    DF = data.frame(LigandCell = G, ReceptorCell = R, CCIscore.sum = CCIsum)

    CellCell.CCI = rbind(CellCell.CCI, DF)

  }

}
InvCellCCI = CellCell.CCI

write.table(CellCell.CCI, file = "CCIcolCELL-Inv.csv", sep = ",")

#######################
CCImerge = subset(BREAST2, subset = Tissue == "MammaryGland.Lactation")
CCImerge = CCImerge[uni.chemo,]
CCImerge.store = CCImerge

unigenes = rownames(CCImerge)
mel = list()
for(i in 1: length(unigenes)){

  expr <- FetchData(object = CCImerge, vars = unigenes[i])

  if(any(expr > 2) >= 1){

    ge = CCImerge[, which(x = expr > 2)]
    ge = table(ge@active.ident)
    ho = melt(ge)
    ho = mutate(ho, Gene = unigenes[i])
    mel = rbind(mel, ho)

  }else{
    print(unigenes[i])
  }
  
}

tocell = table(CCImerge@active.ident)
Ratiomel = list()
for(i in 1:length(tocell)){

  cell = names(tocell[i])
  LS = mel[mel[,1] == cell,]
  LS[,2] = LS[,2]/tocell[i]
  Ratiomel = rbind(Ratiomel, LS)

}
Ratiomel = dcast(Ratiomel, Var1 ~ Gene, value.var = "value")

write.table(Ratiomel, file = "CCIcol-Ratio-LactMAT.csv", sep = ",")

A = list()
B = list()
CCItotal = list()
CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(mel[,3])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(mel[,3])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "|"))

ratiomelt = na.omit(melt(Ratiomel))
ratiomelt = ratiomelt[ratiomelt$value > 0.1,]

CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "&"))
a = as.character(unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[CCI_perlist$merge %in% unique(CCI_perlist$merge),]
CCI_perlist = subset(CCI_perlist, CCI_perlist$merge %in% a)
a = unique(CCI_perlist$merge)
aa = str_split(a, "&", n = 2, simplify = TRUE)
CCI_perlist = data.frame(L = aa[,1], R = aa[,2], merge = unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[!(as.character(CCI_perlist$L) == as.character(CCI_perlist$R)),]

CCItotal = list()
for(i in 1:nrow(CCI_perlist)){

  geneset = CCI_perlist[i,]
  Ligand = as.character(unlist(geneset[1]))
  Receptor = as.character(unlist(geneset[2]))

  LigCCI = ratiomelt[ratiomelt$variable == Ligand, ]
  ResCCI = ratiomelt[ratiomelt$variable == Receptor, ]

  for(k in 1:nrow(LigCCI)){

    Lig = LigCCI[k,3]
    score = lapply(ResCCI$value, function(x){x * Lig})
    Res = data.frame( CCIscore = unlist(score), Recepter.Cell = ResCCI$Var1)
    Res = mutate(Res, Ligand.Cell = LigCCI[k,1], CCI = geneset)

    CCItotal = bind_rows(CCItotal, Res)
  }

}
CCItotal = CCItotal[!(as.character(CCItotal$Recepter.Cell) == as.character(CCItotal$Ligand.Cell)),]

write.table(CCItotal, file = "CCIcol-total-Lact.csv", sep = ",")

CCItop = CCItotal[CCItotal$CCIscore > 0.3,]
write.table(CCItop, file = "CCIcol-Lactover0.3.csv", sep = ",")

Lig = unlist(as.character(unique(CCItotal$Ligand.Cell)))
CellCell.CCI = list()
for(i in 1:length(Lig)){

  G <- Lig[i]
  set = CCItotal[CCItotal$Ligand.Cell == G, ]

  CELLset = unique(set$Recepter.Cell)

  for(k in 1:length(CELLset)){

    R <- CELLset[k]
    REset <- set[set$Recepter.Cell == R, ]
    CCIsum <- sum(REset$CCIscore)

    DF = data.frame(LigandCell = G, ReceptorCell = R, CCIscore.sum = CCIsum)

    CellCell.CCI = rbind(CellCell.CCI, DF)

  }

}
LactCellCCI = CellCell.CCI

write.table(CellCell.CCI, file = "CCIcolCELL-Lact.csv", sep = ",")

#######################
CCImerge = subset(BREAST2, subset = Tissue == "MammaryGland.Virgin")
CCImerge = CCImerge[uni.chemo,]
CCImerge.store = CCImerge

unigenes = rownames(CCImerge)
mel = list()
for(i in 1: length(unigenes)){

  expr <- FetchData(object = CCImerge, vars = unigenes[i])

  if(any(expr > 2) >= 1){

    ge = CCImerge[, which(x = expr > 2)]
    ge = table(ge@active.ident)
    ho = melt(ge)
    ho = mutate(ho, Gene = unigenes[i])
    mel = rbind(mel, ho)

  }else{
    print(unigenes[i])
  }
  
}

tocell = table(CCImerge@active.ident)
Ratiomel = list()
for(i in 1:length(tocell)){

  cell = names(tocell[i])
  LS = mel[mel[,1] == cell,]
  LS[,2] = LS[,2]/tocell[i]
  Ratiomel = rbind(Ratiomel, LS)

}
Ratiomel = dcast(Ratiomel, Var1 ~ Gene, value.var = "value")

write.table(Ratiomel, file = "CCIcol-Ratio-VirMAT.csv", sep = ",")

A = list()
B = list()
CCItotal = list()
CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(mel[,3])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(mel[,3])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "|"))

ratiomelt = na.omit(melt(Ratiomel))
ratiomelt = ratiomelt[ratiomelt$value > 0.1,]

CCI_perlist = Chemo[Chemo$AliasA %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = CCI_perlist[CCI_perlist$AliasB %in% unlist(unique(ratiomelt[,2])),]
CCI_perlist = mutate(CCI_perlist, merge = paste(CCI_perlist$AliasA, CCI_perlist$AliasB, sep = "&"))
a = as.character(unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[CCI_perlist$merge %in% unique(CCI_perlist$merge),]
CCI_perlist = subset(CCI_perlist, CCI_perlist$merge %in% a)
a = unique(CCI_perlist$merge)
aa = str_split(a, "&", n = 2, simplify = TRUE)
CCI_perlist = data.frame(L = aa[,1], R = aa[,2], merge = unique(CCI_perlist$merge))
CCI_perlist = CCI_perlist[!(as.character(CCI_perlist$L) == as.character(CCI_perlist$R)),]

CCItotal = list()
for(i in 1:nrow(CCI_perlist)){

  geneset = CCI_perlist[i,]
  Ligand = as.character(unlist(geneset[1]))
  Receptor = as.character(unlist(geneset[2]))

  LigCCI = ratiomelt[ratiomelt$variable == Ligand, ]
  ResCCI = ratiomelt[ratiomelt$variable == Receptor, ]

  for(k in 1:nrow(LigCCI)){

    Lig = LigCCI[k,3]
    score = lapply(ResCCI$value, function(x){x * Lig})
    Res = data.frame( CCIscore = unlist(score), Recepter.Cell = ResCCI$Var1)
    Res = mutate(Res, Ligand.Cell = LigCCI[k,1], CCI = geneset)

    CCItotal = bind_rows(CCItotal, Res)
  }

}

CCItotal = CCItotal[!(as.character(CCItotal$Recepter.Cell) == as.character(CCItotal$Ligand.Cell)),]
write.table(CCItotal, file = "CCIcol-total-Vir.csv", sep = ",")

CCItop = CCItotal[CCItotal$CCIscore > 0.3,]
write.table(CCItop, file = "CCIcol-Virover0.3.csv", sep = ",")

Lig = unlist(as.character(unique(CCItotal$Ligand.Cell)))
CellCell.CCI = list()
for(i in 1:length(Lig)){

  G <- Lig[i]
  set = CCItotal[CCItotal$Ligand.Cell == G, ]

  CELLset = unique(set$Recepter.Cell)

  for(k in 1:length(CELLset)){

    R <- CELLset[k]
    REset <- set[set$Recepter.Cell == R, ]
    CCIsum <- sum(REset$CCIscore)

    DF = data.frame(LigandCell = G, ReceptorCell = R, CCIscore.sum = CCIsum)

    CellCell.CCI = rbind(CellCell.CCI, DF)

  }

}
VirCellCCI = CellCell.CCI

write.table(CellCell.CCI, file = "CCIcolCELL-Vir.csv", sep = ",")



### graph visualization
library(circlize)
library(grDevices)

grid.col = c(Luminal = "yellow", LuminalProgenitor = "yellow", StemProgenitor = "yellow", DividingEpithelial = "yellow", NK = "green", Tcell = "green", Bcell = "green", Macrophage = "green", Dendriticcell = "green",
  Stromal = "blue", Adipocyte = "blue", SecretoryAlveoli = "yellow", Myeloid = "green", BasalMyoepithelial = "yellow", DuctalLuminal = "yellow", Endothelilal = "blue", Muscle = "blue")
col_fun = c("mistyrose", "pink", "red3")


chordDiagram(VirCellCCI, col = col_fun, grid.col = grid.col, link.border = "pink", link.lty = 2, link.lwd = 1, transparency = 0.15, directional = 1, target.prop.height = mm_h(12), diffHeight = mm_h(5))
chordDiagram(PregCellCCI, col = col_fun, grid.col = grid.col, link.border = "pink", link.lty = 2, link.lwd = 1, transparency = 0.15, directional = 1, target.prop.height = mm_h(12), diffHeight = mm_h(5))
chordDiagram(LactCellCCI, col = col_fun, grid.col = grid.col, link.border = "pink", link.lty = 2, link.lwd = 1, transparency = 0.15, directional = 1, target.prop.height = mm_h(12), diffHeight = mm_h(5))
chordDiagram(InvCellCCI, col = col_fun, grid.col = grid.col, link.border = "pink", link.lty = 2, link.lwd = 1, transparency = 0.15, directional = 1, target.prop.height = mm_h(12), diffHeight = mm_h(5))

