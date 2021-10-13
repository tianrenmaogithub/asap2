library(vegan)
args=commandArgs(TRUE)
taxa=args[1]
meta=args[2]

# function to remove NA columns and rows from feature table and env data
removeNaRowCol=function(df){
  df=df[,colMeans(is.na(df))<0.5]
  df=df[rowMeans(is.na(df))<0.5,]
  df=df[,colSums(is.na(df))==0]
  df=df[rowSums(is.na(df))==0,]
  return(df)
}

# preprocess of feature table removing NA
spec=read.csv(taxa,header=TRUE,sep='\t',row.names=1)
#spec=as.data.frame(t(spec))
print('# samples    # species')
print(dim(spec))
print('Removing rows and columns with NA')
spec=removeNaRowCol(spec)
print('Removing all-zero columns and rows')
spec=spec[,colSums(spec)>0]
spec=spec[rowSums(spec)>0,]
print('# samples    # species')
print(dim(spec))

# preprocess of env data removing NA and env factors with only one level
env=read.csv(meta,header=TRUE,sep='\t',row.names=1)
print('# samples    # env')
print(dim(env))
print('Removing rows and columns with NA')
env=removeNaRowCol(env)
print('# samples    # env')
print(dim(env))
print('Removing env columns with only one level')
env=env[,apply(env,2,function(x) length(unique(x))>1)]
print('# samples    # env')
print(dim(env))

# Remove non-shared samples of feature table and env table
print('Removing non-shared samples')
intersect=intersect(rownames(spec),rownames(env))
print(paste('# Intersection samples: ',length(intersect)))
spec=spec[intersect,]
env=env[intersect,]

# RDA variable selection using ordistep, plot with significant env factors
sink(args[7],append=TRUE)
mod0 <- rda(spec ~ 1, env)  # Model with intercept only
mod1 <- rda(spec ~ ., env)  # Model with all explanatory variables
mod2 <- ordistep(mod0, scope = formula(mod1))
sink()

sink(args[3])
mod2
mod2$anova
sink()

pdf(args[4])
plot(mod2)
dev.off()

# CCA variable selection using ordistep, plot with significant env factors
sink(args[7],append=TRUE)
mod0 <- cca(spec ~ 1, env)  # Model with intercept only
mod1 <- cca(spec ~ ., env)  # Model with all explanatory variables
mod2 <- ordistep(mod0, scope = formula(mod1))
sink()

sink(args[5])
mod2
mod2$anova
sink()

pdf(args[6])
plot(mod2)
dev.off()
