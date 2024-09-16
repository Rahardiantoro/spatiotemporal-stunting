###Title: Spatio-temporal modeling to identify factors associated with stunting in Indonesia using a Modified Generalized Lasso
###Authors: Septian Rahardiantoro1*, Alfidhia Rahman Nasa Juhanda, Anang Kurnia, 
####Aswi Aswi, Bagus Sartono, Dian Handayani, Agus Mohamad Soleh, Yusma Yanti, Susanna Cramb


###Package used
library(sf)
library(spdep)
library(genlasso)
library(RColorBrewer)
library(ggrepel)
library(heatmaply)
library(MASS)
library(dbplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(maps)

###Dataset of case study
dt <- read.csv("20240915 data stunting v4.csv",sep=";")
dt$Province <- as.factor(dt$Province)
dt$Year <- as.factor(dt$Year)

###Indonesian map
ina <- st_read(dsn="idn_adm_bps_20200401_shp",layer="idn_admbnda_adm1_bps_20200401")
ina <- st_make_valid(ina)
ina$ADM1_EN <- ifelse(ina$ADM1_EN == "Daerah Istimewa Yogyakarta", "DI Yogyakarta", ina$ADM1_EN)

###Connection of Provinces
cod <- read.csv("20240915 Coordinates of Provinces.csv",sep=";") #coordinates data
coords <- cod[,c(4,3)]
nbq <- poly2nb(st_geometry(ina), queen=TRUE) #queen method
nbq_new <- nbq
nbq_new[[18]] <- 0L
nbq_new[[8]] <- nbq_new[[8]][which(nbq_new[[8]]!=18)]

nbk2 <- knn2nb(knearneigh(coords, k=2)) #k-nn (k=2) method
nbk3 <- knn2nb(knearneigh(coords, k=3)) #k-nn (k=3) method

##Figure 1. Plot of adjacency between provinces
layout(matrix(c(1,2,3,3),2,2,byrow = TRUE))

plot(st_geometry(ina)) #queen
plot(nbq_new, coords, col='#c70233', lwd=1.5, add=TRUE)
points(coords, pch=19, col='#c70233')
title(sub="(a)")

plot(st_geometry(ina)) #k-nn (k=2)
plot(nbk2, coords, col='#c70233', pch=19, lwd=1.5, add=TRUE)
points(coords, pch=19, col='#c70233')
title(sub="(b)")

plot(st_geometry(ina)) #k-nn (k=3)
plot(nbk3, coords, col='#c70233', pch=19, lwd=1.5, add=TRUE)
points(coords, pch=19, col='#c70233')
title(sub="(c)")

##Figure 2. Boxplot of stunting prevalence in Indonesia
ggplot(dt, aes(x = Year, y = Stunting.Prevalence)) +
  geom_boxplot(fill="lightblue") +
  coord_flip() +
  theme_classic() +
  geom_text_repel(
    data = dt %>%
      group_by(Year) %>%
      mutate(
        Label = case_when(
          Stunting.Prevalence == min(Stunting.Prevalence) ~ as.character(Province),
          Stunting.Prevalence == max(Stunting.Prevalence) ~ as.character(Province),
          (Stunting.Prevalence < (quantile(Stunting.Prevalence, 0.25) - 1.5 * IQR(Stunting.Prevalence)) | 
             Stunting.Prevalence > (quantile(Stunting.Prevalence, 0.75) + 1.5 * IQR(Stunting.Prevalence))-1) ~ as.character(Province),
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(Label)),  # This ensures that only the desired labels are used
    aes(label = Label),
    size = 3.3,
    max.overlaps = Inf)

##Figure 3. Map of Indonesia's stunting diversity
ina19 <- ina; ina20 <- ina; ina21 <- ina; ina22 <- ina; ina23 <- ina 
ina19$long <- cod$longitude; ina19$lat <- cod$latitude
ina20$long <- cod$longitude; ina20$lat <- cod$latitude
ina21$long <- cod$longitude; ina21$lat <- cod$latitude
ina22$long <- cod$longitude; ina22$lat <- cod$latitude
ina23$long <- cod$longitude; ina23$lat <- cod$latitude

###data preparation
ina19$stunting <- dt[dt$Year=="2019",4]
ina20$stunting <- dt[dt$Year=="2020",4]
ina21$stunting <- dt[dt$Year=="2021",4]
ina22$stunting <- dt[dt$Year=="2022",4]
ina23$stunting <- dt[dt$Year=="2023",4]

cols <-function(n) {
  colorRampPalette(rev(c("red4","red2","tomato2","orange","gold1","forestgreen")))(6)                                
}

map1 <- ggplot(data = ina19) +
  geom_sf(aes(fill = stunting, geometry=geometry), color = "black",size = 0.08, lwd=.1) +
  scale_fill_gradientn(colours = cols(6)) +
  theme_classic() + labs(caption = "(a)")

map2 <- ggplot(data = ina20) +
  geom_sf(aes(fill = stunting, geometry=geometry), color = "black",size = 0.08, lwd=.1) +
  scale_fill_gradientn(colours = cols(6)) +
  theme_classic() + labs(caption = "(b)")

map3 <- ggplot(data = ina21) +
  geom_sf(aes(fill = stunting, geometry=geometry), color = "black",size = 0.08, lwd=.1) +
  scale_fill_gradientn(colours = cols(6)) +
  theme_classic() + labs(caption = "(c)")

map4 <- ggplot(data = ina22) +
  geom_sf(aes(fill = stunting, geometry=geometry), color = "black",size = 0.08, lwd=.1) +
  scale_fill_gradientn(colours = cols(6)) +
  theme_classic() + labs(caption = "(d)")

map5 <- ggplot(data = ina23) +
  geom_sf(aes(fill = stunting, geometry=geometry), color = "black",size = 0.08, lwd=.1) +
  scale_fill_gradientn(colours = cols(6)) +
  theme_classic() + labs(caption = "(e)")

###Function to extract legend
g_legend <- function(plot) {
  tmp <- ggplot_gtable(ggplot_build(plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) == 0) {
    stop("No legend found in the plot.")
  }
  return(tmp$grobs[[leg]])
}

###Extract legend from one of the plots
legend1 <- g_legend(map2)

###Hide legends in individual plots for clean appearance
map1 <- map1 + theme(legend.position = "none")
map2 <- map2 + theme(legend.position = "none")
map3 <- map3 + theme(legend.position = "none")
map4 <- map4 + theme(legend.position = "none")
map5 <- map5 + theme(legend.position = "none")

###Combine maps and legend
grid.arrange(
  arrangeGrob(map1, map2, map3, map4, map5,legend1, ncol = 3)
)

####Figure A. Density plot of stunting prevalence
den1 <- ggplot(data.frame(dt$Stunting.Prevalence), aes(x = dt$Stunting.Prevalence)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = diff(range(dt$Stunting.Prevalence))/10, 
                 color = "black", fill = "grey85") +
  geom_density(alpha = .1, fill = "blue") + theme_minimal() +
  labs(x='stunting', caption="(a)")

den2 <- ggplot(data.frame(sqrt(dt$Stunting.Prevalence)), aes(x = sqrt(dt$Stunting.Prevalence))) +
  geom_histogram(aes(y = after_stat(density)), binwidth = diff(range(sqrt(dt$Stunting.Prevalence)))/10, 
                 color = "black", fill = "grey85") +
  geom_density(alpha = .1, fill = "blue") + theme_minimal() +
  labs(x='stunting', caption="(b)")

grid.arrange(
  arrangeGrob(den1, den2,ncol = 2)
)

###Generalized Lasso Preparation
dt2 <- dt[,-c(1:3)]
label <- c("Stunting (Y)","Poverty (X1)","EBF (X2)",
           "LBW (X3)", "High Sch Comp (X4)","Prop. Sanitation Acc (X5)",
           "Unmet need (X6)","GDP (X7)", "Health Svc Rate (X8)",
           "Health Care Rate (X9)")
names(dt2) <- label
str(dt2)

####Transformation Y
datX <- dt2
datX$stunting <- sqrt(dt2$`Stunting (Y)`)
datX <- datX %>% rename("sqrt_stunting"="stunting")
datX <- datX[,-1]

####X standardisation
predictors <- names(datX)[1:(length(datX)-1)]
datX[predictors] <- scale(datX[predictors])

##Figure 4. Correlation plot
par(mfrow = c(1, 1))
####after transformation and standardisation
corrplot::corrplot(cor(datX,use="na.or.complete"), method = "color", 
                   type = "upper", tl.col = "black", 
                   tl.srt = 90, addCoef.col = "black", diag=F)

###X matrix formation
for(i in 1:9) {
  # Construct the variable name
  var_name <- paste0("X", i)
  # Extract the diagonal values
  diagonal_values <- diag(datX[predictors[i]] %>% unlist())
  # Assign the extracted values to the dynamically named variable
  assign(var_name, diagonal_values)
}
Xo <- cbind(X1,X2,X3,X4,X5,X6,X7,X8,X9) %>% as.matrix()
dim(Xo)
#[1]  170 1530

#Remove columns with any NA values
X <- Xo[, !apply(is.na(Xo), 2, any)]
dim(X)
#[1]  170 1428

###centering Y
y <- datX$sqrt_stunting
y <- y-mean(datX$sqrt_stunting) #no intercept
y <- as.matrix(y)

###Contruct D penalty matrices
D_Matrix_Panel <- function(D, years){
  # Create base matrix and vector
  prov <- colnames(D) # Region Names
  D0 <- as.matrix(D) # Base D Matrix
  I1 <- diag(-1,ncol(D))
  I2 <- diag(1,ncol(D))
  o <- diag(0,ncol(D))
  
  # Create Upper Matrix
  if(length(years)==5){
    prov_matrix <- as.matrix(bdiag(D0,D0,D0,D0,D0))}
  else if(length(years)==4){
    prov_matrix <- as.matrix(bdiag(D0,D0,D0,D0))}
  else {
    prov_matrix <- as.matrix(bdiag(D0,D0,D0))}
  
  # Create Lower Matrix
  if(length(years)==5){
    year_matrix <- rbind(cbind(I1,I2,o,o,o), cbind(o,I1,I2,o,o),
                       cbind(o,o,I1,I2,o), cbind(o,o,o,I1,I2))}
  else if(length(years)==4){
    year_matrix <- rbind(cbind(I1,I2,o,o), cbind(o,I1,I2,o),
                         cbind(o,o,I1,I2))
  }
  else{
    year_matrix <- rbind(cbind(I1,I2,o), cbind(o,I1,I2))
  }
  
  D_Matrix <- rbind(prov_matrix,year_matrix)
  
  #Construct new column names for the new D Matrix
  new_column_names <- paste(expand.grid(prov = prov, year = years)$year,
                            expand.grid(prov = prov, year = years)$prov,
                            sep = "_")
  
  colnames(D_Matrix) <- new_column_names
  
  return(D_Matrix)
}

matriksD.knn <- function(w3) {
  nk <- ncol(w3)
  indc <- which(w3 == 1,arr.ind = TRUE)
  indc.sort <- indc[order(indc[,1],decreasing=FALSE),]
  
  for(j in 1:nrow(indc.sort)){
    if(indc.sort[j,1]>indc.sort[j,2]) {
      a <- indc.sort[j,1]
      indc.sort[j,1] <- indc.sort[j,2]
      indc.sort[j,2] <- a
    }
  }
  indc.sort
  
  indc.sort.uniq <- indc.sort[!duplicated(indc.sort), ]
  
  D2 <- NULL
  for(i in 1:nrow(indc.sort.uniq)){
    dd <- rep(0,nk)
    dd[indc.sort.uniq[i,1]] <- c(-1)
    dd[indc.sort.uniq[i,2]] <- c(1)
    D2 <- rbind(D2,dd)
  }
  
  return(D2)
}

#D0 Queens
get_graph = function(shp_files) {
  nb_files = poly2nb(shp_files)
  edges = cbind(rep(1:length(nb_files), sapply(nb_files, length)),  unlist(nb_files))
  index = which(edges[,2] == 0)
  if(length(index) > 0) edges = edges[-index,]
  graph_files = graph.edgelist(edges, directed=FALSE)
  return(graph_files)}

graph_obj = get_graph(ina)

D_dgCMatrix = getDgSparse(graph_obj)
D0 <- as.matrix(D_dgCMatrix)
colnames(D0) <- ina$ADM1_EN
sum(duplicated(D0))

D0 <- D0 %>% as.data.frame() %>% distinct() %>% as.matrix()
dim(D0)
## [1] 33 34

#D queens
Dx1 <- D_Matrix_Panel(D0,years=2019:2023)
dim(Dx1)
Dx3 <- D_Matrix_Panel(D0,years=2020:2022)
dim(Dx3)
Dx5 <- D_Matrix_Panel(D0,years=2019:2022)
dim(Dx5)
D <- as.matrix(bdiag(Dx1,Dx1,Dx3,Dx1,Dx5,Dx1,Dx1,Dx1,Dx1)) #D queens
dim(D)
#[1] 2508 1428 

#D0 knn
w2 <-nb2mat(knn2nb(knearneigh(coords, k=2)), style='B')
w3 <-nb2mat(knn2nb(knearneigh(coords, k=3)), style='B')

D0.k2 <- matriksD.knn(w2)
dim(D0.k2)
#[1] 44 34

D0.k3 <- matriksD.knn(w3)
dim(D0.k3)
#[1] 66 34

sum(duplicated(D0.k2))
sum(duplicated(D0.k3))

Dx1.k2 <- D_Matrix_Panel(D0.k2,years=2019:2023)
Dx3.k2 <- D_Matrix_Panel(D0.k2,years=2020:2022)
Dx5.k2 <- D_Matrix_Panel(D0.k2,years=2019:2022)

Dx1.k3 <- D_Matrix_Panel(D0.k3,years=2019:2023)
Dx3.k3 <- D_Matrix_Panel(D0.k3,years=2020:2022)
Dx5.k3 <- D_Matrix_Panel(D0.k3,years=2019:2022)

D.k2 <- as.matrix(bdiag(Dx1.k2,Dx1.k2,Dx3.k2,Dx1.k2,Dx5.k2,
                        Dx1.k2,Dx1.k2,Dx1.k2,Dx1.k2))
dim(D.k2)
#[1] 2970 1428
D.k3 <- as.matrix(bdiag(Dx1.k3,Dx1.k3,Dx3.k3,Dx1.k3,Dx5.k3,
                        Dx1.k3,Dx1.k3,Dx1.k3,Dx1.k3))
dim(D.k3)
#[1] 3894 1428

###Generalized Ridge Modeling
#Evaluate for lambda <= 0.4
lambdas <- seq(0.1,0.4,len=50)
loocv_q <- c()
loocv_k2 <- c()
loocv_k3 <- c()
for(i in 1:length(lambdas)){
  beta_ridge_q <- ginv(t(X)%*%X+lambdas[i]*t(D)%*%D)%*%t(X)%*%y
  beta_ridge_k2 <- ginv(t(X)%*%X+lambdas[i]*t(D.k2)%*%D.k2)%*%t(X)%*%y
  beta_ridge_k3 <- ginv(t(X)%*%X+lambdas[i]*t(D.k3)%*%D.k3)%*%t(X)%*%y
  Hq <- X%*%ginv(t(X)%*%X+lambdas[i]*t(D)%*%D)%*%t(X)
  Hk2 <- X%*%ginv(t(X)%*%X+lambdas[i]*t(D.k2)%*%D.k2)%*%t(X)
  Hk3 <- X%*%ginv(t(X)%*%X+lambdas[i]*t(D.k3)%*%D.k3)%*%t(X)
  yduga_q <- X%*%beta_ridge_q
  yduga_k2 <- X%*%beta_ridge_k2
  yduga_k3 <- X%*%beta_ridge_k3
  loocv_q0 <- mean(((y - yduga_q) / (1 - diag(Hq)))^2)
  loocv_k20 <- mean(((y - yduga_k2) / (1 - diag(Hk2)))^2)
  loocv_k30 <- mean(((y - yduga_k3) / (1 - diag(Hk3)))^2)
  loocv_q <- c(loocv_q, loocv_q0)
  loocv_k2 <- c(loocv_k2, loocv_k20)
  loocv_k3 <- c(loocv_k3, loocv_k30)
  print(i)
}
#queen
which.min(loocv_q) #13
lambdas[13] #[1] 0.1734694
loocv_q[13] #[1] 2.904429
#k-nn (k=2)
which.min(loocv_k2) #40
lambdas[40] #[1] 0.3387755
loocv_k2[40] #[1] 0.3757478
#k-nn (k=3)
which.min(loocv_k3) #32
lambdas[32] #[1] 0.2897959
loocv_k3[32] #[1] 0.4423482

###RMSE
beta_rq <- ginv(t(X)%*%X+lambdas[13]*t(D)%*%D)%*%t(X)%*%y
beta_rk2 <- ginv(t(X)%*%X+lambdas[40]*t(D.k2)%*%D.k2)%*%t(X)%*%y
beta_rk3 <- ginv(t(X)%*%X+lambdas[32]*t(D.k3)%*%D.k3)%*%t(X)%*%y
predr.q <- X%*%beta_rq
predr.k2 <- X%*%beta_rk2
predr.k3 <- X%*%beta_rk3

rmser.q <- sqrt(mean((y-predr.q)^2)) #[1] 0.06270314
rmser.k2 <- sqrt(mean((y-predr.k2)^2)) #[1] 0.1252489
rmser.k3 <- sqrt(mean((y-predr.k3)^2)) #[1] 0.1485008

Hat.q <- X%*%ginv(t(X)%*%X+lambdas[13]*t(D)%*%D)%*%t(X)
Hat.k2 <- X%*%ginv(t(X)%*%X+lambdas[40]*t(D.k2)%*%D.k2)%*%t(X)
Hat.k3 <- X%*%ginv(t(X)%*%X+lambdas[32]*t(D.k3)%*%D.k3)%*%t(X)

df.q <- sum(diag(Hat.q)) #[1] 156.8579
df.k2 <- sum(diag(Hat.k2)) #[1] 141.2237
df.k3 <- sum(diag(Hat.k3)) #[1] 136.9212

###Generalized Lasso Modeling
####function for ALOCV and GCV
ALOG <- function(genlasso_model, X, D, y, lam.range = NULL, epsilon = 1e-10, lim1 = 1-1e-3, iter=100, seqq=TRUE) {
  cat("---------------------------\nEVALUATING:\n")
  cat("Range lambda of the model:", min(genlasso_model$lambda), "to", max(genlasso_model$lambda), "\n")
  cat("Range df of the model:", min(genlasso_model$df), "to", max(genlasso_model$df), "\n")
  cat("---------------------------\nCONSTRAINT:\n")
  if (!is.null(lam.range)) {
    lam <- genlasso_model$lambda[genlasso_model$lambda >= lam.range[1] & genlasso_model$lambda <= lam.range[2]] %>% rev()  
    if(seqq){
      lam <- exp(seq(log(lam.range[1]), log(lam.range[2]), length=iter))
    }
    cat("lam.range:", min(lam),"to", max(lam), "\n")
  } else {
    lam <- rev(genlasso_model$lambda)
  }
  
  cat("---------------------------\n")
  cat("iterations =",length(lam),"\n")
  cat("l(E) : Length of E (satisfy the condition |u| ≈ lambda):\n")
  
  alpha <- coef(genlasso_model, lambda = lam, type = "both")
  u <- alpha$u
  df <- alpha$df
  alot <- c()
  gt <- c()
  
  for(j in 1:length(lam)) {
    E <- which(abs(abs(u[,j]) - lam[j]) < epsilon)
    cat("j =", j, "| l(E) =", length(E), "| df =", df[j],"| lambda =", lam[j], "\n")
    if (length(E) == 0) {Dbaru <- D} else {Dbaru <- D[-E,]}
    if (length(E) == nrow(D)-1) Dbaru <- matrix(Dbaru, 1, ncol(D))
    A <- X %*% pracma::nullspace(Dbaru)
    Aplus <- MASS::ginv(A)
    H <- A %*% Aplus
    hj <- diag(H)
    trH <- mean(hj)
    hj <- ifelse(hj >= lim1, lim1, hj)
    pred <- predict(genlasso_model, lambda = lam[j], Xnew = X)
    alo <- mean(((y - pred$fit) / (1 - hj))^2)
    alot <- c(alot, alo)
    go <- mean(((y - pred$fit) / (1 - trH))^2)
    gt <- c(gt, go)
  }
  
  best.lam <- lam[which.min(alot)]
  best.df <- df[which.min(alot)]
  
  alo_msg <- sprintf("Best lambda for ALOCV: %f, with DF: %d, at index: %d", best.lam, best.df, which.min(alot))
  cat(alo_msg,"\n")
  
  best.lamg <- lam[which.min(gt)]
  best.dfg <- df[which.min(gt)]
  gcv_msg <- sprintf("Best lambda for GCV: %f, with DF: %d, at index: %d", best.lamg, best.dfg, which.min(gt))
  cat(gcv_msg,"\n")
  
  return(list(ALO = alot, 
              GCV = gt,
              lam = lam,
              bestlam_ALO = best.lam,
              bestlam_GCV = best.lamg,
              results = c(alo_msg, gcv_msg)))
}

####modeling
genlasso_model <- genlasso(y, X, D) #genlasso queen
eval1 <- ALOG(genlasso_model, X, D, y, 
              lam.range = c(0.1, 1.7),
              epsilon = 1e-10,
              lim1 = 1-1e-1,
              iter=100,
              seqq=T
)
#Evaluate when lambda <=4 (for reasonable interpretation)
which.min(eval1$ALO[eval1$lam<=0.4]) #46
eval1$lam[46] #[1] 0.3624891
eval1$ALO[46] #[1] 0.5370916
which.min(eval1$GCV[eval1$lam<=0.4]) #47
eval1$lam[47] #[1] 0.3730128
eval1$GCV[47] #[1] 0.4075174

genlasso_model_knn2 <- genlasso(y,X,D.k2) #genlasso k-nn (k=2)
eval2 <- ALOG(genlasso_model_knn2, X, D.k2, y,
              lam.range = c(0.1, 4),
              epsilon = 1e-10,
              lim1 = 1-1e-1,
              iter=100,
              seqq=T
)
#Evaluate when lambda <=4 (for reasonable interpretation)
which.min(eval2$ALO[eval2$lam<=0.4]) #31
eval2$lam[31] #[1] 0.3058248
eval2$ALO[31] #[1] 1.888772
which.min(eval2$GCV[eval2$lam<=0.4]) #33
eval2$lam[33] #[1] 0.3294865
eval2$GCV[33] #[1] 0.4688366

genlasso_model_knn3 <- genlasso(y,X,D.k3) #genlasso k-nn (k=3)
eval3 <- ALOG(genlasso_model_knn3, X, D.k3, y,
              lam.range = c(0.1, 2),
              epsilon = 1e-10,
              lim1 = 1-1e-1,
              iter=100,
              seqq=T
)
#Evaluate when lambda <=4 (for reasonable interpretation)
which.min(eval3$ALO[eval3$lam<=0.4]) #7
eval3$lam[7] #[1] 0.1199086
eval3$ALO[7] #[1] 3.61982
which.min(eval3$GCV[eval3$lam<=0.4]) #29
eval3$lam[29] #[1] 0.2333287
eval3$GCV[29] #[1] 0.7728948

###RMSE
predq.a <- predict(genlasso_model, lambda = eval1$lam[46], Xnew = X)
predq.g <- predict(genlasso_model, lambda = eval1$lam[47], Xnew = X)
predk2.a <- predict(genlasso_model_knn2, lambda = eval2$lam[31], Xnew = X)
predk2.g <- predict(genlasso_model_knn2, lambda = eval2$lam[33], Xnew = X)
predk3.a <- predict(genlasso_model_knn3, lambda = eval3$lam[7], Xnew = X)
predk3.g <- predict(genlasso_model_knn3, lambda = eval3$lam[29], Xnew = X)

rmseq.a <- sqrt(mean((y-predq.a$fit)^2)) #[1] 0.3330734
rmseq.g <- sqrt(mean((y-predq.g$fit)^2)) #[1] 0.334206
rmsek2.a <- sqrt(mean((y-predk2.a$fit)^2)) #[1] 0.407238
rmsek2.g <- sqrt(mean((y-predk2.g$fit)^2)) #[1] 0.4229131
rmsek3.a <- sqrt(mean((y-predk3.a$fit)^2)) #[1] 0.2744915
rmsek3.g <- sqrt(mean((y-predk3.g$fit)^2)) #[1] 0.4447436

###genlasso result
beta_hats_result3 <- coef(genlasso_model_knn3, 
                          lambda = eval3$lam[7])$beta #knn3 ALOCV

###reformat result
beta_all3 <- matrix(0,42,34)
for(i in 1:42){
  for(j in 1:34){
    #beta_all[i,j] <- beta_hats_result[j+34*(i-1)]
    beta_all3[i,j] <- beta_hats_result3[j+34*(i-1)]
  }
}
result_x1 <- as.matrix(beta_all3[1:5,]) #+
result_x2 <- as.matrix(beta_all3[6:10,]) #+
result_x3 <- as.matrix(beta_all3[11:13,]) #+
result_x4 <- as.matrix(beta_all3[14:18,]) #-
result_x5 <- as.matrix(beta_all3[19:22,]) #-
result_x6 <- as.matrix(beta_all3[23:27,]) #+
result_x7 <- as.matrix(beta_all3[28:32,]) #-
result_x8 <- as.matrix(beta_all3[33:37,]) #-
result_x9 <- as.matrix(beta_all3[38:42,]) #-

english_names <- c("Aceh", "Bali", "Banten", 
                   "Bengkulu", "Yogyakarta", "Jakarta", 
                   "Gorontalo", "Jambi", "West Java", 
                   "Central Java", "East Java", "West Kalimantan", 
                   "South Kalimantan", "Central Kalimantan", 
                   "East Kalimantan", "North Kalimantan", 
                   "Bangka Belitung Islands", "Riau Islands", 
                   "Lampung", "Maluku", "North Maluku", 
                   "West Nusa Tenggara", "East Nusa Tenggara", 
                   "Papua", "West Papua", "Riau", "West Sulawesi", 
                   "South Sulawesi", "Central Sulawesi", 
                   "Southeast Sulawesi", "North Sulawesi", 
                   "West Sumatra", "South Sumatra", 
                   "North Sumatra")
year.name <- c(2019,2020,2021,2022,2024)
year.name2 <- c(2020,2021,2022)
year.name3 <- c(2019,2020,2021,2022)

result_x1ok <- ifelse(result_x1<0,0,result_x1); colnames(result_x1ok) <-english_names; rownames(result_x1ok) <- year.name 
result_x2ok <- ifelse(result_x2<0,0,result_x2); colnames(result_x2ok) <-english_names; rownames(result_x2ok) <- year.name 
result_x3ok <- ifelse(result_x3<0,0,result_x3); colnames(result_x3ok) <-english_names; rownames(result_x3ok) <- year.name2 
result_x4ok <- ifelse(result_x4>0,0,result_x4); colnames(result_x4ok) <-english_names; rownames(result_x4ok) <- year.name 
result_x5ok <- ifelse(result_x5>0,0,result_x5); colnames(result_x5ok) <-english_names; rownames(result_x5ok) <- year.name3 
result_x6ok <- ifelse(result_x6<0,0,result_x6); colnames(result_x6ok) <-english_names; rownames(result_x6ok) <- year.name 
result_x7ok <- ifelse(result_x7>0,0,result_x7); colnames(result_x7ok) <-english_names; rownames(result_x7ok) <- year.name 
result_x8ok <- ifelse(result_x8>0,0,result_x8); colnames(result_x8ok) <-english_names; rownames(result_x8ok) <- year.name 
result_x9ok <- ifelse(result_x9>0,0,result_x9); colnames(result_x9ok) <-english_names; rownames(result_x9ok) <- year.name 

##Figure 5. Heatmap results
#X1 Poverty
heatmaply(result_x1ok, 
          scale_fill_gradient_fun = scale_fill_gradient2(low="red", high="blue"), 
          plot_method="ggplot", 
          Colv = TRUE, 
          Rowv = FALSE, 
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "complete"),
          fontsize_row = 9, 
          fontsize_col = 9, 
          column_text_angle=75)
#X2 EBF
heatmaply(result_x2ok, 
          scale_fill_gradient_fun = scale_fill_gradient2(low="red", high="blue"), 
          plot_method="ggplot", 
          Colv = TRUE, 
          Rowv = FALSE, 
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "complete"),
          fontsize_row = 9, 
          fontsize_col = 9, 
          column_text_angle=75)
#X3 LBW
heatmaply(result_x3ok, 
          scale_fill_gradient_fun = scale_fill_gradient2(low="red", high="blue"), 
          plot_method="ggplot", 
          Colv = TRUE, 
          Rowv = FALSE, 
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "complete"),
          fontsize_row = 9, 
          fontsize_col = 9, 
          column_text_angle=75)
#X4 High School Completion
heatmaply(result_x4ok, 
          scale_fill_gradient_fun = scale_fill_gradient2(low="red", high="blue"), 
          plot_method="ggplot", 
          Colv = TRUE, 
          Rowv = FALSE, 
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "complete"),
          fontsize_row = 9, 
          fontsize_col = 9, 
          column_text_angle=75)
#X5 Proper Sanitation Access
heatmaply(result_x5ok, 
          scale_fill_gradient_fun = scale_fill_gradient2(low="red", high="blue"), 
          plot_method="ggplot", 
          Colv = TRUE, 
          Rowv = FALSE, 
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "complete"),
          fontsize_row = 9, 
          fontsize_col = 9, 
          column_text_angle=75)
#X6 Unmet need
heatmaply(result_x6ok, 
          scale_fill_gradient_fun = scale_fill_gradient2(low="red", high="blue"), 
          plot_method="ggplot", 
          Colv = TRUE, 
          Rowv = FALSE, 
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "complete"),
          fontsize_row = 9, 
          fontsize_col = 9, 
          column_text_angle=75)
#X7 GDP
heatmaply(result_x7ok, 
          scale_fill_gradient_fun = scale_fill_gradient2(low="red", high="blue"), 
          plot_method="ggplot", 
          Colv = TRUE, 
          Rowv = FALSE, 
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "complete"),
          fontsize_row = 9, 
          fontsize_col = 9, 
          column_text_angle=75)
#X8 Calory consumption
heatmaply(result_x8ok, 
          scale_fill_gradient_fun = scale_fill_gradient2(low="red", high="blue"), 
          plot_method="ggplot", 
          Colv = TRUE, 
          Rowv = FALSE, 
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "complete"),
          fontsize_row = 9, 
          fontsize_col = 9, 
          column_text_angle=75)
#X9 Protein Consumption
heatmaply(result_x9ok, 
          scale_fill_gradient_fun = scale_fill_gradient2(low="red", high="blue"), 
          plot_method="ggplot", 
          Colv = TRUE, 
          Rowv = FALSE, 
          distfun = function(x) dist(x, method = "euclidean"),
          hclustfun = function(x) hclust(x, method = "complete"),
          fontsize_row = 9, 
          fontsize_col = 9, 
          column_text_angle=75)




