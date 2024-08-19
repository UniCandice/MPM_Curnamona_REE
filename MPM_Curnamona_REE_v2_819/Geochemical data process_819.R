# Geochemical data processing (Including sampling point conversion and robust principal component analysis)
# Set random seed
set.seed(1)
# Load the required package
library(compositions)
library(factoextra)
library(rrcov)
library(corrplot)
library(sf)          
library(dplyr)       
library(readr)       
library(tidyr)       
library(stats)       
library(progress)  

# Read coordinates of target sampling point
tar_sam <- st_read('C:/Users/14361/AU_MPM/clip_dataset/output/DevNet_MPM/random_point_unlabel_train/geochemical_unlabel/train_random/geochemical_unlabel_all_869_XY.shp')
tar_sam_coords <- st_coordinates(tar_sam)

# Read the original geochemical data
geochem_raw_path <- 'D:/Dataset/ILRRPCA.csv'
geochem_raw <- read_csv(geochem_raw_path)
element_names <- colnames(geochem_raw)[3:ncol(geochem_raw)]

# Sampling points feature conversion
# Calculate feature concentrations for each target sampling point
tar_sam_feature <- list()
pb <- progress_bar$new(total = nrow(tar_sam_coords), format = "  [:bar] :percent in :elapsed")

# Iteratively compute concentration values
for (i in 1:nrow(tar_sam_coords)) {
  concentration_values <- c(tar_sam_coords[i, "X"], tar_sam_coords[i, "Y"])
  for (col in 3:ncol(geochem_raw)) {
    median_value <- median(geochem_raw[[col]], na.rm = TRUE)
    oral_geoc_coords <- geochem_raw[, 1:2]
    geoc_dist <- sqrt((tar_sam_coords[i, "X"] - oral_geoc_coords[[1]])^2 + 
                          (tar_sam_coords[i, "Y"] - oral_geoc_coords[[2]])^2)
    geoc_ids <- which(geoc_dist < 0.015)
    
    if (length(geoc_ids) == 0) {
      concentration_values <- c(concentration_values, median_value)
    } else {
      numerator <- sum(geochem_raw[geoc_ids, col] / (geoc_dist[geoc_ids]^2), na.rm = TRUE)
      denominator <- sum(1 / (geoc_dist[geoc_ids]^2), na.rm = TRUE)
      concentration_values <- c(concentration_values, numerator / denominator)
    }
  }
  
  tar_sam_feature[[i]] <- concentration_values
  pb$tick()
}

# Save the converted data
tar_sam_feature_df <- as.data.frame(do.call(rbind, tar_sam_feature))
colnames(tar_sam_feature_df) <- c("X", "Y", element_names)
write_csv(tar_sam_feature_df, 'C:/Users/14361/AU_MPM/clip_dataset/output/DevNet_MPM/random_point_unlabel_train/geochemical_unlabel/train_random/geochemical_unlabel_random_0p015_all_elements_R_814.csv')

# Robust principal component analysis
# Read the converted data
f=read.csv("C:/Users/14361/AU_MPM/clip_dataset/output/DevNet_MPM/random_point_unlabel_train/geochemical_unlabel/train_random/geochemical_unlabel_random_0p015_all_elements_R_814.csv",header = TRUE)  
sel=c(3:29) 
comp_data=f[,sel]

# Isometric log ratio transform
ilr_data <- ilr(comp_data)
summary(ilr_data)

#Perform robust principal component analysis
ilr_data.mcd=covMcd(ilr_data,cor=TRUE)
summary(ilr_data.mcd)
pca_result=princomp(ilr_data,covmat=ilr_data.mcd,cor=TRUE)

# Extract eigenvalues and visualization results
get_eigenvalue(pca_result)
fviz_eig(pca_result,addlabels = T,ylim=c(0,50))
pca_result.var <- get_pca_var(pca_result)
pca_result.var$cor
pca_result.var$coord  
pca_result.var$cos2
pca_result.var$contrib
fviz_pca_var(pca_result,col.var = "contrib",gradient.cols=c("#00AFBB","#E7B800","#FC4E07"),axes=1:2)
corrplot(pca_result.var$coord,is.corr=F,tl.cex = 0.75)
corrplot(pca_result.var$cos2,is.corr=F,tl.cex = 0.75)
corrplot(pca_result.var$contrib,is.corr=F,tl.cex = 0.75)
fviz_contrib(pca_result, choice = "var", axes = 1,tl.cex = 0.75)
fviz_contrib(pca_result, choice = "var", axes = 1:2,tl.cex = 0.75)
fviz_cos2(pca_result, choice = "var",axes = 1,tl.cex = 0.75)

# Extract the loading matrix and scores of principal components
loadings <- pca_result$loadings
scores <- pca_result$scores
write.csv(pca_result$loadings,"C:/Users/14361/AU_MPM/clip_dataset/output/RPCA/ilr_backilr/random_train/MPM_ilr_RPCA_geochemical_loadings814_0p015_random_11.csv") 
write.csv(pca_result$scores,"C:/Users/14361/AU_MPM/clip_dataset/output/RPCA/ilr_backilr/random_train/MPM_ilr_RPCA_geochemical_Scores814_0p015_random_11.csv") 

# Ready to invert to restore to original space
pca_resultback=pca_result
# construct the orthonormal basis
V=matrix(0,nrow=27,ncol=26) 
for (i in 1:ncol(V)){
  V[1:i,i] = 1/i
  V[i+1,i] = (-1)
  V[,i] = V[,i]*sqrt(i/(i+1))
}
pca_resultback$loadings = V%*%pca_result$loadings
pca_resultback$scores = pca_result$scores%*%t(V)
dimnames(pca_resultback$loadings)[[1]] = names(comp_data)
write.csv(pca_resultback$loadings,"C:/Users/14361/AU_MPM/clip_dataset/output/RPCA/ilr_backilr/random_train/MPM_backilr_RPCA_geochemical_loadings814_0p015_random_11.csv") 
write.csv(pca_resultback$scores,"C:/Users/14361/AU_MPM/clip_dataset/output/RPCA/ilr_backilr/random_train/MPM_backilr_RPCA_geochemical_Scores814_0p015_random_11.csv") 

# Calculate the eigenvalues and interpretive variance ratio after inverse transformation
eigenvalues_back = (pca_resultback$sdev)^2
variance_explained_back = eigenvalues_back / sum(eigenvalues_back) * 100
barplot(eigenvalues_back, main = "Eigenvalues (Transformed Space)", 
        xlab = "Principal Component", ylab = "Eigenvalue",
        names.arg = paste0("PC", 1:length(eigenvalues_back)), 
        col = "steelblue", ylim = c(0, max(eigenvalues_back) * 1.1))
plot(variance_explained_back, type = "b", main = "Variance Explained (Transformed Space)",
     xlab = "Principal Component", ylab = "Variance Explained (%)",
     xaxt = "n", ylim = c(0, max(variance_explained_back) * 1.1), 
     col = "darkred", pch = 19, cex = 1.5)
axis(1, at = 1:length(variance_explained_back), labels = paste0("PC", 1:length(variance_explained_back)))
text(1:length(variance_explained_back), variance_explained_back, 
     labels = paste0(round(variance_explained_back, 2), "%"), 
     pos = 3, cex = 0.5)

# Select the first 15 principal components for detailed analysis
n_components = 15
eigenvalues_back_selected = eigenvalues_back[1:n_components]
variance_explained_back_selected = variance_explained_back[1:n_components]
barplot(eigenvalues_back_selected, main = "Eigenvalues (Transformed Space, First 15 PCs)", 
        xlab = "Principal Component", ylab = "Eigenvalue",
        names.arg = paste0("PC", 1:n_components), 
        col = "steelblue", ylim = c(0, max(eigenvalues_back_selected) * 1.1))
plot(variance_explained_back_selected, type = "b", main = "Variance Explained (Transformed Space, First 15 PCs)",
     xlab = "Principal Component", ylab = "Variance Explained (%)",
     xaxt = "n", ylim = c(0, max(variance_explained_back_selected) * 1.1), 
     col = "darkred", pch = 19, cex = 1.5)
axis(1, at = 1:n_components, labels = paste0("PC", 1:n_components))
text(1:n_components, variance_explained_back_selected, 
     labels = paste0(round(variance_explained_back_selected, 2), "%"), 
     pos = 3, cex = 0.6)

# Calculate converted coordinates, cos², and contributions
pca_resultback.var <- get_pca_var(pca_resultback)
coord_back <- pca_resultback.var$coord
cos_back <- pca_resultback.var$cos2
contrib_back <- pca_resultback.var$contrib
corrplot(coord_back, is.corr = FALSE, method = "shade", 
         main = "Correlation Plot (Transformed Coordinates)", 
         mar = c(0, 0, 2, 0), tl.cex = 0.7)
corrplot(cos_back, is.corr = FALSE, method = "square", 
         main = "Cos2 Plot (Transformed Cosines)", 
         mar = c(0, 0, 2, 0), tl.cex = 0.7)
corrplot(contrib_back, is.corr = FALSE, method = "shade", 
         main = "Correlation Plot (Transformed Contributions)", 
         mar = c(0, 0, 2, 0), tl.cex = 0.7)


