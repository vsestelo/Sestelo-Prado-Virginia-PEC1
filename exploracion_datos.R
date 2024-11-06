## Código para la exploración de los datos

# Matriz de los datos
matriz <- assays(se)[[1]]  

# Ver las primeras filas y columnas de la matriz
head(matriz)

# Resumen estadístico de cada tratamiento
summary(matriz)

# Metadatos de las columnas (información de las muestras)
colData(se)

# Metadatos de las filas (información de los atributos)
rowData(se)

# Metadatos del experimento
metadata(se)

# Dimensiones del objeto SummarizedExperiment
dim(se)

# Nombres de las columnas
colnames(se)

# Nombres de las primeras 6 filas
head(rownames(se))

# Histogramas de las tres primeras muestras
opt <- par(mfrow=c(1,3))
for (i in 1:3)
  hist(matriz[,i], main = colnames(matriz)[i])
par(opt)

# Histograma del tercer metabolito
subse <- se[3,]
msubse <- assays(subse)[[1]] 

hist(msubse, main="Distribución del nivel de 10-Oxodecanoate_1",
     xlab="Nivel de 10-Oxodecanoate_1",
     ylab="Frecuencia",
     breaks=100)

colData(se)[1,]
colData(se)[4,]

# Boxplot de los cinco primeros metabolitos en todas las muestras
subse2 <- se[1:5,]
msubse2 <- assays(subse2)[[1]] 
groupColors <- c(rep("red", 15), rep("blue", 15),rep("green", 15) )
boxplot(msubse2, col=groupColors, main="Niveles de señal de los 5 primeros metabolitos en cada tratamiento",
        xlab="Tratamientos",
        ylab="Nivel de señal", las=2, cex.axis=0.7, cex.main=0.7)

# Boxplot logarítmico
groupColors <- c(rep("red", 15), rep("blue", 15),rep("green", 15) )
boxplot(log2(msubse2), col=groupColors, main="Niveles de señal de los 5 primeros metabolitos en cada tratamiento",
        xlab="Tratamientos",
        ylab="Nivel de señal", las=2, cex.axis=0.7, cex.main=0.7)

# PCA
# Se transforma logarítmicamente la matriz (se suma 1 para evitar log(0))
logX <- log2(matriz+1)

# Se elimina cualquier fila/columna con NA tras la transformación
logX <- na.omit(logX)

# Se realiza el PCA en los datos transpuestos y escalados
pcX <- prcomp(t(logX), scale = TRUE)

# Se calcula la varianza explicada de cada componente
loads <- round(pcX$sdev^2 / sum(pcX$sdev^2) * 100, 1)

# Se representan los resultados de los dos primeros componentes
xlab<-c(paste("PC1",loads[1],"%"))
ylab<-c(paste("PC2",loads[2],"%"))
plot(pcX$x[,1:2],xlab=xlab,ylab=ylab, col=groupColors, 
     main ="Principal components (PCA)")
names2plot<-paste0(substr(colnames(matriz),1,3), 1:4)

text(pcX$x[,1],pcX$x[,2],names2plot, pos=3, cex=.6)

# Clúster jerárquico
clust.euclid.average <- hclust(dist(t(matriz)),method="average")
plot(clust.euclid.average, hang=-1)