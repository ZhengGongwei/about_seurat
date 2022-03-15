Seurat2的数据可以用下面的函数转化
PBMC = UpdateSeuratObject(object = PBMC)

DimPlot(PBMC, label = T)+NoLegend()+labs(x = "UMAP1", y = "UMAP2") +
theme(axis.text.y = element_blank(),
axis.ticks.y = element_blank(),
axis.text.x = element_blank(),
axis.ticks.x = element_blank())+
theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))  #加边框

library(paletteer) #提供了 R 编程环境中提供的数百种其他调色板的组合集合，详情可以查看此包，总有你满意的方案
pal <- paletteer_d("ggsci::nrc_npg")[c(1,3,4,9,5,2,6,8,10)] #有几群细胞需要标记就选几种颜色
DimPlot(PBMC, label = T,cols= pal,pt.size = 1.5,repel = T)+NoLegend()

library(RColorBrewer)#配置自己需要的颜色  
cell_type_cols <- c(brewer.pal(9, "Set1"), "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F","#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B","#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00","#000000","#708090")  
DimPlot(PBMC, label = T,cols= cell_type_cols,pt.size = 1.5,repel = T)+NoLegend()

FeaturePlot(PBMC, features = "S100A8")
mycolor <- c(''lightgrey'', ''blue'',''seagreen2'')#设置颜色
FeaturePlot(PBMC, features = ''S100A8'',cols = mycolor, pt.size = 1.5)+theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#加边框 


install.packages("remotes")
remotes::install_github("lyc-1995/MySeuratWrappers")
library(MySeuratWrappers)
markers <- c(''CD3D'', ''S100A8'', ''S100A9'', ''CD79A'', ''CCL5'', ''NKG7'', ''GZMA'', ''IL32'', ''CD4'', ''CD8A'', ''LTB'', ''FCN1'', ''MS4A1'', ''SPON2'',''FCER1A'',''SERPINF1'', ''TMEM40'', ''CD3E'')
my36colors <-c(''#E5D2DD'', ''#53A85F'', ''#F1BB72'', ''#F3B1A0'', ''#D6E7A3'', ''#57C3F3'', ''#476D87'', ''#E95C59'', ''#E59CC4'', ''#AB3282'', ''#23452F'', ''#BD956A'', ''#8C549C'', ''#585658'', ''#9FA3A8'', ''#E0D4CA'', ''#5F3D69'', ''#C5DEBA'', ''#58A4C3'', ''#E4C755'', ''#F7F398'', ''#AA9A59'', ''#E63863'', ''#E39A35'', ''#C1E6F3'', ''#6778AE'', ''#91D0BE'', ''#B53E2B'', ''#712820'', ''#DCC1DD'', ''#CCE0F5'',  ''#CCC9E6'', ''#625D9E'', ''#68A180'', ''#3A6963'', ''#968175'')
#颜色设置
VlnPlot(PBMC, features = markers,stacked=T,pt.size=0.9,cols = my36colors,direction = "horizontal",x.lab = '''', y.lab = '''')+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())#不显示坐标刻度 

DotPlot(PBMC, features = markers)+coord_flip()+theme_bw()+theme(panel.grid = element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+scale_color_gradientn(values = seq(0,1,0.2),ccolours = c(''#330066'',''#336699'',''#66CC66'',''#FFCC33''))+labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3)) 
