install.packages("forestplot")# 安装包
library(forestplot)# 运行包

# 读取整理好的数据
rs_forest <- read.csv('森林图R格式数据.csv', header = FALSE)

# 获取权重信息
#weights <- rs_forest$V9
normalized_weights<- rs_forest$V9
# 排除 NA 值，只处理有效权重
valid_weights <- weights[!is.na(weights)]

# 标准化权重，以便用于 boxsize 参数
max_weight <- max(valid_weights)
normalized_weights <- weights / max_weight

# 创建包含样本量信息的标签文本
labels <- as.matrix(rs_forest[,1:5])

# 生成森林图
forestplot(
  labeltext = labels, # 设置前四列是文本
  mean = rs_forest$V6,
  lower = rs_forest$V7,
  upper = rs_forest$V8, # 分别设置点估计值和上下置信区间
  #敏感后：T,T,F,F,F,F,F,F,T,T,F,F,F,F,F,F,F,T,T,F,F,F,T,T,F,F,F,T,T,F,F,F,F,F,T
  #敏感前：T,T,F,F,F,F,F,F,F,T,T,F,F,F,F,F,F,F,T,T,F,F,F,T,T,F,F,F,F,T,T,F,F,F,F,F,F,T
  is.summary = c(T,T,F,F,F,F,F,F,F,T,T,F,F,F,F,F,F,F,T,T,F,F,F,T,T,F,F,F,F,T,T,F,F,F,F,F,F,T),
  # 设置加粗
  zero = 0, # 设置无效线
  boxsize = normalized_weights, # 设置点估计的方形大小，根据权重调整
  lineheight = unit(250, 'mm'), # 设置图形行距
  colgap = unit(2, 'mm'), # 设置列间距
  lwd.zero = 1.5, # 设置参考线的粗细
  lwd.ci = 2, # 设置区间估计线的粗细
  col = fpColors(box = '#21a675',summary = '#614099',lines = 'black',zero = '#8ea5c8'),
  xlab = "SMD",
  lwd.xaxis = 2, # 设置X轴的粗细
  lty.ci = 'solid',
  graph.pos = 3 # 设置森林图出现在第几列
)
